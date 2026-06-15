# Timekeepers.jl

Time-series IO and interactive inspection for magnetotelluric field data.

Timekeepers reads logger-native files into `TimeSeries.jl` `TimeArray`s, keeps
native writer paths for supported instruments, and ships **TKApp**, a native
GLMakie window for fast visual inspection and masking of long records.

![TKApp time series view](images/ts.png)

A five-channel LEMI record after a few intervals were masked in TKApp, written
out, and reloaded — masked windows render as gaps in the traces.

Supported native formats:

- LEMI-424 long-period ASCII text files (24-column)
- GEOMAG-02 ASCII text files with semicolon headers
- Metronix ADU (ATS binary + XML) site directories
- Generic LEMI-style `.xyz` exports (7-column `date time Bx By Bz Ex Ey`)

## Installation

```julia
pkg> activate .
pkg> instantiate
```

TKApp uses GLMakie, so you need a desktop session with OpenGL 3.3+ drivers.
Quick smoke test — if this opens a window, `run_tkapp()` will too:

```julia
using GLMakie
display(scatter(1:10))
```

Headless SSH, WSL, containers, and very old GPUs may need extra OpenGL setup.

## Getting test data

Sample data is **not** shipped with the repository; the `data/` directory is
gitignored so you can drop large recordings in without polluting the repo. For
a quick test, download a public LEMI-424 dataset from the British Geological
Survey accession and extract the `.txt` into `data/`:

> https://webapps.bgs.ac.uk/services/ngdc/accessions/index.html#item182849

## Quick Start

Open the explorer window:

```julia
using Timekeepers
run_tkapp()
```

or run the bundled launcher script:

```powershell
julia --project=. examples/tkapp.jl
```

The toolbar covers the workflow: **Load** a file (or site directory), pick a
**Window** length and scroll through the record, **left-drag** to select an
interval, **Mask** / **Unmask** / **Clear** to edit it, and **Write** to
export. Masked rows are written as `NaN` in the original format alongside a
mask file of the intervals. The `View:` menu adds per-channel PSD and
spectrogram panels.

The same workflow from Julia code:

```julia
using Timekeepers, Dates

# Load
ta = first(Timekeepers._load_data_file("data/LEMI090.txt"))   # → TimeArray

# Mask
mask = TimekeeperMask(ta)
mask_interval!(mask, DateTime(2020, 10, 4, 0, 10), DateTime(2020, 10, 4, 0, 20))

# Derive processing-ready outputs
cleaned  = cleaned_timearray(ta, mask)                # NaN in masked rows
segments = good_segments(ta, mask; min_samples = 256) # contiguous good chunks
weights  = sample_weights(mask)                       # for robust processing

# Write (same format as loaded)
write_lemi424("data/LEMI090_clean.txt", cleaned)
write_mask("data/LEMI090_mask.csv", mask)
```

## Fast startup (optional sysimage)

A fresh `julia ... run_tkapp()` spends tens of seconds loading GLMakie and
compiling plotting code before the window appears. For day-to-day use — or for
handing the tool to a colleague — build a **PackageCompiler sysimage** that
bakes Timekeepers, GLMakie, and a warm-up render into one image, dropping
startup to a second or two.

1. Install PackageCompiler in the **global** environment, so it stays out of
   this package's `Project.toml`. Do **not** pass `--project=.` here — that
   would add it as a Timekeepers dependency and force every user to install it:

   ```powershell
   julia -e 'using Pkg; Pkg.activate(); Pkg.add("PackageCompiler")'
   ```

   (If you accidentally added it to the project, remove it with
   `julia --project=. -e 'using Pkg; Pkg.rm("PackageCompiler")'`.)

2. Build the image from the repo root. The output extension is platform
   specific — `timekeepers.dll` on Windows, `timekeepers.so` on Linux,
   `timekeepers.dylib` on macOS. `scripts/warmup.jl` opens and renders the app
   once so the GL/plotting paths are traced into the image:

   ```powershell
   julia -e 'using PackageCompiler; create_sysimage(["Timekeepers"]; sysimage_path = "timekeepers.dll", project = ".", precompile_execution_file = "scripts/warmup.jl")'
   ```

3. Launch with the image — inline or via the bundled script (flags first,
   script last):

   ```powershell
   julia --project=. --sysimage timekeepers.dll -e 'using Timekeepers; run_tkapp()'
   julia --project=. --sysimage timekeepers.dll examples/tkapp.jl
   ```

The sysimage is specific to the OS, CPU, and Julia version it was built on, and
is pinned to the package code at build time — so it is **gitignored**, must be
**rebuilt after changes to `src/`**, and a colleague generally rebuilds it on
their own machine (steps 1–2). `scripts/warmup.jl` travels with the repo; the
image does not.
