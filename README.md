# Timekeepers.jl

Time-series IO and interactive inspection for magnetotelluric field data.

Timekeepers reads logger-native files into `TimeSeries.jl` `TimeArray`s, keeps
native writer paths for supported instruments, and ships **TKApp**, a native
GLMakie window for fast visual inspection and masking of long records.

![TKApp time series view](images/ts.png)

A five-channel LEMI record after a few intervals were marked masked in TKApp,
written out, and reloaded — masked windows render as gaps in the traces.

Supported native formats:

- LEMI-424 long-period ASCII text files (24-column)
- GEOMAG-02 ASCII text files with semicolon headers
- Generic LEMI-style `.xyz` exports (7-column `date time Bx By Bz Ex Ey`)

## Installation

```julia
pkg> activate .
pkg> instantiate
```

### Linux graphics support

TKApp uses GLMakie, so a standard Linux laptop usually works if it has a
normal desktop session and working OpenGL 3.3+ graphics drivers. Intel and AMD
integrated graphics generally work with Mesa drivers; NVIDIA systems generally
work when the NVIDIA driver is installed correctly.

Headless SSH sessions, WSL setups, containers, and very old laptops may need
extra OpenGL/display setup. A quick smoke test is:

```julia
using GLMakie
display(scatter(1:10))
```

If this opens a window, `run_tkapp()` should open too. A 10 Hz one-day GEOMAG
file such as `GEOMAG.TXT` is about 864,000 samples; the data itself fits
comfortably in 16 GB RAM, though spectra and spectrogram views still depend on
CPU/GPU speed.

## Getting test data

Sample data is **not** shipped with the repository. The `data/` directory is
gitignored so you can drop large recordings in without polluting the repo.

For a quick test you can download a public LEMI-424 dataset from the British
Geological Survey accession:

> https://webapps.bgs.ac.uk/services/ngdc/accessions/index.html#item182849

Extract the `.txt` file(s) into the `data/` directory at the repository root:

```
Timekeepers.jl/
├── data/
│   └── LEMI090.txt        # your downloaded file
├── examples/
└── src/
```

## Quick Start

Open the native explorer window:

```julia
using Timekeepers
run_tkapp()
```

Then in the toolbar:

1. **Load…** — pick a `.txt` (LEMI-424) or `.xyz` file from `data/`.
2. **Window** — choose visible window length (1 min … 7 days, or All).
3. **Scroll** slider — pan the visible window through the full record.
4. **Left-drag** on any panel — paints a translucent blue selection band
   across all five channels for the same time range.
5. **Mask** — commits the selection. The selection band disappears, the
   masked rows become NaN in every channel, the traces in that interval
   render in grey, and a soft grey background band is drawn across all panels.
6. **Unmask** — opposite of Mask: clears the mask inside the current selection.
7. **Clear** — removes every mask interval.
8. **Write** — writes two files next to the loaded file:
   - `<stem>_clean.<ext>` — the cleaned series in the **same format** it
     was loaded from (LEMI-424 24-col text, or 7-col `.xyz`). Masked rows
     carry `NaN`.
   - `<stem>_mask.csv` — the mask intervals as `start,stop` datetimes.

You can mask many intervals, unmask parts of them, and only call **Write**
once at the end. Reloading `<stem>_clean.<ext>` shows the same record with
the masked windows as gaps in the lines — useful for verifying the export
before running downstream processing.

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

## Core API

- `run_tkapp()` / `run_tkapp(path)` / `run_tkapp(ta)` — open the GLMakie window
- `TKApp(...)` — build an app without blocking (for embedding/testing)
- `load_lemi424(path)` — returns a five-component `TimeArray` (`bx/by/bz/e1/e2`)
- `read_lemi424(path)` — returns a full `TimekeeperRun`
- `write_lemi424(path, run_or_timearray)` — writes LEMI-424 24-col text
- `read_timekeeper(path)` / `write_timekeeper(path, run)` — format dispatch
- `TimekeeperMask(timearray)` — all-component processing mask
- `mask_interval!` / `unmask_interval!` — edit by `[start, end]` datetimes
- `clear_mask!` / `masked_samples` — bulk ops and status
- `cleaned_timearray(ta, mask)` — preserves timestamps, masked rows become `NaN`
- `good_segments(ta, mask; min_samples)` — contiguous unmasked chunks
- `sample_weights(mask; good, bad)` — vector of robust-processing weights
- `write_cleaned(path, ta, mask)` / `write_mask(path, mask)` — CSV exports

## TKApp design notes

- Five linked panels (`Bx`, `By`, `Bz` in black; `Ex`, `Ey` in blue) share a
  single x-axis. Date ticks appear only on the bottom panel.
- All five panels share a single sample-aligned mask. Selecting and masking
  any panel marks the same time interval on all five — so processing can
  treat the cleaned `TimeArray` as a coherent multi-channel record.
- The horizontal slider, not the mouse, controls the visible time window;
  this leaves left-drag free for selection without accidental rectangle zoom.
- Multiple mask intervals can be accumulated before calling Write.
- A `View:` menu toggles between **Time** (the five linked TS panels alone),
  **Time | Spectra** (each TS panel paired with a log-log PSD panel for the
  same channel) and **Time | Spectrogram** (each TS panel paired with a
  log-frequency / linear-time STFT heatmap). PSDs use Welch's method on the
  unmasked segments inside the currently visible window: a Hann window of
  length `nfft` with 50% overlap, periodograms averaged across all good
  segments. `nfft` is auto-derived from the window length (≈ 7 sub-windows of
  overlap) and a small header shows `nfft / f_min / T_max` next to the right
  column. Spectra and spectrograms recompute on scroll, on window changes,
  and on every Mask / Unmask / Clear; they do **not** recompute while a
  drag-selection is in progress.

## Project layout

```
Timekeepers.jl/
├── data/                  # gitignored, drop test files here
├── examples/
│   ├── tkapp.jl           # `run_tkapp()` launcher
│   └── read_lemi424.jl
├── images/                # README screenshots
├── src/
│   ├── Timekeepers.jl     # module + exports
│   ├── Explorer.jl        # TKApp / GLMakie UI
│   ├── LEMI424.jl         # LEMI-424 reader/writer
│   ├── Masking.jl         # TimekeeperMask and derived series
│   ├── TimeArrayIO.jl     # TimeArray <-> TimekeeperRun glue
│   ├── TimekeeperIO.jl    # format dispatch
│   ├── Types.jl           # core types
│   └── Utilities.jl       # helpers
└── README.md
```
