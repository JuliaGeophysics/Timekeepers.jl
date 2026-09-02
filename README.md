<div align="center">
  <img src="images/timekeepers-logo.svg" alt="Timekeepers.jl" width="120">

  # Timekeepers.jl

  *Time-series I/O and interactive inspection for magnetotelluric field data.*

  [![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://JuliaGeophysics.github.io/Timekeepers.jl/stable)
  [![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://JuliaGeophysics.github.io/Timekeepers.jl/dev)
  [![CI](https://github.com/JuliaGeophysics/Timekeepers.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/JuliaGeophysics/Timekeepers.jl/actions/workflows/CI.yml)
  [![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
</div>

Timekeepers reads logger-native files into [TimeSeries.jl](https://github.com/JuliaStats/TimeSeries.jl)
`TimeArray`s, keeps native writer paths for supported instruments, and ships
**TKApp**, a native GLMakie window for fast visual inspection and masking of
long records.



## Installation

Requires Julia 1.10 or newer.

**As a package** — to use Timekeepers.jl from your own project or scripts. It is
registered in the Julia General registry.

Install it into a dedicated project environment:

```julia
julia> ]  # press ] to enter the Pkg REPL
pkg> activate @timekeepers   # a named shared environment; or `activate .` for the current folder
pkg> add Timekeepers
```

or equivalently, non-interactively:

```bash
julia --project=@timekeepers -e 'using Pkg; Pkg.add("Timekeepers")'
```

> [!TIP]
> As a general Julia best practice, avoid installing packages into your default
> (global) environment. A dedicated per-project environment keeps dependencies
> isolated and reproducible, and avoids slow, unexpected version changes across
> unrelated packages you already have installed.

**From a clone** — to run the bundled examples and helper scripts, or to develop
the package:

```bash
git clone https://github.com/JuliaGeophysics/Timekeepers.jl.git
cd Timekeepers.jl/
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

The `julia --project=.` prefix used in the examples below activates that cloned
environment.

TKApp uses GLMakie, so you also need a desktop session with OpenGL 3.3+
drivers. Quick smoke test — if this opens a window, `run_tkapp()` will too:

```julia
using GLMakie
display(scatter(1:10))
```

Headless SSH, WSL, containers, and very old GPUs may need extra OpenGL setup.
Everything outside the app works headless.

## Quick Start

Open the explorer window:

```julia
using Timekeepers
run_tkapp()
```

or run the bundled launcher script from a clone:

```powershell
julia --project=. examples/tkapp.jl
```


Hovering a PSD panel snaps a cursor to the nearest peak and reads it out over
the trace as `[0.1 Hz / 10s; PSD=1.234e-03]` — frequency, period, power.
**Left-click** pins that frequency across every channel panel and draws faint
guides at `2f0`, `3f0`, … so harmonics of the pinned tone can be told apart from
unrelated peaks; **right-click** clears the pin. The pin is labelled the same
way, `[0.1 Hz / 10s]`, and drawn in its own colour to tell it from the cursor.

The parameter line under the plots reports the spectral configuration in use.
The **i** badge at its left swaps it for a plain-language gloss of each term.

The **View** menu offers `Time`, `Spectra` — the same PSD panels at full width,
with the traces off — and `Time | Spectra` side by side.

Metronix records load either way round: **Load Run…** opens one `meas_*` run
(navigate in and pick any `.ats`; the whole run loads with it), while **Load
Site…** takes the site above them, asking which sampling rate to load when the
site holds more than one.

The same workflow from Julia code:

```julia
using Timekeepers, Dates

# Load
ta = load_lemi424("data/LEMI090.txt")                 # → TimeArray

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

## Documentation

Full documentation, including the format reference, masking workflow and
Metronix site surgery, lives at
**<https://JuliaGeophysics.github.io/Timekeepers.jl/stable>**.

## Contributing

Issues and pull requests are welcome — see [CONTRIBUTING.md](CONTRIBUTING.md).
