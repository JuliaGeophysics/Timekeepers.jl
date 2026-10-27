<div align="center">

  # Timekeepers.jl

  *Timeseries analysis in Electromagnetic Geophysics*

  [![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://JuliaGeophysics.github.io/Timekeepers.jl/dev)
  [![CI](https://github.com/JuliaGeophysics/Timekeepers.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/JuliaGeophysics/Timekeepers.jl/actions/workflows/CI.yml)
  [![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
</div>

Timekeepers reads LEMI-424, GEOMAG-02 and Metronix ADU recordings into
[TimeSeries.jl](https://github.com/JuliaStats/TimeSeries.jl) `TimeArray`s. It
also writes each of these formats. It has three windows:

- **TKDash**, a survey dashboard to select base and remote sites;
- **TKApp**, a GLMakie window to examine, mask and clean long records;
- **TKProc**, a window to estimate transfer functions (experimental).

<p align="center">
  <img src="docs/src/assets/TKDash.png" alt="TKDash showing the sites of a survey on a map" width="32%">
  <img src="docs/src/assets/TK.png" alt="TKApp showing a five-channel Metronix record" width="32%">
  <img src="docs/src/assets/TKproc.png" alt="TKProc showing the apparent resistivity, phase and tipper of a site" width="32%">
</p>

## Installation

Timekeepers needs Julia 1.12 or newer, and is not registered yet. If you
already have Julia and your own project, add it from GitHub:

```julia
pkg> add https://github.com/JuliaGeophysics/Timekeepers.jl
```

### No Julia yet?

1. **Install Julia** with [juliaup](https://julialang.org/install/), the
   official installer:
   - macOS/Linux: `curl -fsSL https://install.julialang.org | sh`
   - Windows: `winget install --name Julia --id 9NJNWW8PVKMN -e -s msstore`

   Open a new terminal and check it with `julia --version`.
2. **Clone this repository** (or download it as a ZIP from the **Code**
   button on GitHub and unzip it):
   ```bash
   git clone https://github.com/JuliaGeophysics/Timekeepers.jl
   cd Timekeepers.jl
   ```
3. **Install its packages**:
   ```bash
   julia --project=. -t auto -e "using Pkg; Pkg.instantiate()"
   ```
   This downloads and compiles everything Timekeepers needs, mostly
   graphics. It takes a few minutes, only the first time.

TKApp, TKDash and TKProc open windows, so they need a desktop session with
OpenGL 3.3 or newer. The readers and writers operate without a display.

The full walkthrough, with more detail for a first-time Julia user, is in
[Getting Started](https://JuliaGeophysics.github.io/Timekeepers.jl/dev/getting_started).

## Quick start

From the `Timekeepers.jl` folder, start Julia with `julia --project=. -t auto`,
then:

```julia
using Timekeepers
run_tkapp()                                   # open the explorer window
```

Or, without starting Julia yourself, run the launcher scripts in `examples/`:

```bash
julia --project=. -t auto examples/tkapp.jl
```

For a full survey, **TKDash** shows when each site recorded. For the site that
you select, it shows the base sites (they recorded with it, near it) and the
remote sites (they recorded with it, far from it). It exports these sites for
each site:

```julia
run_tkdash()                                  # or: julia --project=. examples/tkdash.jl
```

Without a path, `run_tkdash()` opens a folder dialog that asks for the survey
directory; give it a path (`run_tkdash("path/to/survey")`) to skip the dialog.

You can also use code:

```julia
using Timekeepers, Dates

ta   = load_lemi424("data/LEMI090.txt")
mask = TimekeeperMask(ta)
mask_interval!(mask, DateTime(2020, 10, 4, 0, 10), DateTime(2020, 10, 4, 0, 20))

write_lemi424("data/LEMI090_clean.txt", cleaned_timearray(ta, mask))
write_mask("data/LEMI090_mask.csv", mask)
```

## Documentation

<https://JuliaGeophysics.github.io/Timekeepers.jl/dev>

## Contributing

We accept issues and pull requests. Refer to [CONTRIBUTING.md](CONTRIBUTING.md).
