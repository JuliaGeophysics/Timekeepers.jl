<div align="center">

  # Timekeepers.jl

  *Timeseries analysis in Electromagnetic Geophysics*

  [![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://JuliaGeophysics.github.io/Timekeepers.jl/dev)
  [![CI](https://github.com/JuliaGeophysics/Timekeepers.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/JuliaGeophysics/Timekeepers.jl/actions/workflows/CI.yml)
  [![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
</div>

Timekeepers reads LEMI-424, GEOMAG-02 and Metronix ADU recordings into
[TimeSeries.jl](https://github.com/JuliaStats/TimeSeries.jl) `TimeArray`s, writes
each format back, and ships **TKApp**, a GLMakie window for inspecting, masking
and cleaning long records, and **TKDash**, a survey dashboard for choosing
base and remote sites.

![TKApp showing a five-channel Metronix record](docs/src/assets/TK.png)

## Installation

Requires Julia 1.12 or newer. Timekeepers is not yet registered, so add it from
GitHub:

```julia
pkg> add https://github.com/JuliaGeophysics/Timekeepers.jl
```

TKApp needs a desktop session with OpenGL 3.3 or newer; everything else works
headless.

## Quick start

```julia
using Timekeepers
run_tkapp()                                   # open the explorer window
```

For a whole survey, **TKDash** shows when every site recorded and, for the
site you pick, its base sites (recorded with it, close by) and remote
sites (recorded with it, far away), and exports them for each site:

```julia
run_tkdash("path/to/survey")                  # or: julia --project=. examples/tkdash.jl <dir>
```

Or from code:

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

Issues and pull requests are welcome; see [CONTRIBUTING.md](CONTRIBUTING.md).
