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

```julia
pkg> add Timekeepers
```

Requires Julia 1.10 or newer. TKApp uses GLMakie, so you need a desktop session
with OpenGL 3.3+ drivers. Quick smoke test — if this opens a window,
`run_tkapp()` will too:

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
