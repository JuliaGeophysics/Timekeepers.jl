# Getting Started

## Install

Timekeepers.jl needs Julia 1.12 or newer. It is not in the General registry yet, so add it from
GitHub, preferably in its own environment:

```julia
pkg> activate @timekeepers
pkg> add https://github.com/JuliaGeophysics/Timekeepers.jl
```

The first install downloads and precompiles about 290 packages, most of them for the graphics. This
takes a few minutes once.

TKDash and TKApp open windows, so they need a desktop with OpenGL 3.3 or newer. If this opens a
window, you are ready:

```julia
using GLMakie
display(scatter(1:10))
```

## Use

```julia
using Timekeepers

run_tkdash("data/survey")          # 1. see which sites recorded together
run_tkapp("data/LEMI090.txt")      # 2. clean one recording or site
```

- [TKDash](tkdash.md): find the base and remote sites of each site in a survey
- [TKApp](tkapp.md): mask bad intervals and write clean data

## Test data

- LEMI-424: [British Geological Survey accession](https://webapps.bgs.ac.uk/services/ngdc/accessions/index.html#item182849)
- Metronix ADU survey: [BRGM razorback tutorial data](https://github.com/BRGM/razorback-tutorial-data)
