# Getting Started

Four steps, done once. Then you're ready to [use it](#Use).

## 1. Install Julia (1.12 or newer)

- **macOS / Linux**: `curl -fsSL https://install.julialang.org | sh`
- **Windows**: `winget install --name Julia --id 9NJNWW8PVKMN -e -s msstore`

Open a new terminal and check: `julia --version`.

## 2. Download Timekeepers.jl

```bash
git clone https://github.com/JuliaGeophysics/Timekeepers.jl
```

No git? On <https://github.com/JuliaGeophysics/Timekeepers.jl> click **Code → Download ZIP** and
unzip it.

## 3. Install its packages

```bash
cd Timekeepers.jl
julia --project=. -t auto
```

At the `julia>` prompt, type `]`, then `instantiate`, then press Backspace when it finishes:

```julia
(Timekeepers) pkg> instantiate
```

The first time takes several minutes (about 290 packages, mostly graphics).

## 4. Check that windows open

The windows need a desktop with OpenGL 3.3 or newer. Check with:

```julia
using GLMakie
display(scatter(1:10))
```

A window with ten dots means you're ready.

!!! tip "Already have a Julia project?"
    Skip steps 2 and 3. In your project, type `]` and run
    `add https://github.com/JuliaGeophysics/Timekeepers.jl`.

## Use

Start Julia in the Timekeepers.jl folder with `julia --project=. -t auto`, then:

```julia
using Timekeepers

run_tkdash("data/survey")          # 1. see which sites recorded together
run_tkapp("data/LEMI090.txt")      # 2. mask bad intervals, write clean data
run_tkproc("data/survey")          # 3. estimate transfer functions, write EDI files
```

Use your own paths. Each window stays open until you close it.

| Next | |
|:---|:---|
| [TKDash](tkdash.md) | find the base and remote sites of each site |
| [TKApp](tkapp.md) | clean one recording or site |
| [TKProc](tkproc.md) | impedance, tipper and EDI files |

## Test data

- LEMI-424: [British Geological Survey accession](https://webapps.bgs.ac.uk/services/ngdc/accessions/index.html#item182849)
- Metronix ADU survey: the public Metronix tutorial data of BRGM (Alsace, IMAGE FP7 project)
