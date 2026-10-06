# Getting Started

## Install

These steps assume that you have not used Julia much. You do them once.

### 1. Install Julia

Timekeepers.jl needs Julia 1.12 or newer. Install it with juliaup, the official installer:

- **macOS and Linux**: in a terminal, run `curl -fsSL https://install.julialang.org | sh`
- **Windows**: in a terminal, run `winget install --name Julia --id 9NJNWW8PVKMN -e -s msstore`

Close the terminal, open a new one, and check the version with `julia --version`.

### 2. Get Timekeepers.jl

Timekeepers.jl is not in the Julia package registry yet, so get it from GitHub. With git:

```bash
git clone https://github.com/JuliaGeophysics/Timekeepers.jl
```

Without git, open <https://github.com/JuliaGeophysics/Timekeepers.jl>, click **Code**, then
**Download ZIP**, and unzip it.

### 3. Install its packages

In a terminal, go into the folder and start Julia in it:

```bash
cd Timekeepers.jl
julia --project=. -t auto
```

`--project=.` tells Julia to use the packages of this folder. `-t auto` lets the windows stay live
while they work.

Julia shows the `julia>` prompt. Type `]` (a closing square bracket). The prompt changes to
`(Timekeepers) pkg>`: this is the package mode. Type `instantiate` and press Enter:

```julia
(Timekeepers) pkg> instantiate
```

Julia downloads and compiles about 290 packages, most of them for the graphics. This takes several
minutes, only the first time. When it is done, press Backspace to go back to `julia>`.

### 4. Check that windows open

TKDash, TKApp and TKProc open windows, so they need a desktop with OpenGL 3.3 or newer. At the
`julia>` prompt, type:

```julia
using GLMakie
display(scatter(1:10))
```

If a window with ten dots opens, you are ready.

!!! tip "Use Timekeepers.jl in your own project"
    If you have your own Julia project, you can add Timekeepers.jl to it instead of steps 2 and 3.
    Start Julia in your project, type `]`, then:

    ```julia
    pkg> add https://github.com/JuliaGeophysics/Timekeepers.jl
    ```

## Use

Each time, start Julia in the Timekeepers.jl folder as in step 3 (`julia --project=. -t auto`),
then load the package and open a window:

```julia
using Timekeepers

run_tkdash("data/survey")          # 1. see which sites recorded together
run_tkapp("data/LEMI090.txt")      # 2. clean one recording or site
run_tkproc("data/survey/site004")  # 3. estimate its transfer function
```

Replace the paths with the paths of your data. Each window stays open until you close it. Then the
`julia>` prompt comes back.

- [TKDash](tkdash.md): find the base and remote sites of each site in a survey
- [TKApp](tkapp.md): mask bad intervals and write clean data
- [TKProc](tkproc.md): estimate the impedance and the tipper, and write EDI and ModEM files

## Test data

- LEMI-424: [British Geological Survey accession](https://webapps.bgs.ac.uk/services/ngdc/accessions/index.html#item182849)
- Metronix ADU survey: [BRGM razorback tutorial data](https://github.com/BRGM/razorback-tutorial-data)
