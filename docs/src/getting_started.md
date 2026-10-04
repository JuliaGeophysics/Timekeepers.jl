# Getting Started

## Installation

Timekeepers needs Julia 1.12 or newer.

### As a package

Use this procedure to use Timekeepers.jl in your own project or scripts. The
package is not registered yet. Thus, add it from GitHub into a dedicated
project environment:

```julia
julia> ]  # press ] to enter the Pkg REPL
pkg> activate @timekeepers   # a named shared environment; or `activate .` for the current folder
pkg> add https://github.com/JuliaGeophysics/Timekeepers.jl
```

You can also do the same steps with one command:

```bash
julia --project=@timekeepers -e 'using Pkg; Pkg.add(url = "https://github.com/JuliaGeophysics/Timekeepers.jl")'
```

!!! tip
    Do not install packages into your default (global) environment. Use a
    dedicated environment for each project. This keeps the dependencies
    isolated and reproducible. It also prevents slow, unexpected changes to
    the versions of other packages that you installed.

### From a clone

Use this procedure to run the examples in the `examples/` directory, or to
develop the package:

```bash
git clone https://github.com/JuliaGeophysics/Timekeepers.jl.git
cd Timekeepers.jl/
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

These docs use the prefix `julia --project=.`. This prefix activates the
environment of the clone.

## Checking your OpenGL setup

TKApp is a native GLMakie window. Thus, it needs a desktop session with
OpenGL 3.3 or newer. If this code opens a window, [`run_tkapp`](@ref) also
opens a window:

```julia
using GLMakie
display(scatter(1:10))
```

Headless SSH sessions, WSL, containers and very old GPUs usually need more
configuration. For example, use a virtual framebuffer such as `xvfb-run`, or
a forwarded display. All the functions outside the app (readers, writers and
masks) operate correctly without a display.

## The two data containers

Timekeepers has one native container and one interop container. Most
procedures move the data between the two containers.

Each reader returns a [`TimekeeperRun`](@ref). The run holds a `Dict` of
[`TimekeeperChannel`](@ref)s. It also holds the metadata that a writer needs
to make the original file again: the header fields, the position and, for
Metronix, the paths of the XML templates. Keep the run if you will write the
data again.

The app, the mask functions and most processing steps use a
`TimeSeries.TimeArray`. A `TimeArray` is a matrix of samples with one time
axis. Use [`to_timearray`](@ref) and [`from_timearray`](@ref) to move the data
between the two containers.

```julia
using Timekeepers

run = read_timekeeper("data/LEMI090.txt")   # TimekeeperRun
ta  = to_timearray(run)                     # TimeArray, columns bx by bz e1 e2

components(run)        # [:bx, :by, :bz, :e1, :e2, ...]
sampling_rate(run)     # 1.0
start_time(run)        # 2020-10-04T00:00:00
duration_seconds(run)  # 86400.0
```

Each reader also has a `load_*` function that returns a `TimeArray` directly.
Use it when you do not need the run:

```julia
ta = load_lemi424("data/LEMI090.txt")
```

## Reading without knowing the format

[`read_timekeeper`](@ref) finds the format from the path. For a text file, it
also uses the first line that is not blank:

- A `.ats` file, or a directory that contains one, is Metronix.
- A `.txt` file with a `GEOMAG` header is GEOMAG-02.
- All other files are LEMI-424.

```julia
run = read_timekeeper("data/MS_26_250523000000.TXT")
run.source_format   # :geomag
```

To skip the detection, give the `format` argument:

```julia
run = read_timekeeper("data/oddly_named_file.dat"; format = :lemi424)
```

## A first round trip

Each writer does the opposite of its reader. [`write_timekeeper`](@ref) uses
the `source_format` of the run to select the correct writer:

```julia
run = read_timekeeper("data/LEMI090.txt")
write_timekeeper("data/LEMI090_copy.txt", run)
```

The reader keeps the auxiliary columns in the metadata. These columns are the
temperatures, the battery voltage, the GPS fix and the satellite count. The
writer puts each column back in its original position. Thus, the copy has the
same structure as the input, byte for byte.

## Reading a directory of runs

A field deployment usually makes many files. Use the same reader on each
file:

```julia
using Timekeepers

files = filter(f -> endswith(lowercase(f), ".txt"), readdir("data/SITE01"; join = true))
runs  = [read_timekeeper(f) for f in files]
```

For Metronix sites, use the dedicated index [`metronix_site_runs`](@ref). It
puts the runs in groups by sampling rate. Refer to [Metronix Sites](metronix.md).

TKApp can do this work for you. **Load Site…** reads all the runs in a
directory and puts them in order of start time. It then joins the runs and
puts `NaN` in each gap between runs. Thus, a month of hourly files becomes one
continuous record on the screen.

## Next steps

- [Instrument Formats](formats.md): the data that each reader and writer
  uses.
- [Masking & Cleaning](masking.md): how to mark bad intervals and make clean
  data.
- [TKApp Explorer](tkapp.md): the interactive window.
- [TKDash Survey](tkdash.md): the base and remote sites of a survey.
- [API Reference](api.md): all the exported functions.
