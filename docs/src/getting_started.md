# Getting Started

## Installation

```julia
pkg> add Timekeepers
```

Timekeepers requires Julia 1.10 or newer. To work from a clone instead — which
also gives you the `examples/` and `scripts/` directories referenced throughout
these docs:

```bash
git clone https://github.com/JuliaGeophysics/Timekeepers.jl
cd Timekeepers.jl
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

## Checking your OpenGL setup

TKApp is a native GLMakie window, so it needs a desktop session with OpenGL 3.3
or newer. If this opens a window, [`run_tkapp`](@ref) will too:

```julia
using GLMakie
display(scatter(1:10))
```

Headless SSH sessions, WSL, containers and very old GPUs usually need extra
setup — a virtual framebuffer such as `xvfb-run`, or a forwarded display.
Everything outside the app (readers, writers, masking) works fine headless.

## The two data containers

Timekeepers has one native container and one interop container, and most
workflows move between them.

[`TimekeeperRun`](@ref) is what every reader returns. It holds a
`Dict` of [`TimekeeperChannel`](@ref)s plus the run metadata a writer needs to
reproduce the original file — header fields, position, and for Metronix the XML
template paths. Keep a run around when you intend to write the data back out.

`TimeSeries.TimeArray` is what the app, the masking layer and most downstream
processing work with: a matrix of samples with a shared time axis. Move between
the two with [`to_timearray`](@ref) and [`from_timearray`](@ref).

```julia
using Timekeepers

run = read_timekeeper("data/LEMI090.txt")   # TimekeeperRun
ta  = to_timearray(run)                     # TimeArray, columns bx by bz e1 e2

components(run)        # [:bx, :by, :bz, :e1, :e2, ...]
sampling_rate(run)     # 1.0
start_time(run)        # 2020-10-04T00:00:00
duration_seconds(run)  # 86400.0
```

Every reader also has a `load_*` twin that goes straight to a `TimeArray`, for
when you do not need the run:

```julia
ta = load_lemi424("data/LEMI090.txt")
```

## Reading without knowing the format

[`read_timekeeper`](@ref) infers the format from the path and, for text files,
from the first non-blank line: a `.ats` file or a directory containing one is
Metronix, a `.txt` file with a `GEOMAG` header is GEOMAG-02, and anything else
is treated as LEMI-424.

```julia
run = read_timekeeper("data/MS_26_250523000000.TXT")
run.source_format   # :geomag
```

Pass `format` explicitly to skip detection:

```julia
run = read_timekeeper("data/oddly_named_file.dat"; format = :lemi424)
```

## A first round trip

The writers are the inverse of the readers, and [`write_timekeeper`](@ref)
picks the right one from the run's own `source_format`:

```julia
run = read_timekeeper("data/LEMI090.txt")
write_timekeeper("data/LEMI090_copy.txt", run)
```

Auxiliary columns — temperatures, battery voltage, GPS fix, satellite count —
are carried through the reader in metadata and written back in their original
slots, so the copy is byte-comparable in structure to the input.

## Reading a directory of runs

A field deployment usually produces many files. Loop over them with the same
reader:

```julia
using Timekeepers

files = filter(f -> endswith(lowercase(f), ".txt"), readdir("data/SITE01"; join = true))
runs  = [read_timekeeper(f) for f in files]
```

For Metronix sites there is a dedicated index — [`metronix_site_runs`](@ref)
groups the `meas_*` directories by sampling rate — see
[Metronix Sites](metronix.md).

TKApp does this for you: **Load Site…** reads every run in a directory, orders
them by start time and concatenates them with `NaN` filling any gap between
runs, so a month of hourly files becomes one continuous record on screen.

## Next steps

- [Instrument Formats](formats.md) — what each reader and writer handles.
- [Masking & Cleaning](masking.md) — marking bad intervals and deriving clean data.
- [TKApp Explorer](tkapp.md) — the interactive window.
- [API Reference](api.md) — every exported function.
