# Metronix Sites

A Metronix ADU deployment produces a *site* directory holding many `meas_*`
measurement directories, often recorded at several sampling rates in one
campaign — a long 128 Hz run interleaved with short 4096 Hz bursts. This page
covers the two operations that turn such a site into something processable:
splitting it by rate, and amputating masked intervals.

## Anatomy of a site

```text
RK137/
├── meas_2025-04-01_07-00-05/
│   ├── 076_V01_C00_R000_TEx_BL_128H.ats
│   ├── 076_V01_C01_R000_TEy_BL_128H.ats
│   ├── ...
│   └── 076_2025-04-01_07-00-06_..._R000_128H.xml
├── meas_2025-04-01_09-00-05/
└── ...
```

Each `meas_*` directory is one continuous run: an `.ats` binary per channel
plus a shared `.xml` sidecar. [`read_metronix`](@ref) reads one of these; the
functions below work on the site directory above them.

## Indexing a site

[`metronix_site_rates`](@ref) reports the distinct sampling rates present, and
[`metronix_site_runs`](@ref) groups the measurement directories by rate:

```julia
using Timekeepers

is_metronix_site("data/RK137")     # true
metronix_site_rates("data/RK137")  # [128.0, 4096.0]

runs = metronix_site_runs("data/RK137")
length(runs[128.0])                # 14 measurement directories at 128 Hz
```

Both read only the `.ats` headers, so indexing a large site is fast.

## Splitting a mixed-rate site

Most of Timekeepers assumes one sampling rate per site — [`sampling_rate`](@ref)
errors on a run whose channels disagree, and the app plots a single time axis.
Split first.

The bundled script does this from the command line, copying each `meas_*`
directory verbatim into a sibling `<site>.TK<rate>` directory:

```bash
julia --project=. scripts/split_metronix_by_rate.jl data/RK137
```

```text
RK137/  ->  RK137.TK128/     (all 128 Hz meas_ dirs)
            RK137.TK4096/    (all 4096 Hz meas_ dirs)
```

Pass `--move` to move rather than copy when disk space is tight. The script
also defines `split_metronix_site_by_rate`, usable from Julia after `include`:

```julia
include("scripts/split_metronix_by_rate.jl")
dests = split_metronix_site_by_rate("data/RK137")            # copy
dests = split_metronix_site_by_rate("data/RK137"; move = true)
```

The copy is byte-for-byte: the split is a regrouping, not a rewrite.

## Amputating masked intervals

Metronix has no `NaN`, and processing tools expect continuous runs. So a mask
on a Metronix record is not applied by blanking samples — it is applied by
*cutting*, splitting each run at the masked intervals and writing the surviving
stretches as separate `meas_*` directories.

### A whole site on disk

[`write_metronix_site_masked`](@ref) applies a list of `DateTime` intervals
across every run of a single-rate site:

```julia
using Timekeepers, Dates

intervals = [
    (DateTime(2025, 4, 1, 7, 0, 8), DateTime(2025, 4, 1, 7, 0, 10)),
    (DateTime(2025, 4, 1, 9, 32, 0), DateTime(2025, 4, 1, 9, 41, 0)),
]

dest = write_metronix_site_masked("data/RK137.TK128"; intervals = intervals)
# "data/RK137.TK128.W"
```

Runs no interval touches are written whole. Runs that are touched are split;
each surviving stretch becomes `meas_<its own start time>`, with the XML
sidecar regenerated to match — corrected sample counts, start times and
filenames. Segments shorter than `min_samples` are dropped, and segment starts
are trimmed to whole seconds where the sampling rate requires it.

The destination defaults to `<site_dir>.W` (`RK137.TK128` → `RK137.TK128.W`),
so the source is never modified. This is what the app's **Write** button calls
when a Metronix site is loaded.

### A single loaded run

When you already have a run and a mask in memory,
[`write_metronix_site`](@ref) does the same split for that one run:

```julia
run  = read_metronix("data/RK137.TK128/meas_2025-04-01_07-00-05")
ta   = to_timearray(run)
mask = TimekeeperMask(ta)
mask_interval!(mask, DateTime(2025, 4, 1, 7, 30), DateTime(2025, 4, 1, 7, 35))

dest, dirs = write_metronix_site(run; mask = mask)
length(dirs)   # 2 — the record either side of the cut
```

With `mask = nothing` the run is written whole. To write a single measurement
directory with no splitting at all, use [`write_metronix`](@ref).

## The write log

Each call to [`write_metronix_site_masked`](@ref) appends a datetime-stamped
section to a `README.md` in the destination, naming the source site and every
interval that was cut. Repeated writes accumulate sections rather than
overwriting, so the destination carries its own provenance:

```markdown
## Write session 2025-05-20T11:07:33

Source site: `data/RK137.TK128`

Masked intervals:
- 2025-04-01T07:00:08 → 2025-04-01T07:00:10
```

That log is the reason to keep the `.W` directory as the thing you hand to a
processing chain: it is reproducible from the source plus the record of what
was cut.

## Round-trip guarantee

Writing a site with no intervals reproduces the input exactly — the `.ats`
files come back byte-identical, and re-reading gives the same samples and start
times. An unmasked write is a lossless copy, so it is safe to route every site
through the writer whether or not it needed editing.
