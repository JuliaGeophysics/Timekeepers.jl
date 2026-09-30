# Metronix Sites

A Metronix ADU deployment produces a *site* directory holding `meas_*`
measurement directories, often recorded at several sampling rates in one
campaign — a long 128 Hz run interleaved with short 4096 Hz bursts. This page
covers how Timekeepers finds the runs in such a site, loads them one rate at a
time, separates them by sampling rate on disk, and amputates masked intervals.

## Anatomy of a site

```text
DF002/
├── meas_2021-09-25_13-55-25/
│   ├── 081_V01_C00_R000_TEx_BH_131072H.ats
│   ├── ...
│   └── 081_2021-09-25_13-55-25_2021-09-25_14-00-25_R000_131072H.xml
└── meas_2021-09-25_14-02-01/
    ├── 081_V01_C00_R000_TEx_BL_128H.ats          ┐ R000 at 128 Hz
    ├── ...                                        │
    ├── 081_2021-09-25_14-02-02_..._R000_128H.xml  ┘
    ├── 081_V01_C00_R000_TEx_BL_4096H.ats         ┐ R000 at 4096 Hz
    ├── ...                                        │
    ├── 081_2021-09-26_00-00-00_..._R000_4096H.xml ┘
    ├── 081_V01_C00_R001_TEx_BL_4096H.ats         ┐ R001 at 4096 Hz
    ├── ...                                        │
    ├── 081_2021-09-27_00-00-00_..._R001_4096H.xml ┘
    ├── 081_2021-09-28_00-00-00_..._R002_4096H.xml   scheduled, never recorded
    └── Site_meas_2021-09-25_14-02-01.kml
```

A *run* is one run number (`R000`, `R001`, …) at one sampling rate: an `.ats`
binary per channel plus the `.xml` that describes them. A `meas_*` directory
can hold one run or several — the ADU writes a long low-rate run and the
high-rate bursts scheduled during it into the same directory. The run number
and rate in the filenames tie each `.ats` to its XML.

An XML with no `.ats` files beside it is a job the ADU scheduled but never
recorded. Timekeepers skips it when reading, and copies it unchanged when
writing.

## Indexing a site

[`metronix_site_rates`](@ref) reports the distinct sampling rates present, and
[`metronix_site_runs`](@ref) lists the runs of each rate in start-time order,
naming each by the path of its XML:

```julia
using Timekeepers

is_metronix_site("data/DF002")     # true
metronix_site_rates("data/DF002")  # [128.0, 4096.0, 131072.0]

runs = metronix_site_runs("data/DF002")
basename.(runs[4096.0])
# "081_2021-09-26_00-00-00_2021-09-26_02-00-00_R000_4096H.xml"
# "081_2021-09-27_00-00-00_2021-09-27_02-00-00_R001_4096H.xml"
```

Both read only the `.ats` headers, so indexing a large site is fast. Both
accept a single `meas_*` directory as well as a site.

[`read_metronix`](@ref) reads one run, named by its XML, by any of its `.ats`
files, or by a `meas_*` directory that holds only that run:

```julia
run = read_metronix(runs[4096.0][2])
sampling_rate(run)  # 4096.0
```

## Loading in the app

**Load Run…** opens one run: navigate into the `meas_*` directory and pick the
run's `.xml`, or any of its `.ats` files. That run — its channels plus the XML
— loads, and nothing else in the directory does.

**Load Site…** takes the site directory, or a single `meas_*` directory. It
finds every run and groups them by sampling rate. The app holds one rate in
memory at a time — high-rate runs are large — so for a mixed-rate site it asks
which rate to import, listing each with its number of runs and noting any XML
skipped for having no data. Closing that window cancels the import.

Before reading, the loading window checks whether the site is already split
by sampling rate into `<site>.TK` beside it (see
[Separating a site by rate](#Separating-a-site-by-rate)). If every file is
there it says so and copies nothing; otherwise it splits the site, copying
only what is missing. The runs are then read from `<site>.TK`. Picking a
`.TK` directory, or one of its rate directories, loads it directly.

The runs of the chosen rate appear as one record in time order. The rate menu
at the left of the toolbar then switches rates: picking one reads its runs from
disk in place of those on screen. Masks survive the switch — the rate you leave
keeps its masked intervals, and gets them back when you return — and **Write**
writes the whole site into a new `<site>.TK<date>_<time>` directory: the rate
on screen and every other rate you masked cut by their own intervals, the
rest copied as they are.

From Julia, [`load_metronix_site`](@ref) reads one rate of a site the same way,
separating it by rate first; `rate` is required when the site holds more than
one:

```julia
overview = load_metronix_site("data/DF002"; rate = 4096)
```

### Rates above 1 kHz

Up to 1 kHz the runs of a rate share one time grid, the gaps between them
filled with `NaN`. Above 1 kHz several samples fall within each millisecond,
finer than the app's time stamps resolve, so the runs are instead joined end to
end: every sample is kept, the gap between two runs is a jump in the time
axis, and the plot breaks its line there. The spectra treat that gap like a
mask: no FFT segment spans it, and the spectrum of a window covering several
runs is the average of the segments each run contributes. The line under the
plots says how many runs went into it.

A high-rate run is large in memory — two hours at 4096 Hz is about 30 million
samples per channel — so expect a few seconds and a few GB per run loaded.

## Separating a site by rate

[`split_metronix_site`](@ref) copies a site into `<site>.TK`, one directory per
sampling rate, each holding the `meas_*` directories with runs at that rate:

```text
DF002.TK/
├── 128/
│   └── meas_2021-09-25_14-02-01/      R000 at 128 Hz, its XML, the .kml
├── 4096/
│   └── meas_2021-09-25_14-02-01/      R000, R001 at 4096 Hz, their XMLs,
│                                      R002's XML (never recorded), the .kml
└── 131072/
    └── meas_2021-09-25_13-55-25/      R000 at 131072 Hz, its XML, the .kml
```

Each run's `.ats` files and XML go to its rate. The XML of a job that never
recorded goes to the rate in its filename, and the `.kml` to every rate the
`meas_*` directory holds. Every file is copied byte for byte and the source is
left alone; putting the rate directories back together gives the original
site file for file. Loading a site does this for you.

## Amputating masked intervals

Metronix has no `NaN`, and processing tools expect continuous runs. So a mask
on a Metronix record is not applied by blanking samples — it is applied by
*cutting*, splitting each run at the masked intervals and writing the surviving
stretches as separate `meas_*` directories.

### A whole site on disk

[`write_metronix_site_masked`](@ref) applies a list of `DateTime` intervals
across the runs of a site:

```julia
using Timekeepers, Dates

intervals = [
    (DateTime(2021, 9, 26, 0, 30), DateTime(2021, 9, 26, 0, 35)),
]

dest = write_metronix_site_masked("data/DF002.TK";
                                  rate_intervals = Dict(4096.0 => intervals))
# "data/DF002.TK20260930_141205"
```

Every write goes to a new `<site>.TK<date>_<time>` directory beside the site,
so no earlier write is overwritten. It is laid out like `<site>.TK`, one
directory per rate:

```text
DF002.TK20260930_141205/
├── README.md
├── 128/meas_2021-09-25_14-02-01/      untouched, copied as it was
├── 131072/meas_2021-09-25_13-55-25/   untouched, copied as it was
└── 4096/
    ├── meas_2021-09-25_14-02-01/      R001 untouched, R002's XML, the .kml
    ├── meas_2021-09-26_00-00-00/      R000 up to the cut
    └── meas_2021-09-26_00-35-01/      R000 after it
```

A run no interval touches is copied byte for byte into the `meas_*`
directory it came from, together with the `.kml` and the XMLs of jobs at its
rate that never recorded. A run an interval does touch is split: each
unmasked stretch becomes its own `meas_<start>` directory, with `.ats` files
cut from the run, the `.kml`, and a copy of the run's XML. Segments shorter
than `min_samples` are dropped, and segment starts are trimmed to whole
seconds where the sampling rate requires it.

The XML of a stretch is the run's own XML with only the fields that describe
the stretch changed: the recording's start and stop, each channel's start and
sample count, and the `.ats` file size. The file is edited as text, so the
declaration, whitespace, comments and escapes stay exactly as the ADU wrote
them. The stop time follows the ADU's convention — the start plus the whole
seconds recorded, so two hours at 4096 Hz run from `00:00:00` to `02:00:00` —
and the filename carries the same start and stop. The edit is checked against
the same change made through an XML parser, and if they disagree nothing is
written.

`rate_intervals` gives each rate its own intervals, as the app does: masks
drawn on the 4096 Hz bursts do not cut the 128 Hz run recorded through the
same hours. `intervals` applies to every rate not in `rate_intervals`. `only`
restricts the write to the runs named. The source can be the raw site, its
`.TK` copy or an earlier write; each gives the same destination name. This is
what the app's **Write** button calls when a Metronix site is loaded.

### A single loaded run

When you already have a run and a mask in memory,
[`write_metronix_site`](@ref) does the same split for that one run:

```julia
run  = read_metronix("data/RK137/meas_2025-04-01_07-00-05")
ta   = to_timearray(run)
mask = TimekeeperMask(ta)
mask_interval!(mask, DateTime(2025, 4, 1, 7, 30), DateTime(2025, 4, 1, 7, 35))

dest, dirs = write_metronix_site(run; mask = mask)
length(dirs)   # 2 — the record either side of the cut
```

By default the segments go to the run's rate directory in a new write of its
site, `<site>.TK<date>_<time>/<rate>`. With `mask = nothing` the run is written
whole. To write a single measurement directory with no splitting at all, use
[`write_metronix`](@ref).

## The write log

Each call to [`write_metronix_site_masked`](@ref) writes a `README.md` in its
destination, naming the source site, every interval cut at each rate, and the
directories each run was written to, so the destination carries its own
provenance:

```markdown
## Write session 2025-05-20 11:07:33

Amputated (masked) intervals at 4096 Hz:

| # | Start | End | Duration |
|---|-------|-----|----------|
| 1 | 2021-09-26T00:30:00 | 2021-09-26T00:35:00 | 5 minutes |
```

That log is the reason to keep the `.TK<date>_<time>` directory as the thing
you hand to a processing chain: it is reproducible from the source plus the record of what
was cut.

## Round-trip guarantee

Writing a site with no intervals reproduces the input exactly: put the rate
directories back together and every file of every `meas_*` directory comes
back byte-identical — the runs' `.ats` and `.xml` files, the `.kml`, and the
XMLs of jobs that never recorded. An unmasked write is a lossless copy, so it
is safe to route every site through the writer whether or not it needed
editing.
