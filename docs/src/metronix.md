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

The XML is not needed to read a run: sample count, rate, start time and
scaling all come from the `.ats` headers. Data shared as `.ats` files alone
reads with a warning, its run number taken from the filenames and its rate
from the headers, and is written back the same way, without an XML. Tools that
need the XML for sensor and calibration details, such as ProcMT, will not
accept such a run.

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
by sampling rate into `<site>.<rate>` directories beside it (see
[Separating a site by rate](#Separating-a-site-by-rate)). If every rate
directory is there it says so and copies nothing; otherwise it makes the
missing ones. The runs are then read from the directory of the rate
chosen. Picking one of those rate directories stands for the whole site.

The runs of the chosen rate appear as one record in time order. The rate menu
at the left of the toolbar then switches rates: picking one reads its runs from
disk in place of those on screen. Masks survive the switch — the rate you leave
keeps its masked intervals, and gets them back when you return — and **Write**
cuts the rate on screen and every other rate you masked by their own
intervals, each in its rate directory, then shows the cut record.

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

[`split_metronix_site`](@ref) copies a site into one directory per sampling
rate, named `<site>.<rate>` and placed beside it. Each is an ordinary site of
`meas_*` directories holding the runs at that rate, so anything that reads
the original site reads each of them:

```text
DF002/                                 the site as recorded, never changed
DF002.128/
└── meas_2021-09-25_14-02-01/          R000 at 128 Hz, its XML, the .kml
DF002.4096/
└── meas_2021-09-25_14-02-01/          R000, R001 at 4096 Hz, their XMLs,
                                       R002's XML (never recorded), the .kml
DF002.131072/
└── meas_2021-09-25_13-55-25/          R000 at 131072 Hz, its XML, the .kml
```

Each run's `.ats` files and XML go to its rate. The XML of a job that never
recorded goes to the rate in its filename, and the `.kml` to every rate the
`meas_*` directory holds. Every file is copied byte for byte and the site is
left alone; putting the rate directories back together gives the original
site file for file. Loading a site does this for you.

The rate directories are the working copies: masks are cut there, and the
site itself is never changed. So a rate directory that exists is never copied
over. Each is copied in full under a temporary name and renamed when done, so
one that exists is complete. To start a rate again from the recording, delete
its rate directory; the next load makes it afresh.

## Amputating masked intervals

Metronix has no `NaN`, and processing tools expect continuous runs. So a mask
is applied by *cutting*: each run is split at the masked intervals, and its
surviving stretches replace it as runs of their own in its `meas_*`
directory.

### A whole site on disk

[`write_metronix_site_masked`](@ref) cuts a list of `DateTime` intervals out
of the runs of a site, in its rate directories:

```julia
using Timekeepers, Dates

intervals = [
    (DateTime(2021, 9, 26, 0, 30), DateTime(2021, 9, 26, 0, 35)),
]

write_metronix_site_masked("data/DF002"; rate_intervals = Dict(4096.0 => intervals))
# ["data/DF002.4096"]
```

The rate directories are made first if missing, and the call returns those it
cut. A run no interval touches is left as it is. A run an interval does touch
is replaced, in its `meas_*` directory, by its unmasked stretches: the first
keeps the run number, and each later one takes the next run number free at
that rate in that directory, after every run and scheduled job already there.
No other file is renamed and no new `meas_*` directory appears:

```text
DF002.4096/
├── README.md                      where the directory came from
├── mask.csv                       every stretch cut
└── meas_2021-09-25_14-02-01/      R000 up to the cut, R003 after it,
                                   R001 untouched, R002's XML, the .kml
```

Segment starts are trimmed to whole seconds where the sampling rate requires
it. Every `meas_*` directory keeps its `.kml`, and every stretch has its own
XML. That XML is the run's own with only the fields that describe the stretch
changed: the recording's start and stop, each channel's start and sample
count, the `.ats` file size, and the `.ats` names when the stretch is
renumbered. The file is edited as text, so the declaration, whitespace,
comments and escapes stay exactly as the ADU wrote them. The stop time follows
the ADU's convention — the start plus the whole seconds recorded, so two
hours at 4096 Hz run from `00:00:00` to `02:00:00` — and the filename carries
the same start and stop. The edit is checked against the same change made
through an XML parser, and if they disagree nothing is written. Each `.ats`
header that names its XML, as an ADU-07 header does, is pointed at the new
one. The stretches are written aside first and replace the run only once they
are complete.

A stretch too short to be stored as a run is skipped with a warning naming
its start and length. The ADU records start and stop in whole seconds, so
anything under one second would be a run whose stop equals its start; pass
`min_samples` to raise the bar, for instance to the shortest stretch your
processing can use. A `meas_*` directory whose runs are all masked away is
removed, with a warning.

`rate_intervals` gives each rate its own intervals, as the app does: masks
drawn on the 4096 Hz bursts do not cut the 128 Hz run recorded through the
same hours. `intervals` applies to every rate not in `rate_intervals`. `only`
restricts the cut to the runs named. `site_dir` can be the site or one of its
rate directories, which stands for the whole site. Cutting the same intervals
again cuts nothing more. `format = :default` is the layout above;
`format = :MTH5` is reserved for MTH5 output, not available yet. This is what
the app's **Write** button calls when a Metronix site is loaded.

### A single loaded run

When you already have a run and a mask in memory,
[`write_metronix_site`](@ref) does the same cut for that one run:

```julia
run  = read_metronix("data/RK137/meas_2025-04-01_07-00-05")
ta   = to_timearray(run)
mask = TimekeeperMask(ta)
mask_interval!(mask, DateTime(2025, 4, 1, 7, 30), DateTime(2025, 4, 1, 7, 35))

dest, runs = write_metronix_site(run; mask = mask)
# dest: "data/RK137.8/meas_2025-04-01_07-00-05"
length(runs)   # 2 — R000 before the cut, R001 after it
```

By default the stretches replace the run in its rate directory; a run read
from a rate directory is cut where it lies. With `mask = nothing` the run is
written whole. To write a single measurement directory with no cutting at
all, use [`write_metronix`](@ref).

## What was cut: `mask.csv`

Each rate directory cut holds a `mask.csv` listing every stretch removed, one
row each, accumulated over writes and sorted by time — plain comma-separated
text any reader can load:

```text
start_sample,end_sample,start_time,end_time
17,40,2025-04-01T07:00:08.000,2025-04-01T07:00:10.875
57,72,2025-04-01T07:00:13.000,2025-04-01T07:00:14.875
```

Samples are counted from 1 at the start of the run as recorded in the site,
both ends included, so they stay the same however often a run is cut; the
times are those of the first and last sample removed, to the millisecond.
Stretches trimmed to a whole second or skipped as too short are included, so
the file describes the data as written, not only the intervals drawn. A
three-line `README.md` beside it names the site the directory came from.

## Round-trip guarantee

A write with no intervals changes nothing: the rate directories stay copies of
the site, and putting them back together gives every file of every `meas_*`
directory byte-identical — the runs' `.ats` and `.xml` files, the `.kml`, and
the XMLs of jobs that never recorded.
