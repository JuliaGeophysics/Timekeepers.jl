# Metronix Sites

A Metronix site is a directory of `meas_*` directories, often at several sampling rates (say, a
long 128 Hz run plus 4096 Hz bursts). Timekeepers loads one rate at a time and cuts masks out
of copies, so **your recorded site is never changed**.

!!! tip "Tutorial data"
    Each `site*` directory of the public Metronix tutorial data of BRGM (Alsace, IMAGE FP7 project)
    is a Metronix site.

## Quick start

```julia
using Timekeepers, Dates

metronix_site_rates("data/DF002")                    # [128.0, 4096.0, 131072.0]
ta = load_metronix_site("data/DF002"; rate = 4096)   # one rate, all its runs

intervals = [(DateTime(2021, 9, 26, 0, 30), DateTime(2021, 9, 26, 0, 35))]
write_metronix_site_masked("data/DF002"; rate_intervals = Dict(4096.0 => intervals))
# ["data/DF002.4096"]
```

Or open the site in TKApp with **Load Site…**, mask, and press **Write**.

## What's in a site

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

- A *run* is one run number (`R000`, `R001`, …) at one rate: one `.ats` per channel plus an `.xml`.
- An `.xml` with no `.ats` files is a scheduled job that never recorded. It's ignored on read and
  copied as-is on write.
- Without an `.xml`, the run still reads (with a warning), but tools that need sensor data from the
  XML won't accept it.

## Index a site

```julia
is_metronix_site("data/DF002")     # true
runs = metronix_site_runs("data/DF002")
basename.(runs[4096.0])
# "081_2021-09-26_00-00-00_2021-09-26_02-00-00_R000_4096H.xml"
# "081_2021-09-27_00-00-00_2021-09-27_02-00-00_R001_4096H.xml"

run = read_metronix(runs[4096.0][2])
sampling_rate(run)  # 4096.0
```

These read only the `.ats` headers, so they're fast. They accept a site or one `meas_*` directory.

## In the app

- **Load Run…**: go into `meas_*` and pick a run's `.xml` or `.ats`.
- **Load Site…**: pick the site (or a `meas_*`). With several rates, the app asks which one.
  Switch later with the rate menu at the left of the toolbar.
- Each rate keeps its own masks. **Write** cuts every rate you masked.
- Above 1 kHz, runs are joined end to end (each gap is a jump in the time axis) and no FFT segment
  crosses a gap. Up to 1 kHz, gaps are `NaN`.
- High rates are heavy: two hours at 4096 Hz is about 30 million samples per channel.

## Separating a site by rate

[`split_metronix_site`](@ref) copies the site into one directory per rate, next to it. The first
masked write does this for you.

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

Files are copied byte for byte, and an existing rate directory is never overwritten. To start a rate
over, delete its directory; the next masked write makes it again.

## Cut masked intervals out

[`write_metronix_site_masked`](@ref) (what **Write** calls) cuts intervals out of the runs in the
rate directories and returns the directories it cut:

- Untouched runs stay as they are. A cut run is replaced by its good parts: the first keeps its run
  number, the rest get the next free numbers.
- Each part gets its own XML, edited as text so the ADU's formatting is kept.
- Parts shorter than one second (or `min_samples`) are skipped with a warning.
- Cutting the same intervals again does nothing.

```text
DF002.4096/
├── README.md                      where the directory came from
├── mask.csv                       every stretch cut
└── meas_2021-09-25_14-02-01/      R000 up to the cut, R003 after it,
                                   R001 untouched, R002's XML, the .kml
```

| Argument | |
|:---|:---|
| `rate_intervals` | intervals per rate |
| `intervals` | intervals for every rate not in `rate_intervals` |
| `only` | cut only the runs you name |
| `site_dir` | the site or one of its rate directories |
| `format` | `:default`; `:MTH5` is reserved, not available yet |

For one run and a mask in memory, use [`write_metronix_site`](@ref):

```julia
run  = read_metronix("data/RK137/meas_2025-04-01_07-00-05")
ta   = to_timearray(run)
mask = TimekeeperMask(ta)
mask_interval!(mask, DateTime(2025, 4, 1, 7, 30), DateTime(2025, 4, 1, 7, 35))

dest, runs = write_metronix_site(run; mask = mask)
# dest: "data/RK137.8/meas_2025-04-01_07-00-05"
length(runs)   # 2 — R000 before the cut, R001 after it
```

With `mask = nothing` it writes the whole run. For one measurement directory without cuts, use
[`write_metronix`](@ref).

## `mask.csv`

Every rate directory with a cut lists what was removed:

```text
start_sample,end_sample,start_time,end_time
17,40,2025-04-01T07:00:08.000,2025-04-01T07:00:10.875
57,72,2025-04-01T07:00:13.000,2025-04-01T07:00:14.875
```

Sample numbers count from 1 in the run as recorded, both ends included. The file describes the data
as written, including parts moved to a whole second and parts that were too short.

A write with no intervals changes nothing: the rate directories stay byte-for-byte copies of the
site.
