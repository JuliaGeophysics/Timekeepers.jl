# Metronix Sites

A Metronix ADU deployment makes a *site* directory. The site directory holds
`meas_*` measurement directories. One campaign often records at more than one
sampling rate, for example a long 128 Hz run with short 4096 Hz bursts between
its parts. This page tells how Timekeepers:

- finds the runs in a site;
- loads the runs one rate at a time;
- puts the runs of each sampling rate in a separate directory;
- cuts the masked intervals out of the runs.

!!! tip "Tutorial data"
    To try the procedures on this page, download the public Metronix tutorial
    data of BRGM: <https://github.com/BRGM/razorback-tutorial-data>. Each
    `site*` directory in it is a Metronix site.

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

A *run* is one run number (`R000`, `R001`, …) at one sampling rate. A run has
one `.ats` binary file for each channel and one `.xml` file that describes
them. A `meas_*` directory can hold one run or more. The ADU writes a long run
at a low rate and the high-rate bursts in the same directory. The run number
and the rate in the filenames connect each `.ats` file to its XML.

An XML file without `.ats` files next to it is a job that the ADU scheduled
but did not record. Timekeepers ignores this XML when it reads. It copies the
XML without change when it writes.

The reader does not need the XML. The `.ats` headers give the number of
samples, the rate, the start time and the scale. If you have only the `.ats`
files, the reader reads them and gives a warning. It takes the run number from
the filenames and the rate from the headers. The writer then writes the run
in the same way, without an XML. Some processing tools need the XML for the
sensor and calibration data. These tools do not accept such a run.

## Make an index of a site

[`metronix_site_rates`](@ref) gives the different sampling rates in a site.
[`metronix_site_runs`](@ref) gives the runs of each rate in order of start
time. It identifies each run by the path of its XML:

```julia
using Timekeepers

is_metronix_site("data/DF002")     # true
metronix_site_rates("data/DF002")  # [128.0, 4096.0, 131072.0]

runs = metronix_site_runs("data/DF002")
basename.(runs[4096.0])
# "081_2021-09-26_00-00-00_2021-09-26_02-00-00_R000_4096H.xml"
# "081_2021-09-27_00-00-00_2021-09-27_02-00-00_R001_4096H.xml"
```

The two functions read only the `.ats` headers. Thus, the index of a large
site is fast. The two functions accept a site or one `meas_*` directory.

[`read_metronix`](@ref) reads one run. To identify the run, give its XML, one
of its `.ats` files, or a `meas_*` directory that holds only that run:

```julia
run = read_metronix(runs[4096.0][2])
sampling_rate(run)  # 4096.0
```

## Load a site in the app

**Load Run…** opens one run. Go into the `meas_*` directory and select the
`.xml` of the run, or one of its `.ats` files. The app loads that run, with
its channels and its XML. It does not load the other files in the directory.

**Load Site…** accepts the site directory or one `meas_*` directory. It finds
all the runs and puts them in groups by sampling rate. High-rate runs are
large. Thus, the app holds only one rate in memory at a time. For a site with
more than one rate, the app asks which rate to import. The list gives each
rate with its number of runs. It also tells you about each XML that the app
ignored because that XML has no data. If you close that window, the app
cancels the import.

Before it reads the data, the loading window does a check. It looks for the
site already in separate rate directories `<site>.<rate>` next to the site
(refer to [Separating a site by rate](#Separating-a-site-by-rate)). If all the
rate directories are there, the window tells you and copies nothing. If not,
it makes the missing directories. The app then reads the runs from the
directory of the rate that you selected. If you select one of these rate
directories, the app uses the full site.

The runs of the selected rate show as one record in time order. Use the rate
menu at the left of the toolbar to change the rate. When you select a rate,
the app reads its runs from the disk and replaces the runs on the screen.

The masks stay when you change the rate. The rate that you leave keeps its
masked intervals. When you go back to that rate, the intervals come back.
**Write** cuts the rate on the screen and each other rate that has masks. It
uses the intervals of each rate and works in the directory of that rate. It
then shows the record after the cut.

In Julia, [`load_metronix_site`](@ref) reads one rate of a site in the same
way. It first puts the runs in separate rate directories. If the site has
more than one rate, you must give `rate`:

```julia
overview = load_metronix_site("data/DF002"; rate = 4096)
```

### Rates above 1 kHz

Up to 1 kHz, the runs of a rate use one time grid. The gaps between the runs
contain `NaN`.

Above 1 kHz, more than one sample occurs in each millisecond. The time stamps
of the app cannot show such small steps. Thus, the app joins the runs end to
end:

- The app keeps all the samples.
- Each gap between two runs is a jump in the time axis.
- The plot line has a break at each jump.

The spectra use the gap in the same way as a mask. No FFT segment goes across
the gap. For a window that covers more than one run, the spectrum is the
average of the segments from each run. The line below the plots tells how
many runs the spectrum uses.

A high-rate run uses much memory. For example, two hours at 4096 Hz is
approximately 30 million samples for each channel. Thus, each run can take a
few seconds and a few GB to load.

## Separating a site by rate

[`split_metronix_site`](@ref) copies a site into one directory for each
sampling rate. The name of each directory is `<site>.<rate>`, and it is next
to the site. Each directory is a usual site of `meas_*` directories with the
runs at that rate. Thus, a tool that reads the original site can also read
each rate directory:

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

The `.ats` files and the XML of each run go to the directory of its rate. The
XML of a job that did not record goes to the rate in its filename. The `.kml`
goes to each rate that the `meas_*` directory holds. The function copies each
file byte for byte. It does not change the site. If you put the rate
directories together again, you get the original site, file for file. A site
load does this work for you.

The rate directories are the working copies. The cuts occur there, and the
site itself never changes. Thus, the function never copies over a rate
directory that exists. It copies each directory fully under a temporary name.
It then gives the directory its correct name. Thus, a rate directory that
exists is always complete. To start a rate again from the recording, delete
its rate directory. The next load makes the directory again.

## Cut masked intervals out

Metronix has no `NaN`, and the processing tools need continuous runs. Thus,
Timekeepers applies a mask with a *cut*. It cuts each run at the masked
intervals. The good parts of the run then replace the run as separate runs in
the same `meas_*` directory.

### A full site on disk

[`write_metronix_site_masked`](@ref) cuts a list of `DateTime` intervals out
of the runs of a site. It works in the rate directories of the site:

```julia
using Timekeepers, Dates

intervals = [
    (DateTime(2021, 9, 26, 0, 30), DateTime(2021, 9, 26, 0, 35)),
]

write_metronix_site_masked("data/DF002"; rate_intervals = Dict(4096.0 => intervals))
# ["data/DF002.4096"]
```

If the rate directories are missing, the function first makes them. It
returns the directories that it cut.

- If no interval touches a run, the run stays as it is.
- If an interval touches a run, the good parts of the run replace it in its
  `meas_*` directory.
- The first part keeps the run number.
- Each subsequent part gets the next free run number at that rate in that
  directory. This number comes after all the runs and scheduled jobs that are
  there.
- The function does not rename other files, and it does not make new `meas_*`
  directories.

```text
DF002.4096/
├── README.md                      where the directory came from
├── mask.csv                       every stretch cut
└── meas_2021-09-25_14-02-01/      R000 up to the cut, R003 after it,
                                   R001 untouched, R002's XML, the .kml
```

If the sampling rate needs it, the function moves the start of a segment to a
whole second. Each `meas_*` directory keeps its `.kml`, and each part has its
own XML. This XML is the XML of the run. Only the fields that describe the
part change:

- the start and the stop of the recording;
- the start and the number of samples of each channel;
- the size of the `.ats` file;
- the `.ats` names, when the part gets a new run number.

The function edits the file as text. Thus, the declaration, the whitespace,
the comments and the escapes stay as the ADU wrote them. The stop time
follows the convention of the ADU: the start plus the whole seconds recorded.
For example, two hours at 4096 Hz go from `00:00:00` to `02:00:00`. The
filename has the same start and stop. The function compares the edit with the
same change made through an XML parser. If the two are different, it writes
nothing. Some `.ats` headers give the name of their XML, for example the
ADU-07 header. The function changes each such header to the new XML. It
writes the parts in a different location first. It replaces the run only when
the parts are complete.

Some parts are too short to keep as a run. The function ignores each such
part and gives a warning with its start and length. The ADU records the start
and stop in whole seconds. Thus, a part shorter than one second would be a
run with its stop equal to its start. To set a higher limit, give
`min_samples`, for example the shortest part that your processing can use. If
the mask removes all the runs of a `meas_*` directory, the function removes
the directory and gives a warning.

`rate_intervals` gives each rate its own intervals, as the app does. Thus,
masks on the 4096 Hz bursts do not cut the 128 Hz run that recorded in the
same hours. These are the other arguments:

- `intervals` applies to each rate that is not in `rate_intervals`.
- `only` limits the cut to the runs that you name.
- `site_dir` can be the site or one of its rate directories. A rate directory
  stands for the full site.
- `format = :default` gives the layout above. `format = :MTH5` is reserved
  for MTH5 output, which is not available yet.

If you cut the same intervals again, the function cuts nothing more. The
**Write** button of the app calls this function when a Metronix site is
loaded.

### One loaded run

If you have a run and a mask in memory, [`write_metronix_site`](@ref) does the
same cut for that run:

```julia
run  = read_metronix("data/RK137/meas_2025-04-01_07-00-05")
ta   = to_timearray(run)
mask = TimekeeperMask(ta)
mask_interval!(mask, DateTime(2025, 4, 1, 7, 30), DateTime(2025, 4, 1, 7, 35))

dest, runs = write_metronix_site(run; mask = mask)
# dest: "data/RK137.8/meas_2025-04-01_07-00-05"
length(runs)   # 2 — R000 before the cut, R001 after it
```

By default, the parts replace the run in its rate directory. If you read the
run from a rate directory, the function cuts it in that location. With
`mask = nothing`, the function writes the full run. To write one measurement
directory without a cut, use [`write_metronix`](@ref).

## What was cut: `mask.csv`

Each rate directory with a cut holds a `mask.csv` file. The file gives each
part that was removed, one row for each part. The rows of all the writes stay
in the file, in time order. The file is plain text with commas, and all
readers can load it:

```text
start_sample,end_sample,start_time,end_time
17,40,2025-04-01T07:00:08.000,2025-04-01T07:00:10.875
57,72,2025-04-01T07:00:13.000,2025-04-01T07:00:14.875
```

The samples have numbers from 1 at the start of the run, as the site recorded
it. The two ends are included. Thus, the numbers stay the same after more
cuts. The times are the times of the first and the last sample removed, to the
millisecond. The file also includes the parts that the function moved to a
whole second and the parts that it ignored because they were too short. Thus,
the file describes the data as written, not only the intervals that you
drew. A `README.md` of three lines next to the file gives the name of the
site that the directory came from.

## Round-trip guarantee

A write with no intervals changes nothing. The rate directories stay copies
of the site. If you put them together again, each file of each `meas_*`
directory is identical byte for byte. This includes the `.ats` and `.xml`
files of the runs, the `.kml`, and the XML files of the jobs that did not
record.
