```@raw html
---
layout: home

hero:
  name: "Timekeepers.jl"
  tagline: Timeseries analysis in Electromagnetic Geophysics
  actions:
    - theme: brand
      text: Getting Started
      link: /getting_started
    - theme: alt
      text: TKApp Explorer
      link: /tkapp
    - theme: alt
      text: View on GitHub
      link: https://github.com/JuliaGeophysics/Timekeepers.jl
---
```

![TKApp showing a five-channel Metronix record](assets/TK.png)

## What is Timekeepers.jl?

Timekeepers.jl reads recordings in the native format of the logger. It
converts them to [TimeSeries.jl](https://github.com/JuliaStats/TimeSeries.jl)
`TimeArray`s. For each format that it reads, it also has a writer in that
format.

Timekeepers.jl also has **TKApp**, a GLMakie window. Use TKApp to:

- scroll through long records;
- mark bad intervals;
- write the result in the format of the instrument.

Timekeepers.jl is part of the
[JuliaGeophysics ecosystem](https://github.com/JuliaGeophysics). Use it at the
start of a processing chain:

1. Show the raw record on the screen.
2. Remove the noise from the record.
3. Send clean segments to the next step.

## Features

- **Read and write three instrument formats.** Timekeepers reads and writes
  LEMI-424 long-period ASCII, GEOMAG-02 ASCII and Metronix ADU (ATS binary
  with an XML sidecar). The auxiliary columns stay in the data. Thus, the
  acquisition software can read a file that Timekeepers writes. Refer to
  [Instrument Formats](formats.md).
- **Load a full site.** Select a directory, and Timekeepers reads all the runs
  in it. It puts the runs in time order, joins them and fills the gaps. Thus,
  a month of hourly files becomes one continuous series.
- **Use a mask, not a destructive edit.** A [`TimekeeperMask`](@ref) records
  the bad intervals next to the data. From the mask, you can make a series
  with `NaN` in the bad intervals, the contiguous good segments, or a weight
  for each sample. You can save the mask as a small CSV file and use it again
  later. Refer to [Masking & Cleaning](masking.md).
- **Examine the data interactively.** [`run_tkapp`](@ref) opens a native
  window. The window has a time-series view that scrolls, masks that you make
  with the mouse, and optional PSD panels for each channel. Refer to
  [TKApp Explorer](tkapp.md) and [Spectral Views](spectra.md).
- **Edit a Metronix site.** Timekeepers finds all the runs of a site that has
  more than one sampling rate. It loads the site one rate at a time. It then
  cuts the masked intervals out and writes each clean segment in its own
  `meas_*` directory, with a record of each cut. Refer to
  [Metronix Sites](metronix.md).
- **Select base and remote sites for a survey.** [`run_tkdash`](@ref) opens
  TKDash. TKDash shows when each site recorded. For each site, it finds the
  base sites and the remote sites that recorded at the same time. Refer to
  [TKDash Survey](tkdash.md).

## Installation

```julia
pkg> add https://github.com/JuliaGeophysics/Timekeepers.jl
```

Timekeepers.jl needs Julia 1.12 or newer. It is not in the General registry
yet.

GLMakie is a hard dependency. Thus, TKApp needs a desktop session with
OpenGL 3.3 or newer drivers. Refer to
[Getting Started](getting_started.md#Checking-your-OpenGL-setup) for a smoke
test.

## Quick start

To open the explorer on a file, do this command:

```julia
using Timekeepers
run_tkapp("data/LEMI090.txt")
```

You can also do the same work in code:

```julia
using Timekeepers, Dates

# Read a run and convert to a TimeArray
run = read_timekeeper("data/LEMI090.txt")
ta  = to_timearray(run)

# Mark a bad interval
mask = TimekeeperMask(ta)
mask_interval!(mask, DateTime(2020, 10, 4, 0, 10), DateTime(2020, 10, 4, 0, 20))

# Derive processing-ready outputs
cleaned  = cleaned_timearray(ta, mask)                # NaN in masked rows
segments = good_segments(ta, mask; min_samples = 256) # contiguous good chunks
weights  = sample_weights(mask)                       # for robust processing

# Write back out in the source format, plus the mask itself
write_lemi424("data/LEMI090_clean.txt", cleaned)
write_mask("data/LEMI090_mask.csv", mask)
```

## Getting test data

The package does not include sample recordings. For a quick test, download a
public LEMI-424 dataset from the British Geological Survey accession. Then
extract a `.txt` file into your working directory:

> <https://webapps.bgs.ac.uk/services/ngdc/accessions/index.html#item182849>

## Citing

If you use Timekeepers.jl in published work, cite the repository:
<https://github.com/JuliaGeophysics/Timekeepers.jl>.
