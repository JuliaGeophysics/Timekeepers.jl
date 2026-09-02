# Instrument Formats

Timekeepers reads three instrument formats natively and writes all three back
out. Each has the same three-function shape:

| Format | Run reader | `TimeArray` reader | Writer |
|:---|:---|:---|:---|
| LEMI-424 | [`read_lemi424`](@ref) | [`load_lemi424`](@ref) | [`write_lemi424`](@ref) |
| GEOMAG-02 | [`read_geomag`](@ref) | [`load_geomag`](@ref) | [`write_geomag`](@ref) |
| Metronix ADU | [`read_metronix`](@ref) | [`load_metronix`](@ref) | [`write_metronix`](@ref) |

[`read_timekeeper`](@ref) and [`write_timekeeper`](@ref) sit in front of all
three and pick the right one automatically.

## Component naming

Whatever the source format, channels arrive under one set of names, so the rest
of the package — masking, plotting, spectra — does not care where the data came
from:

| Component | Meaning | Units |
|:---|:---|:---|
| `:bx`, `:by`, `:bz` | magnetic field, geographic N / E / down | nT |
| `:e1`, `:e2` | electric field, dipole 1 / 2 | mV/km |
| `:temperature_e`, `:temperature_h` | electronics / sensor temperature | °C |
| `:battery` | supply voltage | V |
| `:elevation` | GPS elevation | m |

[`default_components`](@ref) returns the five signal channels in conventional
plotting order (`bx, by, bz, e1, e2`); [`components`](@ref) returns everything
present, auxiliary channels included.

## LEMI-424

Long-period ASCII text: one whitespace-separated record per sample, six
date/time fields followed by the magnetic, electric, temperature, GPS and
housekeeping columns. [`LEMI424_COLUMNS`](@ref) names all 24 fields in file
order.

```julia
using Timekeepers

run = read_lemi424("data/LEMI090.txt")
ta  = load_lemi424("data/LEMI090.txt")     # straight to a TimeArray
```

Real files are not always exactly 24 columns wide — firmware revisions add
trailing fields, and some exports truncate them. The reader tolerates both:
extra columns are ignored, and missing trailing columns come back as `NaN`
rather than an error.

```julia
run = read_lemi424("data/short_record.txt")
all(isnan, run.channels[:time_diff].data)   # true when the column was absent
```

Pass `include_aux = false` to skip the housekeeping channels when you only want
the five signal components.

The writer reproduces the original 24-field layout. Auxiliary values are held
in the `TimeArray` metadata under `:aux_columns` and written back into their
own slots, so a read/write cycle is lossless:

```julia
ta = load_lemi424("data/LEMI090.txt")
write_lemi424("data/LEMI090_copy.txt", ta)
```

Both writers accept either a [`TimekeeperRun`](@ref) or a `TimeArray`.

### Generic LEMI `.xyz` exports

TKApp additionally loads and writes 7-column LEMI-style exports
(`date time Bx By Bz Ex Ey`) when a file has the `.xyz` extension. There is no
public reader for these; open them through the app or through
[`run_tkapp`](@ref).

## GEOMAG-02

ASCII text opening with a block of `;`-prefixed header lines carrying the
instrument model, sampling interval and station position, then one record per
sample with fractional seconds.

```julia
run = read_geomag("data/MS_26_250523000000.TXT")
sampling_rate(run)                     # 10.0 for a 0.10 s interval
run.metadata[:instrument_model]        # "GEOMAG-02"
```

Format detection keys on the `GEOMAG` token in the header, so
[`read_timekeeper`](@ref) tells these apart from LEMI-424 `.txt` files without
being told.

Beyond the five signal channels the reader picks up `:temperature_h` and
`:temperature_e`, and [`write_geomag`](@ref) emits the header block along with
the data, so the output is readable by the same tools as the input.

## Metronix ADU (ATS)

A Metronix measurement is a *directory*, not a file: one `.ats` binary per
channel plus a shared `.xml` sidecar. The `.ats` header carries the sample
count, sampling rate, start time as a Unix timestamp, the LSB scaling value and
the channel type; the samples themselves are `Int32` counts that the reader
scales by the LSB.

```julia
run = read_metronix("data/RK137/meas_2025-04-01_07-00-05")
components(run)     # [:bx, :by, :bz, :e1, :e2]
sampling_rate(run)  # 128.0
```

Channel types map to component names through
[`METRONIX_CHANNEL_MAP`](@ref) — `Ex → :e1`, `Ey → :e2`, `Hx → :bx`,
`Hy → :by`, `Hz → :bz`. Unrecognised channel types are carried through under a
symbol of their own name rather than dropped.

!!! note "Point at the `meas_*` directory"
    [`read_metronix`](@ref) takes a single measurement directory. The directory
    *above* it — the one holding several `meas_*` runs — is the site directory,
    and is handled by [`metronix_site_runs`](@ref) and the site writers covered
    in [Metronix Sites](metronix.md).

[`write_metronix`](@ref) writes a measurement directory back out, regenerating
the XML sidecar from the template the run was read with. That template path
lives in `run.metadata[:metronix_xml_path]`, so a run must have come from
[`read_metronix`](@ref) to be written as Metronix.

Because a mask usually means "cut this out", Metronix has a second writer that
splits rather than blanks — see [Metronix Sites](metronix.md).

## What metadata survives

Readers keep the run-level information a writer needs, and it is all reachable
from `run.metadata` (or `TimeSeries.meta(ta)` on the `TimeArray` side):

| Key | Present for | Meaning |
|:---|:---|:---|
| `:sample_rate` | all | samples per second |
| `:start_time` | all | timestamp of the first sample |
| `:instrument_model` | GEOMAG, Metronix | instrument identification |
| `:aux_columns` | LEMI-424, GEOMAG | auxiliary column values, for lossless writes |
| `:site`, `:n_files` | site loads | site name and number of runs combined |
| `:metronix_xml_path` | Metronix | XML template used to regenerate the sidecar |
| `:metronix_prefix`, `:metronix_run_token`, `:metronix_freq_token` | Metronix | filename tokens preserved on write |
| `:mask_intervals`, `:masked_samples` | after masking | what was cut, and how much |
