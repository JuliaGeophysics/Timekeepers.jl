# Instrument Formats

Timekeepers reads three instrument formats natively. It also writes the three
formats. Each format has the same three functions:

| Format | Run reader | `TimeArray` reader | Writer |
|:---|:---|:---|:---|
| LEMI-424 | [`read_lemi424`](@ref) | [`load_lemi424`](@ref) | [`write_lemi424`](@ref) |
| GEOMAG-02 | [`read_geomag`](@ref) | [`load_geomag`](@ref) | [`write_geomag`](@ref) |
| Metronix ADU | [`read_metronix`](@ref) | [`load_metronix`](@ref) | [`write_metronix`](@ref) |

[`read_timekeeper`](@ref) and [`write_timekeeper`](@ref) operate for the three
formats. They select the correct function automatically.

## Component naming

All the formats use one set of channel names. Thus, the other parts of the
package (masks, plots and spectra) operate the same for each format:

| Component | Meaning | Units |
|:---|:---|:---|
| `:bx`, `:by`, `:bz` | magnetic field, geographic N / E / down | nT |
| `:e1`, `:e2` | electric field, dipole 1 / 2 | mV/km |
| `:temperature_e`, `:temperature_h` | electronics / sensor temperature | °C |
| `:battery` | supply voltage | V |
| `:elevation` | GPS elevation | m |

[`default_components`](@ref) returns the five signal channels in the usual
plot order (`bx, by, bz, e1, e2`). [`components`](@ref) returns all the
channels in the run, which include the auxiliary channels.

## LEMI-424

LEMI-424 is a long-period ASCII text format. Each sample has one record, with
the fields separated by whitespace. A record has six date and time fields.
Then it has the magnetic, electric, temperature, GPS and housekeeping columns.
[`LEMI424_COLUMNS`](@ref) gives the names of the 24 fields in the file order.

```julia
using Timekeepers

run = read_lemi424("data/LEMI090.txt")
ta  = load_lemi424("data/LEMI090.txt")     # straight to a TimeArray
```

Some real files do not have exactly 24 columns. Some firmware revisions add
fields at the end, and some exports remove fields at the end. The reader
accepts the two cases. It ignores the extra columns. It returns `NaN` for each
missing column at the end, and does not give an error.

```julia
run = read_lemi424("data/short_record.txt")
all(isnan, run.channels[:time_diff].data)   # true when the column was absent
```

If you only need the five signal components, give `include_aux = false`. The
reader then does not read the housekeeping channels.

The writer makes the original layout of 24 fields again. The `TimeArray`
metadata keeps the auxiliary values under `:aux_columns`. The writer puts each
value back in its own position. Thus, a cycle of read and write does not lose
data:

```julia
ta = load_lemi424("data/LEMI090.txt")
write_lemi424("data/LEMI090_copy.txt", ta)
```

The two writers accept a [`TimekeeperRun`](@ref) or a `TimeArray`.

### Generic LEMI `.xyz` exports

TKApp also loads and writes LEMI-style exports that have 7 columns
(`date time Bx By Bz Ex Ey`). The file must have the `.xyz` extension. There
is no public reader for these files. Open them with the app or with
[`run_tkapp`](@ref).

## GEOMAG-02

GEOMAG-02 is an ASCII text format. The file starts with a block of header
lines. Each header line starts with `;`. The header gives the instrument
model, the sampling interval and the position of the station. Then the file
has one record for each sample, with fractional seconds.

```julia
run = read_geomag("data/MS_26_250523000000.TXT")
sampling_rate(run)                     # 10.0 for a 0.10 s interval
run.metadata[:instrument_model]        # "GEOMAG-02"
```

The format detection looks for the `GEOMAG` token in the header. Thus,
[`read_timekeeper`](@ref) can tell a GEOMAG file from a LEMI-424 `.txt` file
without help.

The reader reads the five signal channels. It also reads `:temperature_h` and
`:temperature_e`. [`write_geomag`](@ref) writes the header block and the data.
Thus, the tools that read the input can also read the output.

## Metronix ADU (ATS)

A Metronix run is a set of files in a `meas_*` directory. Each channel has one
`.ats` binary file, and the run has an `.xml` sidecar. A `meas_*` directory
can hold more than one run, for example a 128 Hz run and 4096 Hz bursts. The
run number (`R000`, `R001`, …) and the rate in the filenames identify each
run.

The `.ats` header gives these values:

- the number of samples;
- the sampling rate;
- the start time, as a Unix timestamp;
- the LSB scale value;
- the channel type.

The samples are `Int32` counts. The reader multiplies them by the LSB value.

```julia
run = read_metronix("data/RK137/meas_2025-04-01_07-00-05")
components(run)     # [:bx, :by, :bz, :e1, :e2]
sampling_rate(run)  # 128.0
```

[`METRONIX_CHANNEL_MAP`](@ref) changes the channel types to component names:
`Ex → :e1`, `Ey → :e2`, `Hx → :bx`, `Hy → :by`, `Hz → :bz`. The reader does
not remove a channel type that it does not know. It keeps that channel under
a symbol with the same name.

!!! note "Point at one run"
    [`read_metronix`](@ref) accepts the `.xml` of a run, one of its `.ats`
    files, or a `meas_*` directory that holds one run. A directory that holds
    more than one run is ambiguous, and the error gives a list of the runs.
    For the site directory above the `meas_*` directories, use
    [`metronix_site_runs`](@ref) and the site writers. Refer to
    [Metronix Sites](metronix.md).

[`write_metronix`](@ref) writes a measurement directory. It makes the XML
sidecar again from the template that the reader used. The path of that
template is in `run.metadata[:metronix_xml_path]`. Thus, the run must come
from [`read_metronix`](@ref) before you can write it as Metronix.

A mask on Metronix data usually means "cut this interval out". Thus, Metronix
has a second writer that cuts the data and does not put blanks in it. Refer
to [Metronix Sites](metronix.md).

## Metadata that the readers keep

The readers keep the run information that a writer needs. You can find all of
it in `run.metadata`, or in `TimeSeries.meta(ta)` for a `TimeArray`:

| Key | Present for | Meaning |
|:---|:---|:---|
| `:sample_rate` | all | samples per second |
| `:start_time` | all | timestamp of the first sample |
| `:instrument_model` | GEOMAG, Metronix | instrument identification |
| `:aux_columns` | LEMI-424, GEOMAG | values of the auxiliary columns, for writes that do not lose data |
| `:site`, `:n_files` | site loads | site name and number of runs that were joined |
| `:metronix_xml_path` | Metronix | XML template that makes the sidecar again |
| `:metronix_prefix`, `:metronix_run_token`, `:metronix_freq_token` | Metronix | filename tokens that the writer keeps |
| `:mask_intervals`, `:masked_samples` | after a mask | the intervals that were cut, and the number of samples |
