# Instrument Formats

| Format | Read a run | Read a `TimeArray` | Write |
|:---|:---|:---|:---|
| LEMI-424 | [`read_lemi424`](@ref) | [`load_lemi424`](@ref) | [`write_lemi424`](@ref) |
| GEOMAG-02 | [`read_geomag`](@ref) | [`load_geomag`](@ref) | [`write_geomag`](@ref) |
| Metronix ADU | [`read_metronix`](@ref) | [`load_metronix`](@ref) | [`write_metronix`](@ref) |

Don't know the format? [`read_timekeeper`](@ref) and [`write_timekeeper`](@ref) pick it for you.

## Channel names

Every format uses the same names:

| Component | Meaning | Units |
|:---|:---|:---|
| `:bx`, `:by`, `:bz` | magnetic field, geographic N / E / down | nT |
| `:e1`, `:e2` | electric field, dipole 1 / 2 | mV/km |
| `:temperature_e`, `:temperature_h` | electronics / sensor temperature | °C |
| `:battery` | supply voltage | V |
| `:elevation` | GPS elevation | m |

[`default_components`](@ref) gives the five signal channels (`bx, by, bz, e1, e2`);
[`components`](@ref) gives every channel in the run.

## LEMI-424

```julia
using Timekeepers

run = read_lemi424("data/LEMI090.txt")
ta  = load_lemi424("data/LEMI090.txt")     # straight to a TimeArray
write_lemi424("data/LEMI090_copy.txt", ta) # same 24-field layout, nothing lost
```

- Whitespace-separated ASCII, one line per sample; [`LEMI424_COLUMNS`](@ref) names the 24 fields.
- Extra columns are ignored; missing trailing columns come back as `NaN`.
- `include_aux = false` reads only the five signal channels.
- 7-column `.xyz` exports (`date time Bx By Bz Ex Ey`) open in [`run_tkapp`](@ref) only.

## GEOMAG-02

```julia
run = read_geomag("data/MS_26_250523000000.TXT")
sampling_rate(run)                     # 10.0 for a 0.10 s interval
run.metadata[:instrument_model]        # "GEOMAG-02"
```

- ASCII with a `;` header (model, interval, position), detected by the `GEOMAG` token.
- Reads the five signal channels plus both temperatures; [`write_geomag`](@ref) writes header and
  data back.

## Metronix ADU (ATS)

```julia
run = read_metronix("data/RK137/meas_2025-04-01_07-00-05")
components(run)     # [:bx, :by, :bz, :e1, :e2]
sampling_rate(run)  # 128.0
```

- A run is one `.ats` file per channel plus an `.xml`, inside a `meas_*` directory.
- Pass the run's `.xml`, one of its `.ats` files, or a `meas_*` directory that holds only that
  run. A directory with several runs gives an error listing them.
- [`METRONIX_CHANNEL_MAP`](@ref): `Ex → :e1`, `Ey → :e2`, `Hx → :bx`, `Hy → :by`, `Hz → :bz`;
  unknown types keep their own name.
- [`write_metronix`](@ref) needs a run read by [`read_metronix`](@ref) (it reuses its XML).
- To cut masked intervals out of a whole site, see [Metronix Sites](metronix.md).

## Metadata

In `run.metadata`, or `TimeSeries.meta(ta)`:

| Key | Present for | Meaning |
|:---|:---|:---|
| `:sample_rate` | all | samples per second |
| `:start_time` | all | timestamp of the first sample |
| `:instrument_model` | GEOMAG, Metronix | instrument identification |
| `:aux_columns` | LEMI-424, GEOMAG | auxiliary values, for lossless writes |
| `:site`, `:n_files` | site loads | site name and number of runs joined |
| `:metronix_xml_path` | Metronix | XML template for the writer |
| `:metronix_prefix`, `:metronix_run_token`, `:metronix_freq_token` | Metronix | filename tokens the writer keeps |
| `:mask_intervals`, `:masked_samples` | after a mask | intervals cut and number of samples |
