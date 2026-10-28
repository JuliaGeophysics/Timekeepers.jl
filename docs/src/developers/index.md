# For Developers

Use Timekeepers.jl from code, or extend it. To contribute, see
[CONTRIBUTING.md](https://github.com/JuliaGeophysics/Timekeepers.jl/blob/main/CONTRIBUTING.md).

## In five lines

```julia
using Timekeepers, Dates

run = read_timekeeper("data/LEMI090.txt")      # LEMI-424, GEOMAG-02 or Metronix
ta  = to_timearray(run)                         # a TimeSeries.TimeArray

mask = TimekeeperMask(ta)                       # bad intervals; samples are never deleted
mask_interval!(mask, DateTime(2020, 10, 4, 0, 10), DateTime(2020, 10, 4, 0, 20))

cleaned  = cleaned_timearray(ta, mask)                # NaN in masked rows
segments = good_segments(ta, mask; min_samples = 256) # contiguous good parts
weights  = sample_weights(mask)                       # 1.0 good, 0.0 bad

write_mask("data/LEMI090_mask.csv", mask)       # read back with read_mask
write_lemi424("data/LEMI090_clean.txt", cleaned)      # or write_geomag, in the source format
```

The windows return their state when you close them:

```julia
app = run_tkapp("data/LEMI090.txt")
cleaned_timearray(app)                          # every mask function accepts a TKApp

dash = run_tkdash("data/survey")
write_reference_plan("plan.txt", dash)
```

For Metronix, use [`write_metronix_site_masked`](@ref): it cuts masked intervals out instead of
writing `NaN`. Everything exported is in the [API reference](../api.md).

## Topics

- [Instrument formats](formats.md): readers, writers, channel names
- [Masking & cleaning](masking.md): make, save and apply masks
- [Metronix sites](metronix.md): rates, rate directories, cutting
- [TKDash](tkdash.md), [TKApp](tkapp.md), [TKProc](tkproc.md): the windows, and the same work in code
- [Spectral views](spectra.md): how the spectra are computed

## Code layout

| File | Contents |
|:-----|:---------|
| `src/Timekeepers.jl` | module, includes and exports |
| `src/Types.jl` | `TimekeeperChannel`, `TimekeeperRun` |
| `src/Utilities.jl` | shared helpers |
| `src/TimeArrayIO.jl` | `to_timearray`, `from_timearray` |
| `src/Masking.jl` | `TimekeeperMask` and the cleaning functions |
| `src/Spectra.jl` | Welch PSD and its workspace |
| `src/LEMI424.jl`, `src/GEOMAG.jl`, `src/MetronixATS.jl` | format readers and writers |
| `src/TimekeeperIO.jl` | `read_timekeeper`, `write_timekeeper`, format detection |
| `src/Explorer.jl` | TKApp |
| `src/Survey.jl` | survey scan, overlaps, reference plans |
| `src/Dashboard.jl` | TKDash |
| `src/Calibration.jl` | coil calibrations, dipoles, sensor orientation |
| `src/IGRF.jl` | IGRF-14 declination |
| `src/Processing.jl` | `estimate_tf` and `TransferFunction` |
| `src/SiteProc.jl` | `procset.dat` and `process_site` |
| `src/EDI.jl`, `src/ModEM.jl`, `src/Report.jl` | transfer function output |
| `src/Processor.jl` | TKProc |
| `src/Precompile.jl` | precompile workload |

## Tests

```bash
julia --project=. test/runtests.jl
```
