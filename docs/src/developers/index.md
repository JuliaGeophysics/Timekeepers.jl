# For Developers

This section is for people who want to use Timekeepers.jl from code, or change or extend it. For how to contribute, see
[CONTRIBUTING.md](https://github.com/JuliaGeophysics/Timekeepers.jl/blob/main/CONTRIBUTING.md).

## From code

The windows are built on functions that you can call yourself:

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

`run_tkapp` and `run_tkdash` return the app when the window closes, so its mask or plan stays usable:

```julia
app = run_tkapp("data/LEMI090.txt")
cleaned_timearray(app)                          # every mask function accepts a TKApp

dash = run_tkdash("data/survey")
write_reference_plan("plan.txt", dash)
```

For Metronix, use [`write_metronix_site_masked`](@ref), which cuts the masked intervals out instead
of writing `NaN`. All exported functions are in the [API reference](../api.md).

## In detail

- [Instrument formats](formats.md): file layouts, component names, metadata the readers keep
- [Masking & cleaning](masking.md): mask semantics, outputs, combining masks
- [Metronix sites](metronix.md): site anatomy, rate directories, cutting, `mask.csv`
- [TKDash](tkdash.md): survey scan, window controls, plan file
- [TKApp](tkapp.md): every toolbar control, Metronix loading, write behaviour
- [Spectral views](spectra.md): the Welch estimate, transform length, performance

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
| `src/Precompile.jl` | precompile workload |

## Tests

```bash
julia --project=. test/runtests.jl
```
