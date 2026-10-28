# Masking & Cleaning

A [`TimekeeperMask`](@ref) records *where* the bad intervals are. Your samples are never deleted,
and a mask is a tiny file you can save, review, share and reapply.

## Make a mask

```julia
using Timekeepers, Dates

ta   = load_lemi424("data/LEMI090.txt")
mask = TimekeeperMask(ta)                      # everything starts good

mask_interval!(mask, DateTime(2020, 10, 4, 0, 10), DateTime(2020, 10, 4, 0, 20))
mask_interval!(mask, DateTime(2020, 10, 4, 3, 45), DateTime(2020, 10, 4, 3, 47))

masked_samples(mask)   # 722
mask.intervals         # the two spans, in order
```

[`unmask_interval!`](@ref) undoes a span and [`clear_mask!`](@ref) undoes everything. Limits can be
in either order; overlapping spans merge.

## Get clean data

| You need | Use |
|:---|:---|
| Full time axis, `NaN` in bad rows | `cleaned_timearray(ta, mask)` |
| Bad rows removed (irregular time axis) | `cleaned_timearray(ta, mask; mode = :drop)` |
| Uninterrupted pieces for FFTs | `good_segments(ta, mask; min_samples = 256)` |
| A weight per sample | `sample_weights(mask)` (1.0 good, 0.0 bad; or `good = 1, bad = 0`) |

Set `min_samples` to your FFT length. Prefer `:nan` unless the next step accepts irregular samples.

## Save and reapply

```julia
write_mask("data/LEMI090_mask.csv", mask)

ta2   = load_lemi424("data/LEMI090.txt")      # later, or on another machine
mask2 = read_mask("data/LEMI090_mask.csv", ta2)
mask2.masked == mask.masked   # true
```

The CSV holds `start,stop` times, so it also applies to a decimated copy of the same record.
Merge masks with [`combine_masks`](@ref):

```julia
final = combine_masks(human_mask, despike_mask, gps_dropout_mask)
```

## Write the result

```julia
write_cleaned("data/LEMI090_clean.csv", ta, mask)                     # CSV, timestamp + one column per channel
write_lemi424("data/LEMI090_clean.txt", cleaned_timearray(ta, mask))  # original format, NaN in bad rows
```

Metronix has no `NaN`: use [`write_metronix_site`](@ref) to cut the intervals out instead
([Metronix Sites](metronix.md)).

## From the app

Every function above accepts a [`TKApp`](@ref) instead of `(ta, mask)`:

```julia
app = run_tkapp("data/LEMI090.txt")   # mask a few intervals, then close

segments = good_segments(app; min_samples = 256)
write_mask("data/LEMI090_mask.csv", app)
write_cleaned("data/LEMI090_clean.csv", app)
```
