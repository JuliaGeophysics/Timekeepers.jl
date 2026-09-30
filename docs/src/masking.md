# Masking & Cleaning

Field recordings contain intervals you do not want: a technician walking past
the sensor, a nearby vehicle, a GPS dropout, a battery swap. Timekeepers does
not delete those samples. It records *where* they are in a
[`TimekeeperMask`](@ref) that travels alongside the data, and derives whatever
shape the next processing step needs from the pair.

That separation matters in practice. The mask is a few kilobytes of intervals,
so it can be version-controlled, reviewed, replayed onto a re-read of the same
file, or merged with a colleague's — none of which is possible once the samples
are gone.

## Creating and editing a mask

A mask is built over a `TimeArray`'s time axis and starts out all-good:

```julia
using Timekeepers, Dates

ta   = load_lemi424("data/LEMI090.txt")
mask = TimekeeperMask(ta)

mask_interval!(mask, DateTime(2020, 10, 4, 0, 10), DateTime(2020, 10, 4, 0, 20))
mask_interval!(mask, DateTime(2020, 10, 4, 3, 45), DateTime(2020, 10, 4, 3, 47))

masked_samples(mask)   # 722
mask.intervals         # the two spans, in order
```

[`mask_interval!`](@ref) marks a closed interval bad,
[`unmask_interval!`](@ref) undoes it, and [`clear_mask!`](@ref) resets
everything. Bounds may be given in either order. The derived `intervals` list
is recomputed after each edit, so it always reflects the current flags —
adjacent or overlapping edits merge into single spans automatically.

## Deriving clean data

Three functions turn a mask plus its data into something a processing step can
consume.

### `NaN`-filled series

[`cleaned_timearray`](@ref) keeps the time axis intact and writes `NaN` into
masked rows. Sample spacing stays uniform, and the gaps stay visible when you
plot:

```julia
cleaned = cleaned_timearray(ta, mask)          # mode = :nan, the default
```

Pass `mode = :drop` to remove the rows instead. That gives a shorter series but
a non-uniform time axis, which most spectral methods will not accept — prefer
`:nan` unless you know the consumer handles irregular sampling.

### Contiguous good segments

[`good_segments`](@ref) splits the record at the masked intervals and returns
the surviving stretches. This is what you want when the next step needs
uninterrupted windows:

```julia
segments = good_segments(ta, mask; min_samples = 256)
length(segments)                    # 3
```

`min_samples` discards fragments too short to be useful — set it to the FFT
length you intend to use.

### Per-sample weights

[`sample_weights`](@ref) produces a weight vector for methods that would rather
downweight bad samples than drop them:

```julia
w = sample_weights(mask)                     # 1.0 good, 0.0 bad
w = sample_weights(mask; good = 1, bad = 0)  # integer weights
```

## Saving and replaying a mask

[`write_mask`](@ref) writes the intervals to a two-column `start,stop` CSV, and
[`read_mask`](@ref) rebuilds a mask from that file over any time axis that
covers the same period:

```julia
write_mask("data/LEMI090_mask.csv", mask)

# later, or on another machine
ta2   = load_lemi424("data/LEMI090.txt")
mask2 = read_mask("data/LEMI090_mask.csv", ta2)
mask2.masked == mask.masked   # true
```

Because [`read_mask`](@ref) replays through [`mask_interval!`](@ref) rather
than restoring raw flags, a mask edited against one sample rate can be applied
to a decimated or re-read version of the same record.

[`combine_masks`](@ref) takes the union of several masks over the same axis —
useful when two people review the same record, or when an automated detector
and a human review are merged:

```julia
final = combine_masks(human_mask, despike_mask, gps_dropout_mask)
```

## Writing the result out

[`write_cleaned`](@ref) exports the cleaned series as a delimited text file
with a `timestamp` column followed by one column per component:

```julia
write_cleaned("data/LEMI090_clean.csv", ta, mask)
```

To stay in the instrument's own format instead, clean first and hand the result
to the native writer — masked rows are written as `NaN`:

```julia
write_lemi424("data/LEMI090_clean.txt", cleaned_timearray(ta, mask))
```

For Metronix, blanking is usually the wrong move: the acquisition format has no
`NaN`, and downstream tools expect continuous runs. Use
[`write_metronix_site`](@ref) instead, which splits the record at the masked
intervals into separate `meas_*` directories — see
[Metronix Sites](metronix.md).

## From the app

Every function above accepts a [`TKApp`](@ref) in place of the
`(ta, mask)` pair, so a masking session done by hand in the window can be
picked up in code without unpacking anything:

```julia
app = run_tkapp("data/LEMI090.txt")   # mask a few intervals, then close

segments = good_segments(app; min_samples = 256)
write_mask("data/LEMI090_mask.csv", app)
write_cleaned("data/LEMI090_clean.csv", app)
```

See [TKApp Explorer](tkapp.md) for the interactive side.
