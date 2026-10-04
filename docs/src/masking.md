# Masking & Cleaning

Field recordings contain intervals that you do not want. For example, a
technician walks near the sensor, a vehicle goes past, the GPS signal stops,
or the battery is replaced. Timekeepers does not delete these samples. It
records *where* they are in a [`TimekeeperMask`](@ref). The mask stays with
the data. From the data and the mask, Timekeepers makes the data structure
that the next processing step needs.

This separation is important. A mask is a few kilobytes of intervals. Thus,
you can:

- keep the mask under version control;
- review the mask;
- apply the mask again when you read the same file again;
- merge the mask with the mask of a colleague.

You cannot do these tasks after you delete the samples.

## Make and edit a mask

You make a mask on the time axis of a `TimeArray`. At the start, all the
samples are good:

```julia
using Timekeepers, Dates

ta   = load_lemi424("data/LEMI090.txt")
mask = TimekeeperMask(ta)

mask_interval!(mask, DateTime(2020, 10, 4, 0, 10), DateTime(2020, 10, 4, 0, 20))
mask_interval!(mask, DateTime(2020, 10, 4, 3, 45), DateTime(2020, 10, 4, 3, 47))

masked_samples(mask)   # 722
mask.intervals         # the two spans, in order
```

- [`mask_interval!`](@ref) marks a closed interval as bad.
- [`unmask_interval!`](@ref) marks the interval as good again.
- [`clear_mask!`](@ref) marks all the samples as good.

You can give the two limits in either order. After each edit, the mask
calculates its `intervals` list again. Thus, the list always agrees with the
current flags. Edits that touch or overlap become one interval automatically.

## Make clean data

Three functions use a mask and its data to make data for a processing step.

### Series with `NaN`

[`cleaned_timearray`](@ref) keeps the full time axis and writes `NaN` in the
masked rows. The interval between samples stays the same, and you can see the
gaps in a plot:

```julia
cleaned = cleaned_timearray(ta, mask)          # mode = :nan, the default
```

To remove the masked rows, give `mode = :drop`. The series is then shorter,
but its time axis is not uniform. Most spectral methods do not accept a time
axis that is not uniform. Thus, use `:nan`, unless the next step accepts
irregular samples.

### Contiguous good segments

[`good_segments`](@ref) cuts the record at the masked intervals. It returns
the good parts. Use this function when the next step needs windows without
interruptions:

```julia
segments = good_segments(ta, mask; min_samples = 256)
length(segments)                    # 3
```

`min_samples` removes the parts that are too short. Set it to the FFT length
that you will use.

### Weight for each sample

[`sample_weights`](@ref) makes a vector of weights. Use it for methods that
give bad samples a lower weight and do not remove them:

```julia
w = sample_weights(mask)                     # 1.0 good, 0.0 bad
w = sample_weights(mask; good = 1, bad = 0)  # integer weights
```

## Save a mask and apply it again

[`write_mask`](@ref) writes the intervals to a CSV file with the two columns
`start,stop`. [`read_mask`](@ref) makes a mask from that file again. The time
axis must cover the same period:

```julia
write_mask("data/LEMI090_mask.csv", mask)

# later, or on another machine
ta2   = load_lemi424("data/LEMI090.txt")
mask2 = read_mask("data/LEMI090_mask.csv", ta2)
mask2.masked == mask.masked   # true
```

[`read_mask`](@ref) applies each interval with [`mask_interval!`](@ref). It
does not copy the raw flags. Thus, you can apply a mask that you made at one
sample rate to a decimated version of the same record, or to a new read of
it.

[`combine_masks`](@ref) makes the union of masks on the same axis. Use it when
two persons review the same record, or to merge the result of an automatic
detector with a human review:

```julia
final = combine_masks(human_mask, despike_mask, gps_dropout_mask)
```

## Write the result

[`write_cleaned`](@ref) writes the clean series as a delimited text file. The
file has a `timestamp` column, then one column for each component:

```julia
write_cleaned("data/LEMI090_clean.csv", ta, mask)
```

To keep the format of the instrument, first clean the data. Then give the
result to the native writer. The writer writes the masked rows as `NaN`:

```julia
write_lemi424("data/LEMI090_clean.txt", cleaned_timearray(ta, mask))
```

For Metronix, do not put blanks in the data. The acquisition format has no
`NaN`, and the next tools need continuous runs. Use
[`write_metronix_site`](@ref). It cuts the record at the masked intervals and
writes each part in a separate `meas_*` directory. Refer to
[Metronix Sites](metronix.md).

## From the app

Each function above also accepts a [`TKApp`](@ref) in place of the pair
`(ta, mask)`. Thus, you can make a mask by hand in the window and then
continue in code. You do not have to unpack the app:

```julia
app = run_tkapp("data/LEMI090.txt")   # mask a few intervals, then close

segments = good_segments(app; min_samples = 256)
write_mask("data/LEMI090_mask.csv", app)
write_cleaned("data/LEMI090_clean.csv", app)
```

Refer to [TKApp Explorer](tkapp.md) for the interactive part.
