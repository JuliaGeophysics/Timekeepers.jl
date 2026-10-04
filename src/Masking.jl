# Masking.jl - masks for bad data intervals.
# Author: @pankajkmishra
#
# This file defines TimekeeperMask. A mask is a vector of good or bad flags,
# one for each sample. The mask also keeps a list of the masked time
# intervals, which always agrees with the flags. The file has:
# - the operations that mask and unmask intervals
# - the functions that make a clean series or the contiguous good segments
# - the functions that write a mask and its data to CSV files and read them

"""
    TimekeeperMask(ta::TimeArray)
    TimekeeperMask(timestamps, masked, intervals)

A good or bad flag for each sample of a time series, with the list of the
contiguous bad intervals.

A mask that you make from a `TimeArray` has all its samples good. Edit it with
[`mask_interval!`](@ref), [`unmask_interval!`](@ref) and
[`clear_mask!`](@ref). The `intervals` field always agrees with the flags.

# Fields
- `timestamps::Vector{T}` -- the time axis of the mask.
- `masked::BitVector` -- `true` marks a bad sample.
- `intervals::Vector{Tuple{T, T}}` -- the contiguous masked spans.

To make outputs, use [`cleaned_timearray`](@ref), [`good_segments`](@ref) and
[`sample_weights`](@ref). To save the mask and read it again, use
[`write_mask`](@ref) and [`read_mask`](@ref).
"""
mutable struct TimekeeperMask{T}
    timestamps::Vector{T}
    masked::BitVector
    intervals::Vector{Tuple{T, T}}
end

function TimekeeperMask(ta::TimeArray)
    times = collect(_ta_timestamps(ta))
    return TimekeeperMask(times, falses(length(times)), Tuple{eltype(times), eltype(times)}[])
end

function _ordered_pair(a, b)
    return a <= b ? (a, b) : (b, a)
end

function _refresh_intervals!(mask::TimekeeperMask)
    empty!(mask.intervals)
    active = false
    start_index = 1
    for i in eachindex(mask.masked)
        if mask.masked[i] && !active
            active = true
            start_index = i
        elseif !mask.masked[i] && active
            push!(mask.intervals, (mask.timestamps[start_index], mask.timestamps[i - 1]))
            active = false
        end
    end
    active && push!(mask.intervals, (mask.timestamps[start_index], mask.timestamps[end]))
    return mask
end

"""
    clear_mask!(mask::TimekeeperMask) -> TimekeeperMask

Mark all the samples as good and remove all the intervals. The function
changes `mask` and returns it.
"""
function clear_mask!(mask::TimekeeperMask)
    fill!(mask.masked, false)
    empty!(mask.intervals)
    return mask
end

function _set_interval!(mask::TimekeeperMask, start_time, end_time, value::Bool)
    lo, hi = _ordered_pair(start_time, end_time)
    for i in eachindex(mask.timestamps)
        if lo <= mask.timestamps[i] <= hi
            mask.masked[i] = value
        end
    end
    return _refresh_intervals!(mask)
end

"""
    mask_interval!(mask, start_time, end_time) -> TimekeeperMask

Mark each sample with a timestamp in `[start_time, end_time]` as bad. You can
give the two limits in either order. The function changes `mask` and returns
it.

```julia
mask = TimekeeperMask(ta)
mask_interval!(mask, DateTime(2020, 10, 4, 0, 10), DateTime(2020, 10, 4, 0, 20))
```
"""
mask_interval!(mask::TimekeeperMask, start_time, end_time) = _set_interval!(mask, start_time, end_time, true)

"""
    unmask_interval!(mask, start_time, end_time) -> TimekeeperMask

The opposite of [`mask_interval!`](@ref). Mark each sample in the closed
interval as good again. The function changes `mask` and returns it.
"""
unmask_interval!(mask::TimekeeperMask, start_time, end_time) = _set_interval!(mask, start_time, end_time, false)

"""
    masked_samples(mask::TimekeeperMask) -> Int

The number of samples that are bad now.
"""
masked_samples(mask::TimekeeperMask) = count(mask.masked)

"""
    sample_weights(mask::TimekeeperMask; good = 1.0, bad = 0.0) -> Vector
    sample_weights(app::TKApp; good = 1.0, bad = 0.0) -> Vector

A weight for each sample, for robust processing. The weight is `good` where
the mask is clear and `bad` where the mask is set. Use the weights to give a
mask to a weighted regression without removing samples.
"""
function sample_weights(mask::TimekeeperMask; good = 1.0, bad = 0.0)
    return [m ? bad : good for m in mask.masked]
end

function _assert_mask_matches(ta::TimeArray, mask::TimekeeperMask)
    length(_ta_timestamps(ta)) == length(mask.masked) ||
        error("Mask length $(length(mask.masked)) does not match TimeArray length $(length(_ta_timestamps(ta)))")
    return nothing
end

function _masked_meta(ta::TimeArray, mask::TimekeeperMask)
    metadata = _ta_meta(ta)
    base = metadata isa AbstractDict ? Dict{Symbol, Any}(metadata) : Dict{Symbol, Any}()
    base[:mask_intervals] = copy(mask.intervals)
    base[:masked_samples] = masked_samples(mask)
    return base
end

"""
    cleaned_timearray(ta::TimeArray, mask::TimekeeperMask; mode = :nan) -> TimeArray
    cleaned_timearray(app::TKApp; mode = :nan) -> TimeArray

Apply `mask` to `ta`.

- `mode = :nan` (the default) keeps the full time axis and writes `NaN` in
  the masked rows. Thus, the gaps stay visible and the interval between
  samples stays the same.
- `mode = :drop` removes the masked rows. The time axis is then not uniform.

The metadata of the result holds `:mask_intervals`, `:masked_samples` and
`:cleaning_mode`. If the length of the mask is not the length of `ta`, the
function gives an error.
"""
function cleaned_timearray(ta::TimeArray, mask::TimekeeperMask; mode = :nan)
    _assert_mask_matches(ta, mask)
    times = collect(_ta_timestamps(ta))
    vals = _ta_values(ta)
    names = _symbolize.(_ta_colnames(ta))
    meta = _masked_meta(ta, mask)
    if mode == :nan
        cleaned = Matrix{Float64}(vals)
        for i in eachindex(mask.masked)
            mask.masked[i] && (cleaned[i, :] .= NaN)
        end
        meta[:cleaning_mode] = :nan
        return TimeArray(times, cleaned, names, meta)
    elseif mode == :drop
        keep = .!mask.masked
        meta[:cleaning_mode] = :drop
        return TimeArray(times[keep], vals[keep, :], names, meta)
    else
        error("Unsupported cleaning mode: $mode")
    end
end

"""
    good_segments(ta::TimeArray, mask::TimekeeperMask; min_samples = 1) -> Vector{TimeArray}
    good_segments(app::TKApp; min_samples = 1) -> Vector{TimeArray}

Cut `ta` into the contiguous parts without a mask. The function removes each
part that is shorter than `min_samples`. Use this function to give clean data
to a processing step that needs windows without interruptions. For example,
use `min_samples = 256` for an FFT of 256 points.
"""
function good_segments(ta::TimeArray, mask::TimekeeperMask; min_samples = 1)
    _assert_mask_matches(ta, mask)
    times = collect(_ta_timestamps(ta))
    vals = _ta_values(ta)
    names = _symbolize.(_ta_colnames(ta))
    meta = _masked_meta(ta, mask)
    out = TimeArray[]
    active = false
    start_index = 1
    for i in eachindex(mask.masked)
        if !mask.masked[i] && !active
            active = true
            start_index = i
        elseif mask.masked[i] && active
            stop_index = i - 1
            stop_index - start_index + 1 >= min_samples &&
                push!(out, TimeArray(times[start_index:stop_index], vals[start_index:stop_index, :], names, meta))
            active = false
        end
    end
    if active
        stop_index = length(mask.masked)
        stop_index - start_index + 1 >= min_samples &&
            push!(out, TimeArray(times[start_index:stop_index], vals[start_index:stop_index, :], names, meta))
    end
    return out
end

"""
    combine_masks(first_mask, masks...) -> TimekeeperMask

The union of masks on the same time axis. A sample is bad if it is bad in one
or more inputs. All the masks must have the same length. The function returns
a new mask and does not change the inputs.
"""
function combine_masks(first_mask::TimekeeperMask, masks::TimekeeperMask...)
    combined = TimekeeperMask(copy(first_mask.timestamps), copy(first_mask.masked), copy(first_mask.intervals))
    for mask in masks
        length(mask.masked) == length(combined.masked) || error("Cannot combine masks with different lengths")
        combined.masked .|= mask.masked
    end
    return _refresh_intervals!(combined)
end

function _write_timearray_csv(path::AbstractString, ta::TimeArray; delimiter = ',')
    times = _ta_timestamps(ta)
    vals = _ta_values(ta)
    names = _symbolize.(_ta_colnames(ta))
    open(path, "w") do io
        print(io, "timestamp")
        for name in names
            print(io, delimiter, String(name))
        end
        println(io)
        for i in eachindex(times)
            print(io, times[i])
            for j in eachindex(names)
                print(io, delimiter, vals[i, j])
            end
            println(io)
        end
    end
    return path
end

"""
    write_cleaned(path, ta::TimeArray, mask::TimekeeperMask; mode = :nan, delimiter = ',') -> String
    write_cleaned(path, app::TKApp; mode = :nan, delimiter = ',') -> String

Write [`cleaned_timearray`](@ref) to a delimited text file. The file has a
`timestamp` column, then one column for each component. The function returns
`path`.

To keep the format of the instrument, use [`write_timekeeper`](@ref) or the
writer of the format.
"""
function write_cleaned(path::AbstractString, ta::TimeArray, mask::TimekeeperMask; mode = :nan, delimiter = ',')
    return _write_timearray_csv(path, cleaned_timearray(ta, mask; mode = mode); delimiter = delimiter)
end

"""
    write_mask(path, mask::TimekeeperMask; delimiter = ',') -> String
    write_mask(path, app::TKApp; delimiter = ',') -> String

Write the masked intervals to a file with the two columns `start,stop` and a
header row. The function returns `path`. To read the file, use
[`read_mask`](@ref).
"""
function write_mask(path::AbstractString, mask::TimekeeperMask; delimiter = ',')
    open(path, "w") do io
        print(io, "start", delimiter, "stop")
        println(io)
        for (start_time, stop_time) in mask.intervals
            print(io, start_time, delimiter, stop_time)
            println(io)
        end
    end
    return path
end

function _parse_datetime_token(token::AbstractString)
    clean = replace(strip(token), " UTC" => "", "Z" => "")
    return DateTime(clean)
end

"""
    read_mask(path, ta::TimeArray; delimiter = ',') -> TimekeeperMask

Make a mask on the time axis of `ta` from an interval file that
[`write_mask`](@ref) wrote. The function applies the intervals with
[`mask_interval!`](@ref). Thus, you can apply a mask that you saved for one
series to a different series that covers the same times.
"""
function read_mask(path::AbstractString, ta::TimeArray; delimiter = ',')
    mask = TimekeeperMask(ta)
    open(path, "r") do io
        for (line_number, line) in enumerate(eachline(io))
            line_number == 1 && continue
            isempty(strip(line)) && continue
            parts = split(line, delimiter)
            length(parts) >= 2 || error("Mask line $line_number must contain start and stop")
            mask_interval!(mask, _parse_datetime_token(parts[1]), _parse_datetime_token(parts[2]))
        end
    end
    return mask
end
