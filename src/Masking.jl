# Masking.jl - marking and removing bad data intervals.
# Author: @pankajkmishra
#
# Defines TimekeeperMask, a per-sample bad/good flag vector kept in sync with a
# derived list of masked time intervals. Provides the mask/unmask operations,
# extraction of cleaned series or contiguous good segments, and CSV
# persistence of both a mask and the data it describes.

"""
    TimekeeperMask(ta::TimeArray)
    TimekeeperMask(timestamps, masked, intervals)

Per-sample good/bad flags for a time series, together with the derived list of
contiguous bad intervals.

Constructing from a `TimeArray` gives an all-good mask over its timestamps.
Edit it with [`mask_interval!`](@ref), [`unmask_interval!`](@ref) and
[`clear_mask!`](@ref); the `intervals` field is kept in sync automatically.

# Fields
- `timestamps::Vector{T}` -- the time axis the mask refers to.
- `masked::BitVector` -- `true` marks a bad sample.
- `intervals::Vector{Tuple{T, T}}` -- derived contiguous masked spans.

Derive outputs with [`cleaned_timearray`](@ref), [`good_segments`](@ref) and
[`sample_weights`](@ref); persist with [`write_mask`](@ref) /
[`read_mask`](@ref).
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

Mark every sample good and drop all intervals. Modifies and returns `mask`.
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

Mark every sample with a timestamp in `[start_time, end_time]` as bad. The two
bounds may be given in either order. Modifies and returns `mask`.

```julia
mask = TimekeeperMask(ta)
mask_interval!(mask, DateTime(2020, 10, 4, 0, 10), DateTime(2020, 10, 4, 0, 20))
```
"""
mask_interval!(mask::TimekeeperMask, start_time, end_time) = _set_interval!(mask, start_time, end_time, true)

"""
    unmask_interval!(mask, start_time, end_time) -> TimekeeperMask

Inverse of [`mask_interval!`](@ref): mark every sample in the closed interval
good again. Modifies and returns `mask`.
"""
unmask_interval!(mask::TimekeeperMask, start_time, end_time) = _set_interval!(mask, start_time, end_time, false)

"""
    masked_samples(mask::TimekeeperMask) -> Int

Number of samples currently marked bad.
"""
masked_samples(mask::TimekeeperMask) = count(mask.masked)

"""
    sample_weights(mask::TimekeeperMask; good = 1.0, bad = 0.0) -> Vector
    sample_weights(app::TKApp; good = 1.0, bad = 0.0) -> Vector

Per-sample weights for robust processing: `good` where the mask is clear and
`bad` where it is set. Useful for feeding a mask into a weighted regression
without dropping samples.
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

- `mode = :nan` (default) keeps the time axis intact and writes `NaN` into
  masked rows, so gaps stay visible and sample spacing stays uniform.
- `mode = :drop` removes masked rows entirely, leaving a non-uniform axis.

The returned metadata carries `:mask_intervals`, `:masked_samples` and
`:cleaning_mode`. Errors if the mask length does not match `ta`.
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

Split `ta` into the contiguous unmasked runs, discarding any shorter than
`min_samples`. This is the usual way to hand clean data to a processing step
that needs uninterrupted windows, e.g. `min_samples = 256` for a 256-point FFT.
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

Union of several masks over the same time axis: a sample is bad if it is bad in
any input. All masks must have the same length. Returns a new mask; the inputs
are untouched.
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

Write [`cleaned_timearray`](@ref) to a delimited text file with a `timestamp`
column followed by one column per component. Returns `path`.

To stay in the instrument's own format instead, use [`write_timekeeper`](@ref)
or the format-specific writers.
"""
function write_cleaned(path::AbstractString, ta::TimeArray, mask::TimekeeperMask; mode = :nan, delimiter = ',')
    return _write_timearray_csv(path, cleaned_timearray(ta, mask; mode = mode); delimiter = delimiter)
end

"""
    write_mask(path, mask::TimekeeperMask; delimiter = ',') -> String
    write_mask(path, app::TKApp; delimiter = ',') -> String

Write the masked intervals to a two-column `start,stop` file with a header row.
Returns `path`. Read it back with [`read_mask`](@ref).
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

Rebuild a mask over `ta`'s time axis from an interval file written by
[`write_mask`](@ref). Intervals are applied with [`mask_interval!`](@ref), so a
mask saved against one series can be replayed onto another that covers the same
times.
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

function _set_mask_from_seconds!(mask::TimekeeperMask, offsets::AbstractVector{<:Real}, intervals)
    length(offsets) == length(mask.masked) || error("Offset length does not match mask length")
    clear_mask!(mask)
    for interval in intervals
        length(interval) >= 2 || continue
        lo, hi = _ordered_pair(Float64(interval[1]), Float64(interval[2]))
        for i in eachindex(offsets)
            if lo <= offsets[i] <= hi
                mask.masked[i] = true
            end
        end
    end
    return _refresh_intervals!(mask)
end
