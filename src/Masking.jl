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

mask_interval!(mask::TimekeeperMask, start_time, end_time) = _set_interval!(mask, start_time, end_time, true)

unmask_interval!(mask::TimekeeperMask, start_time, end_time) = _set_interval!(mask, start_time, end_time, false)

masked_samples(mask::TimekeeperMask) = count(mask.masked)

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
        println(io, join(["timestamp"; String.(names)], delimiter))
        for i in eachindex(times)
            row = Vector{String}(undef, length(names) + 1)
            row[1] = string(times[i])
            for j in eachindex(names)
                row[j + 1] = string(vals[i, j])
            end
            println(io, join(row, delimiter))
        end
    end
    return path
end

function write_cleaned(path::AbstractString, ta::TimeArray, mask::TimekeeperMask; mode = :nan, delimiter = ',')
    return _write_timearray_csv(path, cleaned_timearray(ta, mask; mode = mode); delimiter = delimiter)
end

function write_mask(path::AbstractString, mask::TimekeeperMask; delimiter = ',')
    open(path, "w") do io
        println(io, join(["start", "stop"], delimiter))
        for (start_time, stop_time) in mask.intervals
            println(io, join((string(start_time), string(stop_time)), delimiter))
        end
    end
    return path
end

function _parse_datetime_token(token::AbstractString)
    clean = replace(strip(token), " UTC" => "", "Z" => "")
    return DateTime(clean)
end

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
