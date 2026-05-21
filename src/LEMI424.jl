const LEMI424_COLUMNS = [
    :year,
    :month,
    :day,
    :hour,
    :minute,
    :second,
    :bx,
    :by,
    :bz,
    :temperature_e,
    :temperature_h,
    :e1,
    :e2,
    :e3,
    :e4,
    :battery,
    :elevation,
    :latitude,
    :lat_hemisphere,
    :longitude,
    :lon_hemisphere,
    :n_satellites,
    :gps_fix,
    :time_diff,
]

const LEMI424_DATA_COLUMNS = LEMI424_COLUMNS[7:end]
const LEMI424_DEFAULT_COMPONENTS = [:bx, :by, :bz, :e1, :e2]
const LEMI424_DEFAULT_COLUMN_INDEX = Dict(:bx => 7, :by => 8, :bz => 9, :e1 => 12, :e2 => 13)

_lemi424_blank(line::AbstractString) = all(isspace, line)

function _lemi_position(position)
    scaled = Float64(position) / 100.0
    degrees = floor(Int, scaled)
    minutes = (scaled - degrees) * 100.0
    return degrees + minutes / 60.0
end

_lemi_hemisphere_sign(h) = uppercase(String(h)) in ("S", "W") ? -1.0 : 1.0

function _parse_lemi424_line(line::AbstractString, line_number::Integer)
    parts = split(strip(line))
    isempty(parts) && return nothing
    length(parts) == length(LEMI424_COLUMNS) ||
        error("LEMI-424 line $line_number has $(length(parts)) columns; expected $(length(LEMI424_COLUMNS))")

    year = parse(Int, parts[1])
    month = parse(Int, parts[2])
    day = parse(Int, parts[3])
    hour = parse(Int, parts[4])
    minute = parse(Int, parts[5])
    second = parse(Int, parts[6])
    stamp = DateTime(year, month, day, hour, minute, second)

    row = Dict{Symbol, Any}(:date => stamp)
    for (col, raw) in zip(LEMI424_COLUMNS[7:end], parts[7:end])
        if col in (:lat_hemisphere, :lon_hemisphere)
            row[col] = raw
        elseif col in (:n_satellites, :gps_fix)
            row[col] = parse(Int, raw)
        else
            row[col] = parse(Float64, raw)
        end
    end

    row[:latitude_raw] = row[:latitude]
    row[:longitude_raw] = row[:longitude]
    row[:latitude] = _lemi_position(row[:latitude_raw]) * _lemi_hemisphere_sign(row[:lat_hemisphere])
    row[:longitude] = _lemi_position(row[:longitude_raw]) * _lemi_hemisphere_sign(row[:lon_hemisphere])
    return row
end

function _read_lemi424_rows(path::AbstractString)
    rows = Dict{Symbol, Any}[]
    open(path, "r") do io
        for (line_number, line) in enumerate(eachline(io))
            row = _parse_lemi424_line(line, line_number)
            row === nothing || push!(rows, row)
        end
    end
    isempty(rows) && error("LEMI-424 file is empty: $path")
    return rows
end

function _lemi424_channel(rows, component::Symbol, path::AbstractString, site::String)
    data = Float64[row[component] for row in rows]
    start = rows[1][:date]
    header = Dict{String, Any}(
        "format" => "LEMI-424",
        "site" => site,
        "source_file" => path,
        "sample_rate" => 1.0,
    )
    return TimekeeperChannel(component, data, 1.0, start, component_units(component), path, header)
end

function read_lemi424(path::AbstractString; site = _site_from_path(path), include_aux = true)
    abs_path = _as_path(path)
    rows = _read_lemi424_rows(abs_path)
    comps = include_aux ? LEMI424_DATA_COLUMNS : LEMI424_DEFAULT_COMPONENTS
    numeric_comps = [c for c in comps if !(c in (:lat_hemisphere, :lon_hemisphere))]
    channels = Dict{Symbol, TimekeeperChannel}()
    for comp in numeric_comps
        if comp in (:latitude, :longitude, :n_satellites, :gps_fix, :time_diff, :battery, :elevation) || haskey(rows[1], comp)
            channels[comp] = _lemi424_channel(rows, comp, abs_path, String(site))
        end
    end

    latitudes = Float64[row[:latitude] for row in rows]
    longitudes = Float64[row[:longitude] for row in rows]
    elevations = Float64[row[:elevation] for row in rows]
    batteries = Float64[row[:battery] for row in rows]
    metadata = Dict{Symbol, Any}(
        :source_file => abs_path,
        :sample_rate => 1.0,
        :start_time => rows[1][:date],
        :end_time => rows[end][:date],
        :n_samples => length(rows),
        :latitude => median(latitudes),
        :longitude => median(longitudes),
        :elevation => median(elevations),
        :battery_start => first(batteries),
        :battery_end => last(batteries),
        :instrument_model => "LEMI-424",
        :data_logger_manufacturer => "LEMI",
    )
    return TimekeeperRun(String(site), "LEMI-424", :lemi424, channels, metadata)
end

function _count_lemi424_records(path::AbstractString)
    n = 0
    open(path, "r") do io
        for line in eachline(io)
            _lemi424_blank(line) || (n += 1)
        end
    end
    n > 0 || error("LEMI-424 file is empty: $path")
    return n
end

function _parse_lemi424_timearray_line!(values::AbstractMatrix, line::AbstractString, line_number::Integer, row::Integer, slots)
    ntokens = 0
    year = month = day = hour = minute = second = 0
    for raw in eachsplit(line)
        ntokens += 1
        if ntokens == 1
            year = parse(Int, raw)
        elseif ntokens == 2
            month = parse(Int, raw)
        elseif ntokens == 3
            day = parse(Int, raw)
        elseif ntokens == 4
            hour = parse(Int, raw)
        elseif ntokens == 5
            minute = parse(Int, raw)
        elseif ntokens == 6
            second = parse(Int, raw)
        elseif haskey(slots, ntokens)
            values[row, slots[ntokens]] = parse(Float64, raw)
        end
    end
    ntokens == length(LEMI424_COLUMNS) ||
        error("LEMI-424 line $line_number has $ntokens columns; expected $(length(LEMI424_COLUMNS))")
    return DateTime(year, month, day, hour, minute, second)
end

function load_lemi424(path::AbstractString; components = LEMI424_DEFAULT_COMPONENTS, site = _site_from_path(path))
    abs_path = _as_path(path)
    requested = _symbolize.(collect(components))
    missing = [c for c in requested if !(c in LEMI424_DEFAULT_COMPONENTS)]
    isempty(missing) || error("LEMI-424 loader only supports: $(join(LEMI424_DEFAULT_COMPONENTS, ", "))")

    n = _count_lemi424_records(abs_path)
    times = Vector{DateTime}(undef, n)
    values = Matrix{Float64}(undef, n, length(requested))
    slots = Dict(LEMI424_DEFAULT_COLUMN_INDEX[component] => i for (i, component) in enumerate(requested))
    row = 0
    open(abs_path, "r") do io
        for (line_number, line) in enumerate(eachline(io))
            _lemi424_blank(line) && continue
            row += 1
            times[row] = _parse_lemi424_timearray_line!(values, line, line_number, row, slots)
        end
    end

    meta = Dict{Symbol, Any}(
        :site => String(site),
        :instrument => "LEMI-424",
        :source_format => :lemi424,
        :source_file => abs_path,
        :sample_rate => 1.0,
        :start_time => first(times),
        :end_time => last(times),
        :n_samples => n,
        :instrument_model => "LEMI-424",
        :data_logger_manufacturer => "LEMI",
        :reader => :load_lemi424,
        :units => Dict(c => component_units(c) for c in requested),
    )
    return TimeArray(times, values, requested, meta)
end

function _lemi_raw_position(decimal_degrees)
    x = abs(Float64(decimal_degrees))
    degrees = floor(Int, x)
    minutes = (x - degrees) * 60.0
    return degrees * 100.0 + minutes
end

function _value_for_component(run::TimekeeperRun, component::Symbol, i::Integer; default = 0.0)
    haskey(run.channels, component) || return default
    data = run.channels[component].data
    return i <= length(data) ? data[i] : default
end

function _lemi_row(run::TimekeeperRun, i::Integer, t::DateTime)
    lat = get(run.metadata, :latitude, _value_for_component(run, :latitude, i; default = 0.0))
    lon = get(run.metadata, :longitude, _value_for_component(run, :longitude, i; default = 0.0))
    lat_hemi = lat < 0 ? "S" : "N"
    lon_hemi = lon < 0 ? "W" : "E"
    return @sprintf(
        "%04d %02d %02d %02d %02d %02d %.3f %.3f %.3f %.2f %.2f %.3f %.3f %.3f %.3f %.2f %.1f %.5f %s %.5f %s %d %d %.0f",
        year(t),
        month(t),
        day(t),
        hour(t),
        minute(t),
        second(t),
        _value_for_component(run, :bx, i),
        _value_for_component(run, :by, i),
        _value_for_component(run, :bz, i),
        _value_for_component(run, :temperature_e, i; default = NaN),
        _value_for_component(run, :temperature_h, i; default = NaN),
        _value_for_component(run, :e1, i),
        _value_for_component(run, :e2, i),
        _value_for_component(run, :e3, i),
        _value_for_component(run, :e4, i),
        _value_for_component(run, :battery, i; default = get(run.metadata, :battery_end, NaN)),
        _value_for_component(run, :elevation, i; default = get(run.metadata, :elevation, NaN)),
        _lemi_raw_position(lat),
        lat_hemi,
        _lemi_raw_position(lon),
        lon_hemi,
        round(Int, _value_for_component(run, :n_satellites, i; default = 0.0)),
        round(Int, _value_for_component(run, :gps_fix, i; default = 0.0)),
        _value_for_component(run, :time_diff, i; default = 0.0),
    )
end

function write_lemi424(path::AbstractString, run::TimekeeperRun; timestamps = nothing)
    fs = sampling_rate(run)
    isapprox(fs, 1.0; atol = 1e-9) ||
        error("LEMI-424 text writer expects 1 Hz data. Got sample_rate=$fs")
    n = minimum(length(ch.data) for ch in values(run.channels))
    if timestamps === nothing
        st = start_time(run)
        timestamps = [st + Second(round(Int, (i - 1) / fs)) for i in 1:n]
    end
    length(timestamps) >= n ||
        error("write_lemi424: got $(length(timestamps)) timestamps for $n samples")
    open(path, "w") do io
        for i in 1:n
            println(io, _lemi_row(run, i, timestamps[i]))
        end
    end
    return path
end

function write_lemi424(path::AbstractString, ta::TimeArray; kwargs...)
    times = collect(_ta_timestamps(ta))
    (isempty(times) || !(first(times) isa DateTime)) &&
        error("LEMI-424 writer requires a non-empty TimeArray with DateTime timestamps")
    run = from_timearray(ta; instrument = "LEMI-424", source_format = :timearray, kwargs...)
    return write_lemi424(path, run; timestamps = times)
end
