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
const LEMI424_OPTIONAL_DEFAULTS = Dict{Symbol, Any}(
    :temperature_e => NaN,
    :temperature_h => NaN,
    :e3 => NaN,
    :e4 => NaN,
    :battery => NaN,
    :elevation => NaN,
    :latitude => NaN,
    :lat_hemisphere => "N",
    :longitude => NaN,
    :lon_hemisphere => "E",
    :n_satellites => 0,
    :gps_fix => 0,
    :time_diff => NaN,
)

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
    length(parts) >= maximum(values(LEMI424_DEFAULT_COLUMN_INDEX)) ||
        error("LEMI-424 line $line_number has $(length(parts)) columns; expected at least $(maximum(values(LEMI424_DEFAULT_COLUMN_INDEX)))")

    year = parse(Int, parts[1])
    month = parse(Int, parts[2])
    day = parse(Int, parts[3])
    hour = parse(Int, parts[4])
    minute = parse(Int, parts[5])
    second = parse(Int, parts[6])
    stamp = DateTime(year, month, day, hour, minute, second)

    row = Dict{Symbol, Any}(:date => stamp)
    for (col, default) in LEMI424_OPTIONAL_DEFAULTS
        row[col] = default
    end
    max_known = min(length(parts), length(LEMI424_COLUMNS))
    for idx in 7:max_known
        col = LEMI424_COLUMNS[idx]
        raw = parts[idx]
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
    if isfinite(row[:latitude])
        row[:latitude] = _lemi_position(row[:latitude_raw]) * _lemi_hemisphere_sign(row[:lat_hemisphere])
    end
    if isfinite(row[:longitude])
        row[:longitude] = _lemi_position(row[:longitude_raw]) * _lemi_hemisphere_sign(row[:lon_hemisphere])
    end
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
    data = Float64[get(row, component, NaN) for row in rows]
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
        if haskey(rows[1], comp)
            channels[comp] = _lemi424_channel(rows, comp, abs_path, String(site))
        end
    end

    latitudes = Float64[get(row, :latitude, NaN) for row in rows]
    longitudes = Float64[get(row, :longitude, NaN) for row in rows]
    elevations = Float64[get(row, :elevation, NaN) for row in rows]
    batteries = Float64[get(row, :battery, NaN) for row in rows]
    metadata = Dict{Symbol, Any}(
        :source_file => abs_path,
        :sample_rate => 1.0,
        :start_time => rows[1][:date],
        :end_time => rows[end][:date],
        :n_samples => length(rows),
        :latitude => any(isfinite, latitudes) ? median(filter(isfinite, latitudes)) : NaN,
        :longitude => any(isfinite, longitudes) ? median(filter(isfinite, longitudes)) : NaN,
        :elevation => any(isfinite, elevations) ? median(filter(isfinite, elevations)) : NaN,
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
    required = maximum(keys(slots))
    ntokens >= required ||
        error("LEMI-424 line $line_number has $ntokens columns; expected at least $required for requested components")
    return DateTime(year, month, day, hour, minute, second)
end

function _lemi_aux_from_rows(rows)
    f64(sym, default = NaN) = Float64[Float64(get(r, sym, default)) for r in rows]
    str(sym, default) = String[String(get(r, sym, default)) for r in rows]
    return Dict{Symbol, AbstractVector}(
        :temperature_e => f64(:temperature_e),
        :temperature_h => f64(:temperature_h),
        :e3 => f64(:e3),
        :e4 => f64(:e4),
        :battery => f64(:battery),
        :elevation => f64(:elevation),
        :latitude_raw => f64(:latitude_raw),
        :lat_hemisphere => str(:lat_hemisphere, "N"),
        :longitude_raw => f64(:longitude_raw),
        :lon_hemisphere => str(:lon_hemisphere, "E"),
        :n_satellites => f64(:n_satellites, 0.0),
        :gps_fix => f64(:gps_fix, 0.0),
        :time_diff => f64(:time_diff, 0.0),
    )
end

function load_lemi424(path::AbstractString; components = LEMI424_DEFAULT_COMPONENTS,
                     site = _site_from_path(path), include_aux = true)
    abs_path = _as_path(path)
    requested = _symbolize.(collect(components))
    missing = [c for c in requested if !(c in LEMI424_DEFAULT_COMPONENTS)]
    isempty(missing) || error("LEMI-424 loader only supports: $(join(LEMI424_DEFAULT_COMPONENTS, ", "))")

    rows = _read_lemi424_rows(abs_path)
    n = length(rows)
    times = Vector{DateTime}(undef, n)
    values = Matrix{Float64}(undef, n, length(requested))
    @inbounds for (i, row) in enumerate(rows)
        times[i] = row[:date]
        for (j, comp) in enumerate(requested)
            values[i, j] = Float64(get(row, comp, NaN))
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
    include_aux && (meta[:aux_columns] = _lemi_aux_from_rows(rows))
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

_aux_at(::Nothing, ::Symbol, ::Integer, default) = default
function _aux_at(aux::AbstractDict, sym::Symbol, i::Integer, default)
    haskey(aux, sym) || return default
    vec = aux[sym]
    return i <= length(vec) ? vec[i] : default
end

function _lemi_row(run::TimekeeperRun, i::Integer, t::DateTime; aux = nothing)
    if aux !== nothing && haskey(aux, :latitude_raw) && haskey(aux, :lat_hemisphere)
        lat_raw = Float64(_aux_at(aux, :latitude_raw, i, NaN))
        lat_hemi = String(_aux_at(aux, :lat_hemisphere, i, "N"))
    else
        lat = get(run.metadata, :latitude, _value_for_component(run, :latitude, i; default = 0.0))
        lat_raw = _lemi_raw_position(lat)
        lat_hemi = lat < 0 ? "S" : "N"
    end
    if aux !== nothing && haskey(aux, :longitude_raw) && haskey(aux, :lon_hemisphere)
        lon_raw = Float64(_aux_at(aux, :longitude_raw, i, NaN))
        lon_hemi = String(_aux_at(aux, :lon_hemisphere, i, "E"))
    else
        lon = get(run.metadata, :longitude, _value_for_component(run, :longitude, i; default = 0.0))
        lon_raw = _lemi_raw_position(lon)
        lon_hemi = lon < 0 ? "W" : "E"
    end
    isfinite(lat_raw) || (lat_raw = 0.0)
    isfinite(lon_raw) || (lon_raw = 0.0)

    temp_e = Float64(_aux_at(aux, :temperature_e, i,
        _value_for_component(run, :temperature_e, i; default = NaN)))
    temp_h = Float64(_aux_at(aux, :temperature_h, i,
        _value_for_component(run, :temperature_h, i; default = NaN)))
    e3_val = Float64(_aux_at(aux, :e3, i, _value_for_component(run, :e3, i)))
    e4_val = Float64(_aux_at(aux, :e4, i, _value_for_component(run, :e4, i)))
    battery = Float64(_aux_at(aux, :battery, i,
        _value_for_component(run, :battery, i; default = get(run.metadata, :battery_end, NaN))))
    elev = Float64(_aux_at(aux, :elevation, i,
        _value_for_component(run, :elevation, i; default = get(run.metadata, :elevation, NaN))))
    n_sats_f = Float64(_aux_at(aux, :n_satellites, i,
        _value_for_component(run, :n_satellites, i; default = 0.0)))
    gps_f = Float64(_aux_at(aux, :gps_fix, i,
        _value_for_component(run, :gps_fix, i; default = 0.0)))
    tdiff = Float64(_aux_at(aux, :time_diff, i,
        _value_for_component(run, :time_diff, i; default = 0.0)))
    n_sats = isfinite(n_sats_f) ? round(Int, n_sats_f) : 0
    gps = isfinite(gps_f) ? round(Int, gps_f) : 0

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
        temp_e,
        temp_h,
        _value_for_component(run, :e1, i),
        _value_for_component(run, :e2, i),
        e3_val,
        e4_val,
        battery,
        elev,
        lat_raw,
        lat_hemi,
        lon_raw,
        lon_hemi,
        n_sats,
        gps,
        tdiff,
    )
end

function write_lemi424(path::AbstractString, run::TimekeeperRun;
                      timestamps = nothing, aux = nothing)
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
            println(io, _lemi_row(run, i, timestamps[i]; aux = aux))
        end
    end
    return path
end

function write_lemi424(path::AbstractString, ta::TimeArray; kwargs...)
    times = collect(_ta_timestamps(ta))
    (isempty(times) || !(first(times) isa DateTime)) &&
        error("LEMI-424 writer requires a non-empty TimeArray with DateTime timestamps")
    md = _ta_meta(ta)
    aux = md isa AbstractDict ? get(md, :aux_columns, nothing) : nothing
    run = from_timearray(ta; instrument = "LEMI-424", source_format = :timearray, kwargs...)
    return write_lemi424(path, run; timestamps = times, aux = aux)
end
