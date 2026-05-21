const GEOMAG_DEFAULT_COMPONENTS = [:bx, :by, :bz, :e1, :e2]
const GEOMAG_COLUMN_INDEX = Dict(
    :bx => 7,
    :by => 8,
    :bz => 9,
    :e1 => 10,
    :e2 => 11,
    :temperature_h => 12,
    :temperature_e => 13,
)

_geomag_blank(line::AbstractString) = all(isspace, line)
_geomag_header(line::AbstractString) = startswith(strip(line), ";")

function _parse_geomag_latlon(value::AbstractString)
    m = match(r"^\s*(\d+)\s+(\d+)'([0-9.]+)\"?([NSEW])", value)
    m === nothing && return NaN
    deg = parse(Float64, m.captures[1])
    minutes = parse(Float64, m.captures[2])
    seconds = parse(Float64, m.captures[3])
    hemi = m.captures[4]
    sign = hemi in ("S", "W") ? -1.0 : 1.0
    return sign * (deg + minutes / 60.0 + seconds / 3600.0)
end

function _parse_geomag_header!(metadata::Dict{Symbol, Any}, line::AbstractString)
    clean = strip(lstrip(strip(line), ';'))
    if startswith(clean, "MS:")
        metadata[:instrument_model] = strip(split(clean, '#'; limit = 2)[1][4:end])
        metadata[:instrument] = metadata[:instrument_model]
    elseif startswith(clean, "Sampling:")
        m = match(r"Sampling:\s*([0-9.]+)\s*sec", clean)
        if m !== nothing
            dt = parse(Float64, m.captures[1])
            dt > 0 && (metadata[:sample_rate] = 1.0 / dt)
        end
    elseif startswith(clean, "Latitude:")
        lat_match = match(r"Latitude:\s*([^;]+)", clean)
        lon_match = match(r"Longitude:\s*([^;]+)", clean)
        if lat_match !== nothing
            metadata[:latitude] = _parse_geomag_latlon(lat_match.captures[1])
        end
        if lon_match !== nothing
            metadata[:longitude] = _parse_geomag_latlon(lon_match.captures[1])
        end
        alt_match = match(r"Altitude:\s*([^;]+)", clean)
        if alt_match !== nothing
            alt = tryparse(Float64, strip(alt_match.captures[1]))
            metadata[:elevation] = alt === nothing ? NaN : alt
        end
    elseif startswith(clean, "Total Field:")
        for comp in ("X", "Y", "Z")
            m = match(Regex("$comp\\s*=\\s*([+-]?[0-9.]+)nT"), clean)
            if m !== nothing
                metadata[Symbol("total_field_" * lowercase(comp))] = parse(Float64, m.captures[1])
            end
        end
    end
    return metadata
end

function _geomag_time(parts)
    year = parse(Int, parts[1])
    month = parse(Int, parts[2])
    day = parse(Int, parts[3])
    hour = parse(Int, parts[4])
    minute = parse(Int, parts[5])
    second_float = parse(Float64, parts[6])
    whole_second = floor(Int, second_float)
    ms = round(Int, (second_float - whole_second) * 1000)
    return DateTime(year, month, day, hour, minute, whole_second) + Millisecond(ms)
end

function _parse_geomag_values!(values::AbstractMatrix, line::AbstractString, line_number::Integer, row::Integer, slots)
    parts = split(strip(line))
    required = maximum(keys(slots))
    length(parts) >= required ||
        error("GEOMAG line $line_number has $(length(parts)) columns; expected at least $required for requested components")
    for (col, out_col) in slots
        values[row, out_col] = parse(Float64, parts[col])
    end
    return _geomag_time(parts)
end

function _read_geomag_metadata(path::AbstractString)
    metadata = Dict{Symbol, Any}(
        :source_format => :geomag,
        :source_file => _as_path(path),
        :instrument => "GEOMAG",
        :instrument_model => "GEOMAG",
        :data_logger_manufacturer => "LEMI",
        :sample_rate => NaN,
        :latitude => NaN,
        :longitude => NaN,
        :elevation => NaN,
    )
    open(path, "r") do io
        for line in eachline(io)
            _geomag_header(line) || break
            _parse_geomag_header!(metadata, line)
        end
    end
    return metadata
end

function _count_geomag_records(path::AbstractString)
    n = 0
    open(path, "r") do io
        for line in eachline(io)
            (_geomag_blank(line) || _geomag_header(line)) && continue
            n += 1
        end
    end
    n > 0 || error("GEOMAG file is empty: $path")
    return n
end

function _sample_rate_from_times(times::Vector{DateTime}, fallback::Real)
    isfinite(fallback) && fallback > 0 && return Float64(fallback)
    length(times) < 2 && return 1.0
    dt_ms = Dates.value(times[2] - times[1])
    return dt_ms > 0 ? 1000.0 / dt_ms : 1.0
end

function load_geomag(path::AbstractString; components = GEOMAG_DEFAULT_COMPONENTS, site = _site_from_path(path))
    abs_path = _as_path(path)
    requested = _symbolize.(collect(components))
    missing = [c for c in requested if !haskey(GEOMAG_COLUMN_INDEX, c)]
    isempty(missing) || error("GEOMAG loader does not support components: $(join(missing, ", "))")

    n = _count_geomag_records(abs_path)
    times = Vector{DateTime}(undef, n)
    values = Matrix{Float64}(undef, n, length(requested))
    slots = Dict(GEOMAG_COLUMN_INDEX[component] => i for (i, component) in enumerate(requested))
    row = 0
    open(abs_path, "r") do io
        for (line_number, line) in enumerate(eachline(io))
            (_geomag_blank(line) || _geomag_header(line)) && continue
            row += 1
            times[row] = _parse_geomag_values!(values, line, line_number, row, slots)
        end
    end

    metadata = _read_geomag_metadata(abs_path)
    fs = _sample_rate_from_times(times, get(metadata, :sample_rate, NaN))
    metadata[:site] = String(site)
    metadata[:sample_rate] = fs
    metadata[:start_time] = first(times)
    metadata[:end_time] = last(times)
    metadata[:n_samples] = n
    metadata[:reader] = :load_geomag
    metadata[:units] = Dict(
        :bx => "nT",
        :by => "nT",
        :bz => "nT",
        :e1 => "mV",
        :e2 => "mV",
        :temperature_h => "C",
        :temperature_e => "C",
    )
    return TimeArray(times, values, requested, metadata)
end

function read_geomag(path::AbstractString; site = _site_from_path(path), include_aux = true)
    comps = include_aux ? [:bx, :by, :bz, :e1, :e2, :temperature_h, :temperature_e] : GEOMAG_DEFAULT_COMPONENTS
    ta = load_geomag(path; components = comps, site = site)
    metadata = Dict{Symbol, Any}(_ta_meta(ta))
    return from_timearray(
        ta;
        site = String(site),
        instrument = string(get(metadata, :instrument, "GEOMAG")),
        source_format = :geomag,
        units = get(metadata, :units, Dict{Symbol, String}()),
        metadata = metadata,
    )
end

function _geomag_seconds(t::DateTime)
    return second(t) + millisecond(t) / 1000
end

function write_geomag(path::AbstractString, ta::TimeArray)
    times = collect(_ta_timestamps(ta))
    (isempty(times) || !(first(times) isa DateTime)) &&
        error("GEOMAG writer requires a non-empty TimeArray with DateTime timestamps")
    vals = _ta_values(ta)
    names = _symbolize.(_ta_colnames(ta))
    idx(name) = findfirst(==(name), names)
    bx = idx(:bx)
    by = idx(:by)
    bz = idx(:bz)
    e1 = idx(:e1)
    e2 = idx(:e2)
    bx === nothing && error("GEOMAG writer requires :bx")
    by === nothing && error("GEOMAG writer requires :by")
    bz === nothing && error("GEOMAG writer requires :bz")
    e1 === nothing && error("GEOMAG writer requires :e1")
    e2 === nothing && error("GEOMAG writer requires :e2")

    metadata = _ta_meta(ta)
    fs = _sample_rate_from_timearray(ta)
    open(path, "w") do io
        println(io, "; MS:", get(metadata, :instrument_model, "GEOMAG"))
        println(io, "; Date: ", Dates.format(first(times), dateformat"yyyy/mm/dd"), "; Time: ", Dates.format(first(times), dateformat"HH:MM:SS"))
        println(io, "; Sampling: ", @sprintf("%.3g", 1 / fs), " sec")
        println(io, ";")
        println(io, ";    Date        Time    X [nT]   Y [nT]   Z [nT]    Ex[mV]   Ey[mV]   Ts[C] Te[C]")
        println(io, ";")
        @inbounds for i in eachindex(times)
            t = times[i]
            @printf(
                io,
                "%04d %02d %02d  %02d %02d %05.2f %+09.3f %+09.3f %+09.3f %+08.3f %+08.3f %+.1f %+.1f\n",
                year(t),
                month(t),
                day(t),
                hour(t),
                minute(t),
                _geomag_seconds(t),
                vals[i, bx],
                vals[i, by],
                vals[i, bz],
                vals[i, e1],
                vals[i, e2],
                NaN,
                NaN,
            )
        end
    end
    return path
end

function write_geomag(path::AbstractString, run::TimekeeperRun)
    return write_geomag(path, to_timearray(run; components = GEOMAG_DEFAULT_COMPONENTS))
end
