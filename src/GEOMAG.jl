# GEOMAG.jl - GEOMAG text format reader and writer.
# Author: @pankajkmishra
#
# A GEOMAG file starts with a block of header lines. Each header line starts
# with ';'. The header gives the instrument model, the sampling interval and
# the position of the station. Then the file has one record for each sample:
# six date and time fields (fractional seconds), then the magnetic, electric
# and temperature columns. The reader parses the records in one pass into a
# matrix that is already made. It returns a TimeArray (load_geomag) or a
# TimekeeperRun (read_geomag). write_geomag writes the same layout, with the
# header block

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

const GEOMAG_DATE_TOKENS = 6                                   # year, month, day, hour, minute, second
const GEOMAG_MAX_COLUMN = maximum(values(GEOMAG_COLUMN_INDEX))

_geomag_blank(line::AbstractString) = all(isspace, line)
_geomag_header(line::AbstractString) = startswith(strip(line), ";")

#---------- column slots -----

"""
    _geomag_slot_table(components) -> Vector{Int}

The map from the position of a token (from 1) to a column of the destination
matrix (0 = the reader does not need the token). The function accepts the
components to read, in the output order. It returns a plain `Vector{Int}`.
Thus, the loop for each line finds a value in an array and does not calculate
a `Dict` hash for each token.
"""
function _geomag_slot_table(components)
    table = zeros(Int, GEOMAG_MAX_COLUMN)
    for (i, c) in enumerate(components)
        table[GEOMAG_COLUMN_INDEX[c]] = i
    end
    return table
end

"""
    _geomag_required_columns(components) -> Int

The smallest number of tokens that a data line must have to give all the
requested components. The function accepts the components and returns the
highest column index of them.
"""
_geomag_required_columns(components) =
    isempty(components) ? 0 : maximum(GEOMAG_COLUMN_INDEX[c] for c in components)

#---------- header parsing -----

"""
    _parse_geomag_latlon(value) -> Float64

Parse a GEOMAG header coordinate in the form `60 35'14.4"N`. The function
accepts the text of the field. It returns decimal degrees with a sign. If the
text does not agree with the form, it returns `NaN`.
"""
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

"""
    _parse_geomag_header!(metadata, line) -> metadata

Add the data of one header line (it starts with `;`) to the metadata
dictionary. The function accepts the dictionary and the line. It finds the
instrument model, the sampling interval, the position and the total-field
entries. It ignores all other data. It returns the changed dictionary.
"""
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

#---------- line parsing -----

"""
    _geomag_timestamp(yr, mo, dy, hr, mi, sec_float) -> DateTime

Make the timestamp of a sample. The function accepts the five integer date
and time fields and the field of fractional seconds. It returns the instant,
with the fraction of a second as milliseconds.
"""
function _geomag_timestamp(yr::Int, mo::Int, dy::Int, hr::Int, mi::Int, sec_float::Float64)
    whole_second = floor(Int, sec_float)
    ms = round(Int, (sec_float - whole_second) * 1000)
    return DateTime(yr, mo, dy, hr, mi, whole_second) + Millisecond(ms)
end

"""
    _parse_geomag_record!(vals, line, line_number, row, value_slot, required;
                          aux=nothing, aux_slot=nothing) -> DateTime

Parse one data line directly into the row `row` of `vals`, without a
temporary vector of tokens. The function accepts:
- the destination matrix;
- the text of the line;
- its line number in the file, from 1 (for error messages);
- the destination row;
- the table from token to column, from [`_geomag_slot_table`](@ref);
- the minimum number of tokens;
- an optional auxiliary destination and table.

It changes the matrices and returns the timestamp of the sample. If the line
does not have an auxiliary column, that column keeps the `NaN` that the
caller put in it.
"""
function _parse_geomag_record!(vals::AbstractMatrix, line::AbstractString,
                               line_number::Integer, row::Integer,
                               value_slot::Vector{Int}, required::Int;
                               aux::Union{AbstractMatrix, Nothing} = nothing,
                               aux_slot::Union{Vector{Int}, Nothing} = nothing)
    yr = mo = dy = hr = mi = 0
    sec_float = 0.0
    ntokens = 0
    for raw in eachsplit(line)                             # lazy: no SubString vector
        ntokens += 1
        if ntokens <= GEOMAG_DATE_TOKENS
            if ntokens == GEOMAG_DATE_TOKENS
                sec_float = parse(Float64, raw)            # seconds with a fraction
            else
                v = parse(Int, raw)
                ntokens == 1 ? (yr = v) :
                ntokens == 2 ? (mo = v) :
                ntokens == 3 ? (dy = v) : (ntokens == 4 ? (hr = v) : (mi = v))
            end
        elseif ntokens > GEOMAG_MAX_COLUMN
            break                                          # ignore the extra fields at the end
        else
            @inbounds slot = value_slot[ntokens]
            if slot != 0
                @inbounds vals[row, slot] = parse(Float64, raw)
            elseif aux !== nothing && aux_slot !== nothing
                @inbounds a = aux_slot[ntokens]
                a == 0 || (@inbounds aux[row, a] = parse(Float64, raw))
            end
        end
    end
    ntokens >= required ||
        error("GEOMAG line $line_number has $ntokens columns; expected at least $required for requested components")
    return _geomag_timestamp(yr, mo, dy, hr, mi, sec_float)
end

#---------- file reading -----

"""
    _read_geomag_metadata(path) -> Dict{Symbol, Any}

Read the header block at the start of a GEOMAG file. The function accepts the
path and returns the metadata dictionary. It stops at the first line that is
not a header line. Thus, it never reads the data.
"""
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

"""
    _count_geomag_records(path) -> Int

Count the data lines in a GEOMAG file. The function ignores blank lines and
header lines. It accepts the path and returns the number of records. If the
file has no records, it gives an error.
"""
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

"""
    _sample_rate_from_times(times, fallback) -> Float64

Find the sample rate in Hz. The function accepts the timestamps and a fallback
value from the header. If the fallback is finite and positive, it returns the
fallback. If not, it returns the rate from the first gap between timestamps.
The default is `1.0`.
"""
function _sample_rate_from_times(times::Vector{DateTime}, fallback::Real)
    isfinite(fallback) && fallback > 0 && return Float64(fallback)
    length(times) < 2 && return 1.0
    dt_ms = Dates.value(times[2] - times[1])
    return dt_ms > 0 ? 1000.0 / dt_ms : 1.0
end

"""
    _fill_geomag_records!(times, vals, io, value_slot, required, aux, aux_slot) -> Int

Parse each data line of `io` into `times` and `vals`. The function accepts:
- the destination vector of timestamps and the destination matrix;
- an open stream;
- the table from token to column;
- the minimum number of tokens;
- the optional auxiliary destination and table.

It returns the number of records that it filled. It is separate from
[`load_geomag`](@ref). Thus, the row counter is a plain local variable. If the
`open` closure captured and changed the counter, Julia would box it.
"""
function _fill_geomag_records!(times::Vector{DateTime}, vals::AbstractMatrix, io::IO,
                               value_slot::Vector{Int}, required::Int,
                               aux::Union{AbstractMatrix, Nothing},
                               aux_slot::Union{Vector{Int}, Nothing})
    row = 0
    for (line_number, line) in enumerate(eachline(io))
        (_geomag_blank(line) || _geomag_header(line)) && continue
        row += 1
        @inbounds times[row] = _parse_geomag_record!(vals, line, line_number, row,
            value_slot, required; aux = aux, aux_slot = aux_slot)
    end
    return row
end

#---------- readers -----

"""
    load_geomag(path; components, site, include_aux=true) -> TimeArray

Load a GEOMAG text file as a `TimeArray`. The function accepts:
- the file path;
- the components to read (keys of `GEOMAG_COLUMN_INDEX`);
- an optional site name;
- a flag to add the other temperature columns to the metadata as
  `:aux_columns`.

It returns the array. The sample rate comes from the header. If the header
does not give it, the rate comes from the first gap between timestamps.
"""
function load_geomag(path::AbstractString; components = GEOMAG_DEFAULT_COMPONENTS,
                     site = _site_from_path(path), include_aux = true)
    abs_path = _as_path(path)
    requested = _symbolize.(collect(components))
    missing = [c for c in requested if !haskey(GEOMAG_COLUMN_INDEX, c)]
    isempty(missing) || error("GEOMAG loader does not support components: $(join(missing, ", "))")

    n = _count_geomag_records(abs_path)
    times = Vector{DateTime}(undef, n)
    vals = Matrix{Float64}(undef, n, length(requested))
    value_slot = _geomag_slot_table(requested)
    required = _geomag_required_columns(requested)

    aux_components = include_aux ? Symbol[c for c in (:temperature_h, :temperature_e) if !(c in requested)] : Symbol[]
    aux_values = isempty(aux_components) ? nothing : fill(NaN, n, length(aux_components))
    aux_slot = isempty(aux_components) ? nothing : _geomag_slot_table(aux_components)

    open(io -> _fill_geomag_records!(times, vals, io, value_slot, required, aux_values, aux_slot),
         abs_path, "r")

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
    if aux_values !== nothing
        metadata[:aux_columns] = Dict{Symbol, AbstractVector}(
            c => aux_values[:, i] for (i, c) in enumerate(aux_components)
        )
    end
    return TimeArray(times, vals, requested, metadata)
end

"""
    read_geomag(path; site, include_aux=true) -> TimekeeperRun

Read a GEOMAG text file into a run with one channel for each component. The
function accepts the file path, an optional site name and a flag for the
temperature channels. It returns a `TimekeeperRun`.
"""
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

#---------- writers -----

"""
    _geomag_seconds(t) -> Float64

The seconds in the minute of a timestamp. The function accepts the instant. It
returns the whole seconds and the fraction in milliseconds, as in the GEOMAG
time field.
"""
function _geomag_seconds(t::DateTime)
    return second(t) + millisecond(t) / 1000
end

"""
    write_geomag(path, ta) -> String

Write a `TimeArray` as a GEOMAG text file. The function accepts the
destination path and the array, and returns the path. If the array does not
have one of `:bx`, `:by`, `:bz`, `:e1` and `:e2`, the function writes 0 for
it. Thus, the fixed column layout of the format stays the same. The
temperatures come from the metadata `:aux_columns`. If the metadata does not
have them, they are `NaN`.
"""
function write_geomag(path::AbstractString, ta::TimeArray)
    times = collect(_ta_timestamps(ta))
    (isempty(times) || !(first(times) isa DateTime)) &&
        error("GEOMAG writer requires a non-empty TimeArray with DateTime timestamps")
    vals = _ta_values(ta)
    names = _symbolize.(_ta_colnames(ta))
    # The format has a fixed layout of five columns. If the array does not have
    # a channel, the function writes 0 for it. The logger does the same for an
    # input that is not connected, and the LEMI-424 writer also does this
    column(name) = (j = findfirst(==(name), names); j === nothing ? nothing : view(vals, :, j))
    bx, by, bz, e1, e2 = column.((:bx, :by, :bz, :e1, :e2))
    value(col, i) = col === nothing ? 0.0 : col[i]

    metadata = _ta_meta(ta)
    aux = metadata isa AbstractDict ? get(metadata, :aux_columns, nothing) : nothing
    th_vec = aux isa AbstractDict ? get(aux, :temperature_h, nothing) : nothing
    te_vec = aux isa AbstractDict ? get(aux, :temperature_e, nothing) : nothing
    fs = _sample_rate_from_timearray(ta)
    open(path, "w") do io
        println(io, "; MS:", (metadata isa AbstractDict ? get(metadata, :instrument_model, "GEOMAG") : "GEOMAG"))
        println(io, "; Date: ", Dates.format(first(times), dateformat"yyyy/mm/dd"), "; Time: ", Dates.format(first(times), dateformat"HH:MM:SS"))
        println(io, "; Sampling: ", @sprintf("%.3g", 1 / fs), " sec")
        println(io, ";")
        println(io, ";    Date        Time    X [nT]   Y [nT]   Z [nT]    Ex[mV]   Ey[mV]   Ts[C] Te[C]")
        println(io, ";")
        @inbounds for i in eachindex(times)
            t = times[i]
            th = th_vec !== nothing && i <= length(th_vec) ? Float64(th_vec[i]) : NaN
            te = te_vec !== nothing && i <= length(te_vec) ? Float64(te_vec[i]) : NaN
            @printf(
                io,
                "%04d %02d %02d  %02d %02d %05.2f %+09.3f %+09.3f %+09.3f %+08.3f %+08.3f %+.1f %+.1f\n",
                year(t),
                month(t),
                day(t),
                hour(t),
                minute(t),
                _geomag_seconds(t),
                value(bx, i),
                value(by, i),
                value(bz, i),
                value(e1, i),
                value(e2, i),
                th,
                te,
            )
        end
    end
    return path
end

"""
    write_geomag(path, run) -> String

Write a run as a GEOMAG text file. The function accepts the destination path
and the run, and returns the path. It first keeps only the
`GEOMAG_DEFAULT_COMPONENTS` of the run.
"""
function write_geomag(path::AbstractString, run::TimekeeperRun)
    return write_geomag(path, to_timearray(run; components = GEOMAG_DEFAULT_COMPONENTS))
end
