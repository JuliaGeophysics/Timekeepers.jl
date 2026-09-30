# LEMI424.jl - LEMI-424 text format reader and writer.
# Author: @pankajkmishra
#
# The LEMI-424 logger writes one whitespace-separated record per second: six
# date/time fields followed by magnetic, electric, temperature, GPS and
# housekeeping columns. Records are parsed in a single forward pass into
# column-oriented Float64 storage, from which both a TimekeeperRun
# (read_lemi424) and a TimeArray (load_lemi424) are built. Files with extra or
# missing trailing columns are tolerated. The writer reproduces the original
# 24-field layout, reusing auxiliary columns so a round trip is lossless.

"""
    LEMI424_COLUMNS

Names of the 24 whitespace-separated fields in a LEMI-424 record, in file
order: six date/time fields, then the magnetic, electric, temperature, GPS and
housekeeping columns. Used by the reader to index tokens and by the writer to
reproduce the original layout.
"""
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

const LEMI424_HEMISPHERE_COLUMNS = (:lat_hemisphere, :lon_hemisphere)

"""
Numeric data columns, in the fixed slot order used by every parsed record
matrix. Excludes the six date/time tokens and the two hemisphere letters.
"""
const LEMI424_NUMERIC_COLUMNS =
    Symbol[c for c in LEMI424_DATA_COLUMNS if !(c in LEMI424_HEMISPHERE_COLUMNS)]

"""
Map component -> its column index in a record matrix. Lookup is done once per
read, never inside the per-line loop.
"""
const LEMI424_SLOT_INDEX = Dict(c => i for (i, c) in enumerate(LEMI424_NUMERIC_COLUMNS))

"""
Map 1-based whitespace-token position -> record matrix column (0 = not a
numeric data column). A plain `Vector{Int}` so the hot loop indexes an array
instead of hashing a `Dict` key per token.
"""
const LEMI424_TOKEN_SLOT = let table = zeros(Int, length(LEMI424_COLUMNS))
    for (slot, comp) in enumerate(LEMI424_NUMERIC_COLUMNS)
        table[findfirst(==(comp), LEMI424_COLUMNS)] = slot
    end
    table
end

"""
Value written to a record matrix cell whose column is absent from the file.
Mirrors the legacy per-row defaults: counters default to zero, everything else
to `NaN`.
"""
const LEMI424_SLOT_DEFAULT = Float64[
    (c === :n_satellites || c === :gps_fix) ? 0.0 : NaN for c in LEMI424_NUMERIC_COLUMNS
]

const LEMI424_LAT_HEMISPHERE_TOKEN = findfirst(==(:lat_hemisphere), LEMI424_COLUMNS)
const LEMI424_LON_HEMISPHERE_TOKEN = findfirst(==(:lon_hemisphere), LEMI424_COLUMNS)
const LEMI424_MIN_COLUMNS = maximum(values(LEMI424_DEFAULT_COLUMN_INDEX))

const _LEMI_HEMI_N = "N"                                       # interned, one object per file
const _LEMI_HEMI_S = "S"
const _LEMI_HEMI_E = "E"
const _LEMI_HEMI_W = "W"

_lemi424_blank(line::AbstractString) = all(isspace, line)

#---------- record storage -----

"""
    LEMI424Records(times, values, lat_hemisphere, lon_hemisphere)

Column-oriented store for one parsed LEMI-424 file.

Fields:
- `times`: sample timestamps, one per record.
- `values`: `n x length(LEMI424_NUMERIC_COLUMNS)` matrix of raw numeric fields,
  indexed by `LEMI424_SLOT_INDEX`. Latitude/longitude are held here in the raw
  `DDMM.mmmm` form, not decimal degrees.
- `lat_hemisphere`, `lon_hemisphere`: hemisphere letters as written in the file.
"""
struct LEMI424Records
    times::Vector{DateTime}
    values::Matrix{Float64}
    lat_hemisphere::Vector{String}
    lon_hemisphere::Vector{String}
end

"""
    _lemi424_alloc(n) -> LEMI424Records

Allocate a record store for `n` samples with every numeric cell preset to its
absent-column default. Takes the record count; returns the empty store.
"""
function _lemi424_alloc(n::Integer)
    vals = Matrix{Float64}(undef, n, length(LEMI424_NUMERIC_COLUMNS))
    @inbounds for j in eachindex(LEMI424_SLOT_DEFAULT)
        fill!(view(vals, :, j), LEMI424_SLOT_DEFAULT[j])    # column-major fill
    end
    return LEMI424Records(
        Vector{DateTime}(undef, n),
        vals,
        fill(_LEMI_HEMI_N, n),
        fill(_LEMI_HEMI_E, n),
    )
end

#---------- field helpers -----

"""
    _lemi_position(position) -> Float64

Convert a raw LEMI `DDMM.mmmm` coordinate to decimal degrees (unsigned).
Takes the raw numeric field; returns degrees.
"""
function _lemi_position(position)
    scaled = Float64(position) / 100.0
    degrees = floor(Int, scaled)
    minutes = (scaled - degrees) * 100.0
    return degrees + minutes / 60.0
end

"""
    _lemi_hemisphere_sign(h) -> Float64

Sign implied by a hemisphere letter. Takes the letter; returns `-1.0` for south
or west, `+1.0` otherwise. Allocation-free.
"""
function _lemi_hemisphere_sign(h::AbstractString)
    (h == "S" || h == "s" || h == "W" || h == "w") && return -1.0
    return 1.0
end

"""
    _lemi_hemisphere(raw) -> String

Intern a hemisphere token so a whole file shares one `String` per letter.
Takes the raw token; returns the interned letter, or a fresh copy of the token
verbatim when it is not one of N/S/E/W.
"""
function _lemi_hemisphere(raw::AbstractString)
    raw == "N" && return _LEMI_HEMI_N
    raw == "S" && return _LEMI_HEMI_S
    raw == "E" && return _LEMI_HEMI_E
    raw == "W" && return _LEMI_HEMI_W
    return String(raw)                                     # non-canonical, keep as written
end

#---------- line parsing -----

"""
    _parse_lemi424_record!(rec, line, line_number, row) -> nothing

Parse one data line straight into row `row` of `rec`, with no intermediate
allocation. Takes the record store, the line text, its 1-based file line number
(for error messages) and the destination row; mutates `rec` and returns
nothing. Columns past the known layout are ignored; columns absent from the
line keep the defaults set by [`_lemi424_alloc`](@ref).
"""
function _parse_lemi424_record!(rec::LEMI424Records, line::AbstractString,
                                line_number::Integer, row::Integer)
    yr = mo = dy = hr = mi = sec = 0
    lat_hemi = _LEMI_HEMI_N
    lon_hemi = _LEMI_HEMI_E
    ntokens = 0
    for raw in eachsplit(line)                             # lazy, no SubString vector
        ntokens += 1
        if ntokens <= 6
            v = parse(Int, raw)
            ntokens == 1 ? (yr = v) :
            ntokens == 2 ? (mo = v) :
            ntokens == 3 ? (dy = v) :
            ntokens == 4 ? (hr = v) :
            ntokens == 5 ? (mi = v) : (sec = v)
        elseif ntokens > length(LEMI424_COLUMNS)
            break                                          # trailing extras ignored
        elseif ntokens == LEMI424_LAT_HEMISPHERE_TOKEN
            lat_hemi = _lemi_hemisphere(raw)
        elseif ntokens == LEMI424_LON_HEMISPHERE_TOKEN
            lon_hemi = _lemi_hemisphere(raw)
        else
            @inbounds slot = LEMI424_TOKEN_SLOT[ntokens]
            slot == 0 || (@inbounds rec.values[row, slot] = parse(Float64, raw))
        end
    end
    ntokens >= LEMI424_MIN_COLUMNS ||
        error("LEMI-424 line $line_number has $ntokens columns; expected at least $LEMI424_MIN_COLUMNS")

    @inbounds rec.times[row] = DateTime(yr, mo, dy, hr, mi, sec)
    @inbounds rec.lat_hemisphere[row] = lat_hemi
    @inbounds rec.lon_hemisphere[row] = lon_hemi
    return nothing
end

#---------- file reading -----

"""
    _count_lemi424_records(path) -> Int

Count non-blank lines in a LEMI-424 file. Takes the path; returns the record
count and errors when the file holds none.
"""
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

"""
    _fill_lemi424_records!(rec, io) -> Int

Parse every non-blank line of `io` into `rec`. Takes the record store and an
open stream; returns the number of records filled. Kept separate from
[`_read_lemi424_records`](@ref) so the row counter is a plain local rather than
a variable captured and mutated by the `open` closure, which would box it.
"""
function _fill_lemi424_records!(rec::LEMI424Records, io::IO)
    row = 0
    for (line_number, line) in enumerate(eachline(io))
        _lemi424_blank(line) && continue
        row += 1
        _parse_lemi424_record!(rec, line, line_number, row)
    end
    return row
end

"""
    _read_lemi424_records(path) -> LEMI424Records

Read a whole LEMI-424 file in two passes: count records, then parse directly
into preallocated columns. Takes the path; returns the filled record store.
"""
function _read_lemi424_records(path::AbstractString)
    rec = _lemi424_alloc(_count_lemi424_records(path))
    open(io -> _fill_lemi424_records!(rec, io), path, "r")
    return rec
end

"""
    _lemi424_column(rec, component) -> Vector{Float64}

Copy one numeric column out of a record store. Takes the store and a component
name; returns a fresh vector owned by the caller.
"""
_lemi424_column(rec::LEMI424Records, component::Symbol) =
    rec.values[:, LEMI424_SLOT_INDEX[component]]

"""
    _lemi424_geographic(rec) -> (Vector{Float64}, Vector{Float64})

Convert the raw coordinate columns to signed decimal degrees. Takes the record
store; returns `(latitude, longitude)`, leaving non-finite raw values as-is.
"""
function _lemi424_geographic(rec::LEMI424Records)
    lat_slot = LEMI424_SLOT_INDEX[:latitude]
    lon_slot = LEMI424_SLOT_INDEX[:longitude]
    n = length(rec.times)
    lat = Vector{Float64}(undef, n)
    lon = Vector{Float64}(undef, n)
    @inbounds for i in 1:n
        raw_lat = rec.values[i, lat_slot]
        lat[i] = isfinite(raw_lat) ?
            _lemi_position(raw_lat) * _lemi_hemisphere_sign(rec.lat_hemisphere[i]) : raw_lat
        raw_lon = rec.values[i, lon_slot]
        lon[i] = isfinite(raw_lon) ?
            _lemi_position(raw_lon) * _lemi_hemisphere_sign(rec.lon_hemisphere[i]) : raw_lon
    end
    return lat, lon
end

"""
    _lemi424_median(v) -> Float64

Median over the finite entries of `v`. Takes a vector; returns `NaN` when no
entry is finite.
"""
_lemi424_median(v::Vector{Float64}) =
    any(isfinite, v) ? median(filter(isfinite, v)) : NaN

"""
    _lemi424_aux(rec) -> Dict{Symbol, AbstractVector}

Build the auxiliary-column dictionary carried in reader metadata. Takes the
record store; returns numeric columns as `Vector{Float64}` and hemisphere
letters as `Vector{String}` (the element types downstream padding dispatches
on). Coordinates appear here in raw `DDMM.mmmm` form.
"""
function _lemi424_aux(rec::LEMI424Records)
    col(c) = _lemi424_column(rec, c)
    return Dict{Symbol, AbstractVector}(
        :temperature_e => col(:temperature_e),
        :temperature_h => col(:temperature_h),
        :e3 => col(:e3),
        :e4 => col(:e4),
        :battery => col(:battery),
        :elevation => col(:elevation),
        :latitude_raw => col(:latitude),
        :lat_hemisphere => rec.lat_hemisphere,
        :longitude_raw => col(:longitude),
        :lon_hemisphere => rec.lon_hemisphere,
        :n_satellites => col(:n_satellites),
        :gps_fix => col(:gps_fix),
        :time_diff => col(:time_diff),
    )
end

#---------- readers -----

"""
    read_lemi424(path; site, include_aux=true) -> TimekeeperRun

Read a LEMI-424 text file into a run of per-component channels. Takes the file
path, an optional site name and whether to include auxiliary channels
(temperatures, GPS, battery); returns a `TimekeeperRun` at 1 Hz. Channels are
created for every known data column, NaN-filled when the column is absent from
the file.
"""
function read_lemi424(path::AbstractString; site = _site_from_path(path), include_aux = true)
    abs_path = _as_path(path)
    rec = _read_lemi424_records(abs_path)
    latitudes, longitudes = _lemi424_geographic(rec)
    site_str = String(site)

    comps = include_aux ? LEMI424_DATA_COLUMNS : LEMI424_DEFAULT_COMPONENTS
    header = Dict{String, Any}(
        "format" => "LEMI-424",
        "site" => site_str,
        "source_file" => abs_path,
        "sample_rate" => 1.0,
    )
    start = first(rec.times)
    channels = Dict{Symbol, TimekeeperChannel}()
    for comp in comps
        comp in LEMI424_HEMISPHERE_COLUMNS && continue
        data = comp === :latitude ? latitudes :
               comp === :longitude ? longitudes : _lemi424_column(rec, comp)
        channels[comp] = TimekeeperChannel(comp, data, 1.0, start, component_units(comp),
                                           abs_path, copy(header))
    end

    batteries = _lemi424_column(rec, :battery)
    metadata = Dict{Symbol, Any}(
        :source_file => abs_path,
        :sample_rate => 1.0,
        :start_time => start,
        :end_time => last(rec.times),
        :n_samples => length(rec.times),
        :latitude => _lemi424_median(latitudes),
        :longitude => _lemi424_median(longitudes),
        :elevation => _lemi424_median(_lemi424_column(rec, :elevation)),
        :battery_start => first(batteries),
        :battery_end => last(batteries),
        :instrument_model => "LEMI-424",
        :data_logger_manufacturer => "LEMI",
    )
    return TimekeeperRun(site_str, "LEMI-424", :lemi424, channels, metadata)
end

"""
    load_lemi424(path; components, site, include_aux=true) -> TimeArray

Load a LEMI-424 text file as a `TimeArray`. Takes the file path, the components
to keep (any of `LEMI424_DEFAULT_COMPONENTS`), an optional site name and
whether to attach auxiliary columns to the metadata; returns a 1 Hz `TimeArray`
whose metadata carries `:aux_columns` when `include_aux` is set.
"""
function load_lemi424(path::AbstractString; components = LEMI424_DEFAULT_COMPONENTS,
                     site = _site_from_path(path), include_aux = true)
    abs_path = _as_path(path)
    requested = _symbolize.(collect(components))
    missing = [c for c in requested if !(c in LEMI424_DEFAULT_COMPONENTS)]
    isempty(missing) || error("LEMI-424 loader only supports: $(join(LEMI424_DEFAULT_COMPONENTS, ", "))")

    rec = _read_lemi424_records(abs_path)
    n = length(rec.times)
    vals = Matrix{Float64}(undef, n, length(requested))
    @inbounds for (j, comp) in enumerate(requested)
        copyto!(view(vals, :, j), view(rec.values, :, LEMI424_SLOT_INDEX[comp]))
    end

    meta = Dict{Symbol, Any}(
        :site => String(site),
        :instrument => "LEMI-424",
        :source_format => :lemi424,
        :source_file => abs_path,
        :sample_rate => 1.0,
        :start_time => first(rec.times),
        :end_time => last(rec.times),
        :n_samples => n,
        :instrument_model => "LEMI-424",
        :data_logger_manufacturer => "LEMI",
        :reader => :load_lemi424,
        :units => Dict(c => component_units(c) for c in requested),
    )
    include_aux && (meta[:aux_columns] = _lemi424_aux(rec))
    return TimeArray(rec.times, vals, requested, meta)
end

#---------- writers -----

"""
    _lemi_raw_position(decimal_degrees) -> Float64

Convert signed decimal degrees back to the raw LEMI `DDMM.mmmm` form. Takes the
decimal value; returns its unsigned raw encoding (the hemisphere letter carries
the sign in the file).
"""
function _lemi_raw_position(decimal_degrees)
    x = abs(Float64(decimal_degrees))
    degrees = floor(Int, x)
    minutes = (x - degrees) * 60.0
    return degrees * 100.0 + minutes
end

"""
    _value_for_component(run, component, i; default=0.0) -> Float64

Read sample `i` of one channel. Takes the run, the component name, the sample
index and a fallback; returns the sample, or `default` when the channel is
absent or shorter than `i`.
"""
function _value_for_component(run::TimekeeperRun, component::Symbol, i::Integer; default = 0.0)
    haskey(run.channels, component) || return default
    data = run.channels[component].data
    return i <= length(data) ? data[i] : default
end

"""
    _aux_at(aux, sym, i, default)

Read sample `i` of one auxiliary column. Takes the aux dictionary (or
`nothing`), the column name, the sample index and a fallback; returns the
value, or `default` when the column is missing or too short.
"""
_aux_at(::Nothing, ::Symbol, ::Integer, default) = default
function _aux_at(aux::AbstractDict, sym::Symbol, i::Integer, default)
    haskey(aux, sym) || return default
    vec = aux[sym]
    return i <= length(vec) ? vec[i] : default
end

"""
    _lemi_row(run, i, t; aux=nothing) -> String

Format one LEMI-424 text record. Takes the run, the sample index, the
timestamp and optional auxiliary columns; returns the 24-field line. Auxiliary
columns are preferred when present so a read/write round trip preserves the
original GPS and housekeeping fields.
"""
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

"""
    write_lemi424(path, run; timestamps=nothing, aux=nothing) -> String

Write a run as a LEMI-424 text file. Takes the destination path, the run, an
optional timestamp vector (synthesised at 1 Hz when omitted) and optional
auxiliary columns; returns the path. Errors unless the run is 1 Hz. The record
count is the shortest channel.
"""
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

"""
    write_lemi424(path, ta; kwargs...) -> String

Write a `TimeArray` as a LEMI-424 text file. Takes the destination path and the
array, forwarding `kwargs` to `from_timearray`; returns the path. Auxiliary
columns found in the metadata are reused so a round trip keeps the original GPS
and housekeeping fields.
"""
function write_lemi424(path::AbstractString, ta::TimeArray; kwargs...)
    times = collect(_ta_timestamps(ta))
    (isempty(times) || !(first(times) isa DateTime)) &&
        error("LEMI-424 writer requires a non-empty TimeArray with DateTime timestamps")
    md = _ta_meta(ta)
    aux = md isa AbstractDict ? get(md, :aux_columns, nothing) : nothing
    run = from_timearray(ta; instrument = "LEMI-424", source_format = :timearray, kwargs...)
    return write_lemi424(path, run; timestamps = times, aux = aux)
end
