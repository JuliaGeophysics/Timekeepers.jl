# Calibration.jl - sensor responses for transfer function estimation.
# Author: @pankajkmishra
#
# A transfer function relates fields in physical units: the electric field in
# mV/km and the magnetic field in nT. The recorders store voltages. This file
# gives, for each channel, the complex factor that changes a Fourier
# coefficient of the record into the field:
# - a Metronix induction coil: its calibration table (a magnitude and a phase
#   for each frequency, with the chopper on and off), from a text file for the
#   type and the serial number of the coil
# - a dipole: a length that you give, or the length from d/dipoles.dat of the
#   survey, or (Metronix) from the electrode positions in the .ats header.
#   Without any of them, each electrode is 50 m from the centre. LEMI-424 and
#   GEOMAG files store the electric channels in mV and need the length
# - other channels: the field already, in the units of the channel
#
# A survey directory holds the coil calibration files in s/ and the dipole
# table in d/dipoles.dat
#
# It also gives the azimuth of each sensor. The processing turns the fields
# to geographic north and east with it

const _ATS_OFF_CHOPPER = 37         # UInt8, 1 = chopper on
const _ATS_OFF_SENSOR_TYPE = 40     # 6 characters, e.g. "MFS07e"
const _ATS_OFF_SENSOR_SERIAL = 46   # Int16
const _ATS_OFF_POSITIONS = 48       # six Float32: x1 y1 z1 x2 y2 z2, metres, x north
const _ATS_OFF_DIPOLE = 72          # Float32, metres
const _ATS_OFF_ANGLE = 76           # Float32, degrees east of north

"""
    SensorCalibration

The calibration of one induction coil. A Metronix table gives the response
normalized by frequency: `magnitude` in V/(nT·Hz) and `phase` in degrees. Thus,
the response of the coil is `magnitude · f · exp(i·phase)` in V/nT. The file
has one table with the chopper on and one with the chopper off.

# Fields
- `sensor::String` -- the sensor type, e.g. `"MFS07e"`.
- `serial::Int` -- the serial number.
- `path::String` -- the calibration file.
- `on`, `off` -- `(f, magnitude, phase)` tables in increasing frequency, empty
  if the file has no table for that state.

See also [`read_calibration`](@ref), [`sensor_response`](@ref).
"""
struct SensorCalibration
    sensor::String
    serial::Int
    path::String
    on::NTuple{3, Vector{Float64}}
    off::NTuple{3, Vector{Float64}}
end

function Base.show(io::IO, c::SensorCalibration)
    span(t) = isempty(t[1]) ? "none" : @sprintf("%.3g–%.3g Hz", first(t[1]), last(t[1]))
    print(io, "SensorCalibration(\"$(c.sensor)\", $(c.serial), chopper on: $(span(c.on)), off: $(span(c.off)))")
end

const _EMPTY_TABLE = (Float64[], Float64[], Float64[])

"""
    read_calibration(path; sensor = "", serial = 0) -> SensorCalibration

Read a Metronix calibration text file. It has a `Chopper On` and a `Chopper
Off` block, each with rows of frequency, magnitude and phase.

The function examines each table. A coil has a smooth response. Thus, a jump
by a factor of more than 5 in the normalized magnitude from one row to the next
is a part of a table that does not belong to the coil (some files hold rows of
a different model at their low end). The function keeps the longest smooth
part of each table and gives a warning about the rows that it removed.
"""
function read_calibration(path::AbstractString; sensor::AbstractString = "", serial::Integer = 0)
    tables = Dict{Symbol, Vector{NTuple{3, Float64}}}(:on => NTuple{3, Float64}[], :off => NTuple{3, Float64}[])
    state = :on
    for line in eachline(path)
        s = lowercase(strip(line))
        occursin("chopper on", s) && (state = :on; continue)
        occursin("chopper off", s) && (state = :off; continue)
        parts = split(s)
        length(parts) == 3 || continue
        vals = tryparse.(Float64, parts)
        any(isnothing, vals) && continue
        f, m, p = vals
        (f > 0 && m > 0) || continue
        push!(tables[state], (f, m, p))
    end
    isempty(tables[:on]) && isempty(tables[:off]) && error("No calibration rows in $path")
    clean(rows, label) = _smooth_table(rows, path, label)
    return SensorCalibration(String(sensor), Int(serial), abspath(path),
                             clean(tables[:on], "chopper on"), clean(tables[:off], "chopper off"))
end

# Sort a table by frequency, merge repeated frequencies and keep the longest
# part without a jump in the normalized magnitude
function _smooth_table(rows::Vector{NTuple{3, Float64}}, path, label)
    isempty(rows) && return _EMPTY_TABLE
    sort!(rows; by = first)
    rows = [rows[i] for i in eachindex(rows) if i == 1 || rows[i][1] > rows[i - 1][1] * (1 + 1.0e-9)]
    parts = UnitRange{Int}[]
    lo = 1
    for i in 2:length(rows)
        r = rows[i][2] / rows[i - 1][2]
        if r > 5 || r < 1 / 5
            push!(parts, lo:(i - 1))
            lo = i
        end
    end
    push!(parts, lo:length(rows))
    keep = parts[argmax(length.(parts))]
    if length(keep) < length(rows)
        dropped = length(rows) - length(keep)
        @warn "Calibration $(basename(path)) ($label): removed $dropped rows that jump from the smooth response" maxlog = 4
    end
    sel = rows[keep]
    return (first.(sel), getindex.(sel, 2), last.(sel))
end

# Linear interpolation in log frequency of a log magnitude and a phase. Out of
# the table the end value stays
function _interp_table(t, f::Real)
    fs, ms, ps = t
    lf = log(f)
    f <= fs[1] && return (ms[1], ps[1], false)
    f >= fs[end] && return (ms[end], ps[end], f > fs[end] * (1 + 1.0e-9))
    j = searchsortedlast(fs, f)
    a = (lf - log(fs[j])) / (log(fs[j + 1]) - log(fs[j]))
    m = exp((1 - a) * log(ms[j]) + a * log(ms[j + 1]))
    p = (1 - a) * ps[j] + a * ps[j + 1]
    return (m, p, true)
end

_covers(t, f) = !isempty(t[1]) && t[1][1] * (1 - 1.0e-9) <= f <= t[1][end] * (1 + 1.0e-9)

"""
    sensor_response(cal::SensorCalibration, f; chopper = true) -> ComplexF64

The response of the coil at the frequency `f` in Hz, in V/nT. The function
uses the table for the chopper state. If that table does not cover `f` but the
other table does, the function uses the other table. Out of both tables, the
normalized response keeps its value at the nearest end. That is the response
of a coil far below its corner frequency, where its output is proportional to
`f`.
"""
function sensor_response(cal::SensorCalibration, f::Real; chopper::Bool = true)
    want, other = chopper ? (cal.on, cal.off) : (cal.off, cal.on)
    t = !_covers(want, f) && _covers(other, f) ? other :
        isempty(want[1]) ? other : want
    m, p, _ = _interp_table(t, f)
    return m * f * cis(deg2rad(p))
end

#---------- finding the calibration files -----

const _CALIBRATION_CACHE = Dict{String, SensorCalibration}()

function _cached_calibration(path::AbstractString, sensor, serial)
    key = abspath(path)
    return get!(_CALIBRATION_CACHE, key) do
        read_calibration(key; sensor = sensor, serial = serial)
    end
end

# The directories to search for calibration files: the ones that you give, or
# the s/ directory of the survey and directories with "cal" in their name, in
# the site directory and in up to three directories above it
function _calibration_dirs(site_dir::AbstractString, calibration)
    calibration === nothing || return [abspath(String(d)) for d in (calibration isa AbstractString ? [calibration] : calibration)]
    dirs = String[]
    d = abspath(site_dir)
    for _ in 1:4
        isdir(d) || break
        push!(dirs, d)
        for n in readdir(d)
            full = joinpath(d, n)
            isdir(full) && (n == "s" || occursin("cal", lowercase(n))) && push!(dirs, full)
        end
        parent = dirname(d)
        parent == d && break
        d = parent
    end
    return unique(dirs)
end

"""
    find_calibration(dirs, sensor, serial) -> Union{String, Nothing}

The calibration file of a coil in `dirs`. The name of a Metronix file is the
sensor type and the serial number in three digits, e.g. `MFS07e160.txt`. If no
file has that name, a file with the same serial number and a sensor family
that agrees (the type of an unknown sensor agrees with each family) is the
file, but only if it is the only one.
"""
function find_calibration(dirs, sensor::AbstractString, serial::Integer)
    number = lpad(serial, 3, '0')
    family = lowercase(strip(sensor))
    want = family * number * ".txt"
    unknown = isempty(family) || startswith(family, "unkn")
    candidates = String[]
    for d in dirs
        isdir(d) || continue
        for n in readdir(d)
            full = joinpath(d, n)
            isfile(full) || continue
            ln = lowercase(n)
            ln == want && return full
            stem, ext = splitext(ln)
            (ext == ".txt" && endswith(stem, number)) || continue
            prefix = stem[1:(end - length(number))]
            occursin(r"^[a-z]+\d*[a-z]?$", prefix) || continue
            (unknown || startswith(prefix, family) || startswith(family, prefix)) && push!(candidates, full)
        end
    end
    unique!(candidates)
    return length(candidates) == 1 ? candidates[1] : nothing
end

#---------- the dipole table of a survey -----

# The distance of each electrode from the centre of a site when no header and
# no table gives it
const DEFAULT_ELECTRODE_DISTANCE = 50.0

"""
    read_dipoles(path) -> Dict{String, NTuple{4, Float64}}

Read a dipole table, `d/dipoles.dat` of a survey. The first line is the
header. Each other line holds a site name and the distances in metres from the
centre of the site to the N, S, E and W electrodes:

    site       N      S      E      W
    site002   24.5   48.5   50.9   51.5

Ex is N + S long and Ey is E + W. `-` is a distance that is not known. The
result maps a site name to `(N, S, E, W)`, with `NaN` for `-`.
"""
function read_dipoles(path::AbstractString)
    table = Dict{String, NTuple{4, Float64}}()
    header = true
    for (i, line) in enumerate(readlines(path))
        t = split(line)
        isempty(t) && continue
        length(t) == 5 || error("$path, line $i: give a site and four distances (N S E W), got \"$(strip(line))\"")
        if header
            header = false
            uppercase.(t[2:5]) == ["N", "S", "E", "W"] && continue
        end
        d = map(t[2:5]) do x
            x == "-" && return NaN
            v = tryparse(Float64, x)
            (v === nothing || v < 0) && error("$path, line $i: \"$x\" is not a distance in metres")
            v
        end
        table[String(t[1])] = Tuple(d)
    end
    return table
end

# d/dipoles.dat in the site directory or in up to three directories above it
function _dipole_file(site_dir::AbstractString)
    d = abspath(site_dir)
    for _ in 1:4
        f = joinpath(d, "d", "dipoles.dat")
        isfile(f) && return f
        parent = dirname(d)
        parent == d && break
        d = parent
    end
    return nothing
end

# The row of a site in the dipole table of its survey: (N, S, E, W), the file
# and the name of the row, or `nothing` for the row if the table has none. The
# row of a rate directory "site.128" is the row of "site"
function _table_row(site_dir::AbstractString, names)
    file = _dipole_file(site_dir)
    file === nothing && return nothing, nothing, ""
    table = read_dipoles(file)
    for n in names, key in (n, replace(n, r"\.\d+(?:\.\d+)?$" => ""))
        haskey(table, key) && return table[key], file, key
    end
    return nothing, file, ""
end

# The length in metres of the dipole `comp` of a site from the dipole table
# of its survey: (length, file), the length NaN if the table does not give it
function _table_dipole(site_dir::AbstractString, names, comp::Symbol)
    row, file, _ = _table_row(site_dir, names)
    row === nothing && return NaN, file
    return (comp === :e1 ? row[1] + row[2] : row[3] + row[4]), file
end

# The length in metres of an electric channel and where it comes from, first
# found first: the dipole that you give, the row of the site in d/dipoles.dat
# of the survey (a distance not known there is 50 m), the length from the
# header (`header_len`, 0 if none), or 50 m from the centre to each electrode
function _dipole_length(run::TimekeeperRun, comp::Symbol, site_dir::AbstractString, dipole, header_len::Real)
    haskey(dipole, comp) && return dipole[comp], "given to the processing", ""
    row, file, key = _table_row(site_dir, unique([run.site, basename(rstrip(site_dir, ['/', '\\']))]))
    if row !== nothing
        sides = comp === :e1 ? (("N", row[1]), ("S", row[2])) : (("E", row[3]), ("W", row[4]))
        unknown = [name for (name, d) in sides if !isfinite(d)]
        isempty(unknown) ||
            @warn "$(file) gives no $(join(unknown, ", ")) distance for $key: it is $(DEFAULT_ELECTRODE_DISTANCE) m" maxlog = 4
        ds = [isfinite(d) ? d : DEFAULT_ELECTRODE_DISTANCE for (_, d) in sides]
        text = join((@sprintf("%s %.1f m%s", name, d, isfinite(raw) ? "" : " (default)")
                     for ((name, raw), d) in zip(sides, ds)), " + ")
        return sum(ds), "$text from the centre, row $key of $(file)", file
    end
    header_len > 0 && return Float64(header_len), "from the electrode positions in the .ats header", ""
    len = 2 * DEFAULT_ELECTRODE_DISTANCE
    where = file === nothing ? "no d/dipoles.dat in the survey" : "no row for $(run.site) in $(file)"
    @warn "No dipole length for $comp of $(run.site) ($where): each electrode is " *
          "$(DEFAULT_ELECTRODE_DISTANCE) m from the centre, $(len) m" maxlog = 4
    return len, "default, $(DEFAULT_ELECTRODE_DISTANCE) m from the centre to each electrode ($where)", ""
end

#---------- channel responses -----

"""
    _ChannelResponse

The response of one channel for the processing. A Fourier coefficient of the
record divided by `response(f)` is the field in physical units (mV/km or nT).
`azimuth` is the direction of the sensor in degrees east of north.
"""
struct _ChannelResponse
    kind::Symbol              # :coil, :dipole or :unit
    gain::Float64             # dipole length in km, or a fixed factor
    cal::Union{Nothing, SensorCalibration}
    chopper::Bool
    azimuth::Float64
    note::String
    file::String              # the calibration file or the dipole table, "" if none
end

function _response(r::_ChannelResponse, f::Real)
    r.kind === :coil && return 1000.0 * sensor_response(r.cal, f; chopper = r.chopper)  # mV per nT
    return complex(r.gain)
end

# The sensor fields of an .ats header
function _ats_sensor(hbytes::Vector{UInt8})
    length(hbytes) >= _ATS_OFF_ANGLE + 4 || return nothing
    raw = hbytes[(_ATS_OFF_SENSOR_TYPE + 1):(_ATS_OFF_SENSOR_TYPE + 6)]
    sensor = strip(String(filter(b -> b != 0x00, raw)))
    pos = [Float64(_ats_get(Float32, hbytes, _ATS_OFF_POSITIONS + 4k)) for k in 0:5]
    return (chopper = hbytes[_ATS_OFF_CHOPPER + 1] == 0x01,
            sensor = String(sensor),
            serial = Int(_ats_get(Int16, hbytes, _ATS_OFF_SENSOR_SERIAL)),
            pos = pos,
            dipole = Float64(_ats_get(Float32, hbytes, _ATS_OFF_DIPOLE)),
            angle = Float64(_ats_get(Float32, hbytes, _ATS_OFF_ANGLE)))
end

const _NOMINAL_AZIMUTH = Dict(:e1 => 0.0, :e2 => 90.0, :bx => 0.0, :by => 90.0, :bz => 0.0)

# The response of one channel of a run. `dipole` holds lengths in metres for
# the electric channels of text formats, `azimuths` the directions that you
# give
function _channel_response(run::TimekeeperRun, comp::Symbol; calibration = nothing,
                           dipole = Dict{Symbol, Float64}(), azimuths = Dict{Symbol, Float64}())
    ch = run.channels[comp]
    az0 = get(azimuths, comp, get(_NOMINAL_AZIMUTH, comp, 0.0))
    electric = comp in (:e1, :e2)
    site_dir = String(get(run.metadata, :site_dir, dirname(ch.source_file)))
    if run.source_format === :metronix && haskey(ch.header, "ats_header_bytes")
        info = _ats_sensor(ch.header["ats_header_bytes"])
        info === nothing && error("Short .ats header in $(ch.source_file)")
        if electric
            dx, dy, dz = info.pos[4] - info.pos[1], info.pos[5] - info.pos[2], info.pos[6] - info.pos[3]
            header_len = hypot(dx, dy, dz)
            az = header_len > 0 ? mod(rad2deg(atan(dy, dx)), 360.0) : az0
            header_len > 0 || (header_len = info.dipole)
            len, source, table = _dipole_length(run, comp, site_dir, dipole, header_len)
            haskey(azimuths, comp) && (az = azimuths[comp])
            return _ChannelResponse(:dipole, len / 1000, nothing, false, az,
                                    @sprintf("dipole %.1f m at %.1f°, %s", len, az, source), table)
        end
        dirs = _calibration_dirs(site_dir, calibration)
        file = find_calibration(dirs, info.sensor, info.serial)
        file === nothing && error("No calibration file for $(info.sensor) $(info.serial) ($comp of $(run.site)). " *
                                  "Searched: $(join(dirs, ", ")). Put the file in s/ of the survey or give the directory with `calibration`")
        cal = _cached_calibration(file, info.sensor, info.serial)
        return _ChannelResponse(:coil, 1.0, cal, info.chopper, az0,
                                "$(info.sensor) $(info.serial), chopper $(info.chopper ? "on" : "off"), calibration $(file)", file)
    end
    if electric && run.source_format in (:lemi424, :geomag)
        len, source, table = _dipole_length(run, comp, site_dir, dipole, 0.0)
        return _ChannelResponse(:dipole, len / 1000, nothing, false, az0, @sprintf("dipole %.1f m, %s", len, source), table)
    end
    return _ChannelResponse(:unit, 1.0, nothing, false, az0, "in $(ch.units)", "")
end
