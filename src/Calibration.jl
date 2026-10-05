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
# - a Metronix dipole: its length from the electrode positions in the .ats
#   header
# - a LEMI-424 or GEOMAG dipole: a length that you give, because these files
#   store the electric channels in mV
# - other channels: the field already, in the units of the channel
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
# directories with "cal" in their name in the site directory and in up to three
# directories above it
function _calibration_dirs(site_dir::AbstractString, calibration)
    calibration === nothing || return [abspath(String(d)) for d in (calibration isa AbstractString ? [calibration] : calibration)]
    dirs = String[]
    d = abspath(site_dir)
    for _ in 1:4
        isdir(d) || break
        push!(dirs, d)
        for n in readdir(d)
            full = joinpath(d, n)
            isdir(full) && occursin("cal", lowercase(n)) && push!(dirs, full)
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
    if run.source_format === :metronix && haskey(ch.header, "ats_header_bytes")
        info = _ats_sensor(ch.header["ats_header_bytes"])
        info === nothing && error("Short .ats header in $(ch.source_file)")
        if electric
            dx, dy, dz = info.pos[4] - info.pos[1], info.pos[5] - info.pos[2], info.pos[6] - info.pos[3]
            len = hypot(dx, dy, dz)
            az = len > 0 ? mod(rad2deg(atan(dy, dx)), 360.0) : az0
            haskey(dipole, comp) && (len = dipole[comp])
            len > 0 || (len = info.dipole)
            len > 0 || error("No dipole length for $comp of $(run.site): the .ats header has no electrode positions; give `dipole`")
            haskey(azimuths, comp) && (az = azimuths[comp])
            return _ChannelResponse(:dipole, len / 1000, nothing, false, az,
                                    @sprintf("dipole %.1f m at %.1f°", len, az))
        end
        site_dir = get(run.metadata, :site_dir, dirname(ch.source_file))
        dirs = _calibration_dirs(site_dir, calibration)
        file = find_calibration(dirs, info.sensor, info.serial)
        file === nothing && error("No calibration file for $(info.sensor) $(info.serial) ($comp of $(run.site)). " *
                                  "Searched: $(join(dirs, ", ")). Give the directory with `calibration`")
        cal = _cached_calibration(file, info.sensor, info.serial)
        return _ChannelResponse(:coil, 1.0, cal, info.chopper, az0,
                                "$(info.sensor) $(info.serial), chopper $(info.chopper ? "on" : "off"), $(basename(file))")
    end
    if electric && run.source_format in (:lemi424, :geomag)
        if haskey(dipole, comp)
            return _ChannelResponse(:dipole, dipole[comp] / 1000, nothing, false, az0,
                                    @sprintf("dipole %.1f m", dipole[comp]))
        end
        @warn "No dipole length for $comp of $(run.site) ($(run.source_format) stores mV): " *
              "the transfer function uses 1 km. Give `dipole = Dict(:e1 => L1, :e2 => L2)` in metres" maxlog = 2
        return _ChannelResponse(:dipole, 1.0, nothing, false, az0, "dipole unknown, 1 km")
    end
    return _ChannelResponse(:unit, 1.0, nothing, false, az0, "in $(ch.units)")
end
