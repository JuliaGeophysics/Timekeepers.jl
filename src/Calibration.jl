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
#   Without any of them, each electrode is 50 m from the centre. GEOMAG files
#   store the electric channels in mV. LEMI-424 files store them in µV: the
#   logger writes µV/m for the line length L set on it, and L is taken as 1 m
#   (the default of the logger). Both need the dipole length
# - other channels: the field already, in the units of the channel
#
# A survey directory holds the coil calibration files in s/, the dipole table
# in d/dipoles.dat and the sensor orientations in d/orientations.dat
#
# It also gives the azimuth of each sensor, from d/orientations.dat when the
# survey has one (also a reversed or swapped channel: see
# `read_orientations`). The processing turns the fields to geographic north
# and east with it

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
    isempty(strip(site_dir)) && return nothing
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

# The names under which a site can have a row in a table of its survey: the
# name of the run, the name of its directory, and the name of the site that the
# survey scan gives it (the folder above a LEMI card folder, see
# `_site_name_from_dir`)
_table_names(run::TimekeeperRun, site_dir::AbstractString) =
    unique([run.site, basename(rstrip(site_dir, ['/', '\\'])), _site_name_from_dir(site_dir)])

# The line (:ns or :ew) that a channel at `azimuth` records
_line_of(azimuth::Real) = abs(cosd(azimuth)) >= abs(sind(azimuth)) ? :ns : :ew

# The length in metres of the dipole `comp` of a site from the dipole table
# of its survey: (length, file), the length NaN if the table does not give it.
# `line` is the line that the channel records (:ns = N + S, :ew = E + W)
function _table_dipole(site_dir::AbstractString, names, comp::Symbol; line::Symbol = comp === :e1 ? :ns : :ew)
    row, file, _ = _table_row(site_dir, names)
    row === nothing && return NaN, file
    return (line === :ns ? row[1] + row[2] : row[3] + row[4]), file
end

# The length in metres of an electric channel and where it comes from, first
# found first: the dipole that you give, the row of the site in d/dipoles.dat
# of the survey (a distance not known there is 50 m), the length from the
# header (`header_len`, 0 if none), or 50 m from the centre to each electrode.
# `line` is the line that the channel records: E1 records the N-S line unless
# d/orientations.dat turns it to the E-W line
function _dipole_length(run::TimekeeperRun, comp::Symbol, site_dir::AbstractString, dipole, header_len::Real;
                        line::Symbol = comp === :e1 ? :ns : :ew)
    haskey(dipole, comp) && return dipole[comp], "given to the processing", ""
    row, file, key = _table_row(site_dir, _table_names(run, site_dir))
    if row !== nothing
        sides = line === :ns ? (("N", row[1]), ("S", row[2])) : (("E", row[3]), ("W", row[4]))
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

#---------- the orientation table of a survey -----

"""
    read_orientations(path) -> (declination, rows)

Read the orientation table of a survey, `d/orientations.dat`: the
declination, and for each site the polarity of each channel.

    # comments start with #
    declination 10.49
    site   Ex   Ey   Hx   Hy
    1      +1   +1   +1   +1
    2      +1   -1   +1   +1

- `+1` is a channel wired as planned, `-1` a channel wired the other way round
  in the field: the processing reverses it. [`flip_channel!`](@ref) changes
  the sign. `-` is `+1`.
- Ex is the N-S line (red at geographic north), Ey the E-W line (red at
  east).
- The x axis of a fluxgate (LEMI-424, GEOMAG) comes from its static field:
  `declination − atan2(By, Bx)`; y is 90° more and z is down. Metronix coils
  are taken as Hx north and Hy east, and Metronix electric lines in the
  directions of the electrode positions of their .ats headers.

The result is the declination (`NaN` if the file has none; the processing
then takes IGRF-14 at the position of the run) and a `Dict` from a site name
to its row as text (`Dict{Symbol, String}` keyed by the channels :e1 (Ex),
:e2 (Ey), :bx (Hx), :by (Hy)). Hz has no column: a vertical sensor cannot
be wired the other way round. See also [`write_orientations`](@ref),
[`flip_channel!`](@ref).
"""
function read_orientations(path::AbstractString)
    decl = NaN
    rows = Dict{String, Dict{Symbol, String}}()
    cols = Symbol[]
    for (i, raw) in enumerate(readlines(path))
        line = strip(replace(raw, r"#.*$" => ""))
        isempty(line) && continue
        t = split(line)
        if lowercase(t[1]) == "declination"
            length(t) == 2 && (v = tryparse(Float64, t[2])) !== nothing ||
                error("$path, line $i: write the declination as \"declination 10.5\"")
            decl = v
        elseif lowercase(t[1]) == "site"
            names = join(values(ORIENTATION_NAMES), " ")
            cols = [get(_ORIENTATION_BY_NAME, lowercase(n), :bad) for n in t[2:end]]
            all(!=(:bad), cols) || error("$path, line $i: the columns are site $names")
        else
            isempty(cols) && error("$path, line $i: the header line (site $(join(values(ORIENTATION_NAMES), " "))) must come before the sites")
            length(t) == length(cols) + 1 ||
                error("$path, line $i: give a value for each of $(join((ORIENTATION_NAMES[c] for c in cols), ", "))")
            row = _new_orientation_row()
            for (c, v) in zip(cols, t[2:end])
                _parse_polarity(v) === :bad && error("$path, line $i: \"$v\" is not +1 or -1 ($(ORIENTATION_NAMES[c]))")
                row[c] = String(v)
            end
            rows[String(t[1])] = row
        end
    end
    return decl, rows
end

const ORIENTATION_COLUMNS = (:e1, :e2, :bx, :by)
# the names of the channels in the table and on the screen
const ORIENTATION_NAMES = (e1 = "Ex", e2 = "Ey", bx = "Hx", by = "Hy")
const _ORIENTATION_BY_NAME = Dict(lowercase(v) => k for (k, v) in pairs(ORIENTATION_NAMES))

_new_orientation_row() = Dict{Symbol, String}(c => "+1" for c in ORIENTATION_COLUMNS)

# +1 or -1 from the text of the table ("-" is +1), or :bad
function _parse_polarity(v::AbstractString)
    s = strip(v)
    s in ("+1", "1", "-") && return 1.0
    s == "-1" && return -1.0
    return :bad
end

"""
    write_orientations(path, declination, rows)

Write an orientation table (see [`read_orientations`](@ref)), with the sites in
the order of their names. `declination` can be `NaN` (no line).
"""
function write_orientations(path::AbstractString, declination::Real, rows::AbstractDict)
    mkpath(dirname(abspath(path)))
    open(path, "w") do io
        println(io, "# Polarity of the channels of each site (Timekeepers.jl, read_orientations).")
        println(io, "# +1 wired as planned, -1 wired the other way round (the processing reverses it).")
        println(io, "# Ex = N-S line, Ey = E-W line, geographic. Fluxgate axes from the static field")
        println(io, "# and the declination; Metronix coils Hx north, Hy east, lines from the .ats electrode positions.")
        isfinite(declination) && @printf(io, "declination %.2f\n", declination)
        println(io, rpad("site", 10), join((rpad(ORIENTATION_NAMES[c], 6) for c in ORIENTATION_COLUMNS)))
        for name in sort(collect(keys(rows)); by = n -> (something(tryparse(Int, n), typemax(Int)), n))
            r = rows[name]
            println(io, rpad(name, 10), " ", join((rpad(get(r, c, "+1"), 6) for c in ORIENTATION_COLUMNS)))
        end
    end
    return path
end

"""
    orientation_file(site_dir) -> Union{String, Nothing}

`d/orientations.dat` in `site_dir` or in up to three directories above it.
"""
function orientation_file(site_dir::AbstractString)
    # a record without a directory (made in memory) has no survey
    isempty(strip(site_dir)) && return nothing
    d = abspath(site_dir)
    for _ in 1:4
        f = joinpath(d, "d", "orientations.dat")
        isfile(f) && return f
        parent = dirname(d)
        parent == d && break
        d = parent
    end
    return nothing
end

_wrap180(a) = (x = mod(a + 180.0, 360.0) - 180.0; x == -180.0 ? 180.0 : x)
_fmt_az(a) = @sprintf("%.1f", _wrap180(a))

"""
    flip_channel!(path, site, comp) -> Dict{Symbol, String}

Change the polarity of one channel (:e1 Ex, :e2 Ey, :bx Hx or :by Hy) of `site` in
the orientation table `path` (made if it is not there), `+1` ↔ `-1`, and
return the new row.
"""
function flip_channel!(path::AbstractString, site::AbstractString, comp::Symbol)
    comp in ORIENTATION_COLUMNS || error("flip_channel!: no channel $comp")
    decl, rows = isfile(path) ? read_orientations(path) : (NaN, Dict{String, Dict{Symbol, String}}())
    row = get!(_new_orientation_row, rows, String(site))
    row[comp] = _parse_polarity(row[comp]) == -1.0 ? "+1" : "-1"
    write_orientations(path, decl, rows)
    return row
end

# The direction of the x axis of a fluxgate from its static field: the field
# is at the declination from true north, and at atan2(By, Bx) from x
function _auto_azimuth(run::TimekeeperRun, declination::Real; sx::Real = 1.0, sy::Real = 1.0)
    (haskey(run.channels, :bx) && haskey(run.channels, :by) && isfinite(declination)) || return NaN
    X, Y = get(run.metadata, :total_field_x, nothing), get(run.metadata, :total_field_y, nothing)
    X isa Real && Y isa Real && hypot(X, Y) > 0 && return _wrap180(declination - atand(sy * Y, sx * X))
    bx = filter(isfinite, run.channels[:bx].data); by = filter(isfinite, run.channels[:by].data)
    (isempty(bx) || isempty(by)) && return NaN
    return _wrap180(declination - atand(sy * median(by), sx * median(bx)))
end

_has_static_field(run::TimekeeperRun) = run.source_format in (:lemi424, :geomag)

"""
    site_orientation(run; site_dir) -> NamedTuple

The orientation of the channels of `run`: `(e1, e2, bx, by, signs, notes,
file)`. The azimuths are in degrees east of true north: for
a fluxgate the x axis from its static field and the declination (`NaN` keeps
the header or nominal direction). `signs[c]` is the polarity of each channel
from the orientation table of the survey (refer to
[`read_orientations`](@ref)), `notes` a short note for each channel and
`file` the table (`""` if none).
"""
function site_orientation(run::TimekeeperRun; site_dir::AbstractString = _run_site_dir(run))
    signs = Dict{Symbol, Float64}(c => 1.0 for c in ORIENTATION_COLUMNS)
    notes = Dict{Symbol, String}()
    file = orientation_file(site_dir)
    decl, rows = file === nothing ? (NaN, Dict{String, Dict{Symbol, String}}()) : read_orientations(file)
    key = findfirst(n -> haskey(rows, n), _table_names(run, site_dir))
    if key !== nothing
        row = rows[_table_names(run, site_dir)[key]]
        for c in ORIENTATION_COLUMNS
            signs[c] = _parse_polarity(row[c])
            signs[c] < 0 && (notes[c] = "flipped (-1)")
        end
    end
    bx = by = NaN
    if _has_static_field(run)
        # a table without a declination: IGRF-14 at the run's position and time
        decl_note = "declination"
        if !isfinite(decl)
            lat, lon = get(run.metadata, :latitude, NaN), get(run.metadata, :longitude, NaN)
            if lat isa Real && lon isa Real && isfinite(lat) && isfinite(lon)
                decl = try
                    igrf_declination(lat, lon, start_time(run))
                catch
                    NaN
                end
                decl_note = "IGRF-14 declination"
            end
        end
        a = _auto_azimuth(run, decl; sx = signs[:bx], sy = signs[:by])
        if isfinite(a)
            bx, by = a, _wrap180(a + 90)
            for (c, v) in ((:bx, bx), (:by, by))
                notes[c] = @sprintf("%.1f° (%s %.2f° and the static field)", v, decl_note, decl) *
                           (signs[c] < 0 ? ", flipped (-1)" : "")
            end
        end
    end
    return (; e1 = NaN, e2 = NaN, bx, by, signs, notes, file = something(file, ""))
end

_run_site_dir(run::TimekeeperRun) =
    String(get(run.metadata, :site_dir,
               isempty(run.channels) ? "" : dirname(first(values(run.channels)).source_file)))

#---------- channel responses -----

"""
    _ChannelResponse

The response of one channel for the processing. A Fourier coefficient of the
record divided by `response(f)` is the field in physical units (mV/km or nT).
`azimuth` is the direction of the sensor in degrees east of north. `scale`
multiplies the record before that: the unit of the record (1e-3 for the µV
of a LEMI-424 line), the sign of a line voltage (`_E_POLARITY`) and the
polarity of the channel (-1: wired the other way round).
"""
struct _ChannelResponse
    kind::Symbol              # :coil, :dipole or :unit
    gain::Float64             # dipole length in km, or a fixed factor
    cal::Union{Nothing, SensorCalibration}
    chopper::Bool
    azimuth::Float64
    note::String
    file::String              # the calibration file or the dipole table, "" if none
    scale::Float64
    delay::Float64            # seconds that the record lags the field (filters of the logger)
end
_ChannelResponse(kind, gain, cal, chopper, azimuth, note, file) =
    _ChannelResponse(kind, gain, cal, chopper, azimuth, note, file, 1.0, 0.0)
_ChannelResponse(kind, gain, cal, chopper, azimuth, note, file, scale) =
    _ChannelResponse(kind, gain, cal, chopper, azimuth, note, file, scale, 0.0)

function _response(r::_ChannelResponse, f::Real)
    lag = r.delay == 0 ? complex(1.0) : cis(-2π * f * r.delay)
    r.kind === :coil && return 1000.0 * sensor_response(r.cal, f; chopper = r.chopper) * lag / r.scale  # mV per nT
    return r.gain * lag / r.scale
end

# The delay of each channel of a LEMI-424 against the time of its record, from
# the typical response in appendix 3 of the LEMI-424 manual (flat to 0.1 %
# below 0.05 Hz; the fluxgate channels lag the electric ones by about 32 ms,
# 2.9° of impedance phase at 4 s)
const _LEMI424_DELAY = Dict(:bx => 0.031658, :by => 0.033291, :bz => 0.031442,
                            :e1 => -0.000899, :e2 => -0.000873, :e3 => -0.000913, :e4 => -0.000863)

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

# The unit of the electric channels of a text format, as a factor to mV: the
# LEMI-424 logger writes µV/m for its line length L, and L is taken as 1 m
const _E_TO_MV = Dict(:lemi424 => 1.0e-3, :geomag => 1.0)

const _E_POLARITY = -1.0

# The response of one channel of a run. `dipole` holds lengths in metres for
# the electric channels of text formats, `azimuths` the directions that you
# give. The direction of a channel, first found first: `azimuths`, the row of
# the site in d/orientations.dat of the survey, the .ats header (Metronix), the
# nominal direction. A channel with polarity -1 in the table (wired the other
# way round) gets a negative scale
function _channel_response(run::TimekeeperRun, comp::Symbol; calibration = nothing,
                           dipole = Dict{Symbol, Float64}(), azimuths = Dict{Symbol, Float64}(), orientation = nothing)
    ch = run.channels[comp]
    site_dir = String(get(run.metadata, :site_dir, dirname(ch.source_file)))
    orient = comp in ORIENTATION_COLUMNS ? (orientation === nothing ? site_orientation(run; site_dir) : orientation) :
             (; e1 = NaN, e2 = NaN, bx = NaN, by = NaN, signs = Dict{Symbol, Float64}(), notes = Dict{Symbol, String}(), file = "")
    table_az = comp === :bz ? NaN : get(orient, comp, NaN)
    from_table(a) = isfinite(table_az) ? table_az : a
    az0 = get(azimuths, comp, from_table(get(_NOMINAL_AZIMUTH, comp, 0.0)))
    sign = get(orient.signs, comp, 1.0)          # -1: wired the other way round in the field
    orient_note = haskey(azimuths, comp) ? "" :
                  haskey(orient.notes, comp) ? ", orientation $(orient.notes[comp]) from $(orient.file)" : ""
    electric = comp in (:e1, :e2)
    if run.source_format === :metronix && haskey(ch.header, "ats_header_bytes")
        info = _ats_sensor(ch.header["ats_header_bytes"])
        info === nothing && error("Short .ats header in $(ch.source_file)")
        if electric
            dx, dy, dz = info.pos[4] - info.pos[1], info.pos[5] - info.pos[2], info.pos[6] - info.pos[3]
            header_len = hypot(dx, dy, dz)
            az = header_len > 0 ? mod(rad2deg(atan(dy, dx)), 360.0) : az0
            header_len > 0 || (header_len = info.dipole)
            az = haskey(azimuths, comp) ? azimuths[comp] : from_table(az)
            len, source, table = _dipole_length(run, comp, site_dir, dipole, header_len; line = _line_of(az))
            return _ChannelResponse(:dipole, len / 1000, nothing, false, az,
                                    @sprintf("dipole %.1f m at %.1f°, %s, E = -V/L", len, az, source) * orient_note, table,
                                    _E_POLARITY * sign)
        end
        dirs = _calibration_dirs(site_dir, calibration)
        file = find_calibration(dirs, info.sensor, info.serial)
        file === nothing && error("No calibration file for $(info.sensor) $(info.serial) ($comp of $(run.site)). " *
                                  "Searched: $(join(dirs, ", ")). Put the file in s/ of the survey or give the directory with `calibration`")
        cal = _cached_calibration(file, info.sensor, info.serial)
        return _ChannelResponse(:coil, 1.0, cal, info.chopper, az0,
                                "$(info.sensor) $(info.serial), chopper $(info.chopper ? "on" : "off"), calibration $(file)" *
                                orient_note, file, sign)
    end
    if electric && haskey(_E_TO_MV, run.source_format)
        len, source, table = _dipole_length(run, comp, site_dir, dipole, 0.0; line = _line_of(az0))
        unit = run.source_format === :lemi424 ? "µV (LEMI-424, L = 1 m on the logger)" : "mV"
        return _ChannelResponse(:dipole, len / 1000, nothing, false, az0,
                                @sprintf("dipole %.1f m at %.1f°, %s, record in %s, E = -V/L", len, az0, source, unit) * orient_note,
                                table, _E_POLARITY * sign * _E_TO_MV[run.source_format],
                                run.source_format === :lemi424 ? get(_LEMI424_DELAY, comp, 0.0) : 0.0)
    end
    delay = run.source_format === :lemi424 ? get(_LEMI424_DELAY, comp, 0.0) : 0.0
    return _ChannelResponse(:unit, 1.0, nothing, false, az0,
                            "in $(ch.units)" * (delay == 0 ? "" : @sprintf(", LEMI-424 delay %.1f ms removed", 1000delay)) * orient_note,
                            "", sign, delay)
end
