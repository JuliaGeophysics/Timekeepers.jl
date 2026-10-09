# Survey.jl - survey scan and selection of reference sites.
# Author: @pankajkmishra
#
# A survey is a directory of sites that recorded in the same weeks. For each
# site, the sites that recorded at the same time and rate are its base sites
# when they are near, and its remote sites when they are far. To find them,
# the file needs only when, how fast and where each site recorded. Thus, it
# reads only the headers: the .ats headers of a Metronix site, and the first
# and the last lines of LEMI-424 and GEOMAG files.
#
# From this index, the file measures the time that two sites recorded together
# at a common rate, and the distance between them. It puts the base and remote
# sites of each site in order. It writes them as a reference plan that TKApp
# can read. The TKDash window (Dashboard.jl) shows all of this

"""
    SurveyRun

One continuous recording that [`scan_survey`](@ref) found. The scan reads only
the headers.

# Fields
- `path::String` -- the XML (Metronix) or the data file of the run.
- `sample_rate::Float64` -- samples per second.
- `start::DateTime` -- time of the first sample.
- `stop::DateTime` -- time immediately after the last sample.
- `n_samples::Int` -- samples for each channel.
- `components::Vector{Symbol}` -- the recorded channels (`:e1`, `:e2`, `:bx`,
  `:by`, `:bz`). For a Metronix run, these are all its `.ats` files. For
  LEMI-424 and GEOMAG, these are the columns that are not zero or `NaN` at the
  two ends of the file.
"""
struct SurveyRun
    path::String
    sample_rate::Float64
    start::DateTime
    stop::DateTime
    n_samples::Int
    components::Vector{Symbol}
end

"""
    SurveySite

A site that [`scan_survey`](@ref) found: its runs and its position.

# Fields
- `name::String` -- the name of the site directory.
- `path::String` -- the site directory.
- `format::Symbol` -- `:metronix`, `:lemi424` or `:geomag`.
- `latitude::Float64`, `longitude::Float64` -- decimal degrees, `NaN` if the
  headers give no position.
- `elevation::Float64` -- metres, `NaN` if not known.
- `runs::Vector{SurveyRun}` -- in order of start time.
"""
struct SurveySite
    name::String
    path::String
    format::Symbol
    latitude::Float64
    longitude::Float64
    elevation::Float64
    runs::Vector{SurveyRun}
end

"""
    Survey

All the sites in a survey directory, as [`scan_survey`](@ref) returns them. Use
a position or a site name as the index: `survey["site004"]`.
"""
struct Survey
    root::String
    sites::Vector{SurveySite}
end

Base.length(s::Survey) = length(s.sites)
Base.iterate(s::Survey, i = 1) = i > length(s.sites) ? nothing : (s.sites[i], i + 1)
Base.getindex(s::Survey, i::Integer) = s.sites[i]
function Base.getindex(s::Survey, name::AbstractString)
    i = findfirst(x -> x.name == name, s.sites)
    i === nothing && throw(KeyError(name))
    return s.sites[i]
end

"""
    site_components(site) -> Vector{Symbol}

All the channels that a run of a [`SurveySite`](@ref) recorded, with the
electric channels first.
"""
function site_components(site::SurveySite)
    all = unique(c for r in site.runs for c in r.components)
    order = Dict(c => i for (i, c) in enumerate((:e1, :e2, :e3, :e4, :bx, :by, :bz)))
    return sort!(all; by = c -> (get(order, c, 99), string(c)))
end

"""
    has_magnetic(site) -> Bool

Tells if a [`SurveySite`](@ref) recorded the two horizontal magnetic channels,
Hx and Hy. A base or remote site must supply these channels. A site with only
electric channels (a telluric site) can have base and remote sites, but it
can never be one.
"""
has_magnetic(site::SurveySite) = (c = site_components(site); :bx in c && :by in c)

# "Ex Ey Hx Hy Hz", the names that the instruments use
_channel_names(comps) = join((get(Dict(:e1 => "Ex", :e2 => "Ey", :bx => "Hx", :by => "Hy", :bz => "Hz"),
                                  c, string(c)) for c in comps), " ")

function Base.show(io::IO, site::SurveySite)
    rates = join(_fs_label.(site_rates(site)), ", ")
    print(io, "SurveySite(\"$(site.name)\", :$(site.format), $(length(site.runs)) runs at $rates)")
end

function Base.show(io::IO, ::MIME"text/plain", s::Survey)
    print(io, "Survey of $(length(s.sites)) sites in $(s.root)")
    for site in s.sites
        pos = isfinite(site.latitude) ? @sprintf("%.4f°, %.4f°", site.latitude, site.longitude) : "no position"
        @printf(io, "\n  %-16s %-9s %3d runs  %6.1f h  %-22s %s", site.name, site.format,
                length(site.runs), recording_seconds(site) / 3600, pos,
                _channel_names(site_components(site)))
    end
end

_fs_label(fs::Real) = fs >= 1 ? (isinteger(fs) ? "$(Int(fs)) Hz" : "$(fs) Hz") :
                        "$(round(1 / fs; sigdigits = 4)) s"

#---------- scanning -----

"""
    scan_survey(root; include_split = false, maxdepth = 4, progress = nothing) -> Survey

Find all the sites in `root` and make an index of their runs. The function
reads only the headers. It does not read the samples.

- A site is a directory that holds Metronix `meas_*` directories (or `.ats`
  files), or LEMI-424 or GEOMAG files.
- The function looks into the directories that are not sites, down to
  `maxdepth` levels.
- The position of each site comes from its headers: the `.ats` header of a
  Metronix run, the GPS columns of a LEMI-424 record, or the header block of
  a GEOMAG file.

The function ignores the rate directories that a split of a Metronix site
makes (`site002.128` next to `site002`). Thus, it does not count a site two
times. To make an index of them as separate sites, give `include_split = true`.

The function uses the first and the last line of a LEMI-424 or GEOMAG file as
the start and the end of the recording. Thus, it does not see gaps in one
file.

`progress`, when given, is a function that receives a line of text for each
directory that the function looks into and for each site that it finds.

See also [`site_references`](@ref), [`reference_plan`](@ref), [`run_tkdash`](@ref).
"""
function scan_survey(root::AbstractString; include_split::Bool = false, maxdepth::Integer = 4,
                     progress = nothing)
    root = _norm_path(root)
    isdir(root) || error("Not a directory: $root")
    sites = SurveySite[]
    say = progress === nothing ? (_ -> nothing) : progress
    _scan_dir!(sites, root, 0, maxdepth, include_split, say, root)
    sort!(sites; by = s -> s.name)
    # two sites with the same name in different subdirectories keep their
    # relative paths
    names = [s.name for s in sites]
    for (i, s) in enumerate(sites)
        count(==(s.name), names) > 1 || continue
        sites[i] = SurveySite(relpath(s.path, root), s.path, s.format, s.latitude,
                              s.longitude, s.elevation, s.runs)
    end
    return Survey(root, sites)
end

function _scan_dir!(sites, dir, depth, maxdepth, include_split, say = _ -> nothing, root = dir)
    say("Looking in $(dir == root ? basename(root) : relpath(dir, root))")
    site = try
        _scan_site(dir)
    catch err
        @warn "Skipping $(basename(dir)): its headers could not be read" exception = err
        nothing
    end
    if site !== nothing
        push!(sites, site)
        say("Found $(site.name) ($(site.format), $(length(site.runs)) run$(length(site.runs) == 1 ? "" : "s")) · " *
            "$(length(sites)) site$(length(sites) == 1 ? "" : "s") so far")
        return sites
    end
    depth >= maxdepth && return sites
    names = readdir(dir; sort = true)
    for name in names
        startswith(name, ".") && continue
        full = joinpath(dir, name)
        isdir(full) || continue
        !include_split && _is_split_dir(name, names) && continue
        _scan_dir!(sites, full, depth + 1, maxdepth, include_split, say, root)
    end
    return sites
end

# "site002.128" next to "site002": a rate directory from the split of a site
function _is_split_dir(name::AbstractString, siblings)
    m = match(r"^(.+)\.\d+(?:\.\d+)?$", name)
    return m !== nothing && m.captures[1] in siblings
end

# The site that `dir` is, or nothing if `dir` holds no recordings itself
function _scan_site(dir::AbstractString)
    is_metronix_site(dir) && return _scan_metronix_site(dir)
    files = _survey_data_files(dir)
    isempty(files) && return nothing
    runs = SurveyRun[]
    pos = (NaN, NaN, NaN)
    fmt = :lemi424
    for f in files
        scanned = try
            _scan_text_run(f)
        catch err
            @warn "Skipping $(basename(f)): not a LEMI-424 or GEOMAG record" exception = err
            nothing
        end
        scanned === nothing && continue
        run, p, fmt = scanned
        push!(runs, run)
        isfinite(pos[1]) || (pos = p)
    end
    isempty(runs) && return nothing
    sort!(runs; by = r -> r.start)
    return SurveySite(_site_name_from_dir(dir), _norm_path(dir), fmt, pos..., runs)
end

# The LEMI-424 and GEOMAG files of a site directory. A clean copy
# (`*_clean.txt`) or a joined `<site>.txt` next to raw files contains the same
# data again. Thus, the function ignores it. If such a file is alone, it is the
# data of the site
function _survey_data_files(dir::AbstractString)
    site = lowercase(_site_name_from_dir(dir))
    files = [joinpath(dir, n) for n in readdir(dir; sort = true)
             if isfile(joinpath(dir, n)) && lowercase(splitext(n)[2]) in _SITE_DATA_EXTS]
    copy(f) = occursin(r"_clean\.[^.]+$"i, basename(f)) || lowercase(splitext(basename(f))[1]) == site
    raw = filter(!copy, files)
    return isempty(raw) ? files : raw
end

const _ATS_OFF_LATITUDE = 96        # Int32, milliseconds of arc
const _ATS_OFF_LONGITUDE = 100      # Int32, milliseconds of arc
const _ATS_OFF_ELEVATION = 104      # Int32, centimetres

# (latitude, longitude, elevation) from an .ats header. NaN if not set
function _ats_position(hbytes::Vector{UInt8})
    length(hbytes) >= _ATS_OFF_ELEVATION + 4 || return (NaN, NaN, NaN)
    lat = _ats_get(Int32, hbytes, _ATS_OFF_LATITUDE) / 3.6e6
    lon = _ats_get(Int32, hbytes, _ATS_OFF_LONGITUDE) / 3.6e6
    elev = _ats_get(Int32, hbytes, _ATS_OFF_ELEVATION) / 100
    (lat == 0 && lon == 0) && return (NaN, NaN, NaN)
    (abs(lat) <= 90 && abs(lon) <= 180) || return (NaN, NaN, NaN)
    return (Float64(lat), Float64(lon), Float64(elev))
end

function _scan_metronix_site(dir::AbstractString)
    index, _ = _metronix_site_index(dir)
    isempty(index) && return nothing
    runs = SurveyRun[]
    for r in index
        start = Dates.unix2datetime(r.start_unix)
        stop = start + Millisecond(round(Int, 1000 * r.n_samples / r.rate))
        comps = unique(Symbol[let ct = _read_ats_header(a)["channel_type"]
                                  get(METRONIX_CHANNEL_MAP, ct, _symbolize(ct))
                              end for a in r.ats])
        push!(runs, SurveyRun(_run_id(r), r.rate, start, stop, r.n_samples, comps))
    end
    pos = (NaN, NaN, NaN)
    for r in index, ats in r.ats
        pos = _ats_position(_read_ats_header(ats)["header_bytes"])
        isfinite(pos[1]) && break
    end
    # the name of one meas_ directory is the name of the site above it
    name = _has_ats(dir) ? basename(dirname(_norm_path(dir))) : basename(_norm_path(dir))
    return SurveySite(name, _norm_path(dir), :metronix, pos..., runs)
end

_text_data_line(line) = !(all(isspace, line) || startswith(strip(line), ";"))

# The first `k` data lines of a LEMI-424 or GEOMAG file
function _first_data_lines(path::AbstractString, k::Integer)
    lines = String[]
    open(path, "r") do io
        for line in eachline(io)
            _text_data_line(line) || continue
            push!(lines, line)
            length(lines) == k && break
        end
    end
    return lines
end

# The last data line. The function reads the end of the file. Thus, it never
# reads the body
function _last_data_line(path::AbstractString)
    sz = filesize(path)
    return open(path, "r") do io
        chunk = min(sz, 8192)
        while true
            seek(io, sz - chunk)
            lines = split(read(io, String), '\n')
            chunk < sz && popfirst!(lines)               # the first line is not complete
            i = findlast(_text_data_line, lines)
            i === nothing || return String(lines[i])
            chunk == sz && return nothing
            chunk = min(sz, 4chunk)
        end
    end
end

function _text_line_time(line::AbstractString, fmt::Symbol)
    t = split(line)
    fmt == :geomag &&
        return _geomag_timestamp(parse.(Int, t[1:5])..., parse(Float64, t[6]))
    return DateTime(parse.(Int, t[1:6])...)
end

# The channels of a LEMI-424 or GEOMAG record that hold a value (finite and not
# zero) on one or more of `lines`, the first and the last lines of the file.
# The logger writes zeros for an input that is not connected
function _text_components(lines, fmt::Symbol)
    cols = fmt == :geomag ? GEOMAG_COLUMN_INDEX : LEMI424_DEFAULT_COLUMN_INDEX
    comps = Symbol[]
    for c in (:e1, :e2, :bx, :by, :bz)
        k = cols[c]
        live = any(lines) do line
            t = split(line)
            v = length(t) >= k ? tryparse(Float64, t[k]) : nothing
            v !== nothing && isfinite(v) && v != 0
        end
        live && push!(comps, c)
    end
    return comps
end

# A data line starts with the year, month, day, hour, minute and second
const _RECORD_START = r"^\s*\d{4}\s+\d{1,2}\s+\d{1,2}\s+\d{1,2}\s+\d{1,2}\s+\d"

# (run, (latitude, longitude, elevation), format) of one LEMI-424 or GEOMAG
# file, or nothing if the file holds no records with timestamps
function _scan_text_run(path::AbstractString)
    fmt = endswith(lowercase(path), ".txt") ? _detect_format(path) : :lemi424
    head = _first_data_lines(path, 2)
    # calibration tables and other text without timestamped records
    (isempty(head) || !occursin(_RECORD_START, head[1])) && return nothing
    tail = _last_data_line(path)
    t0 = _text_line_time(head[1], fmt)
    t1 = _text_line_time(tail, fmt)
    pos = (NaN, NaN, NaN)
    fs = NaN
    if fmt == :geomag
        meta = _read_geomag_metadata(path)
        fs = meta[:sample_rate]
        pos = (meta[:latitude], meta[:longitude], meta[:elevation])
    else
        tok = split(head[1])
        if length(tok) >= 21
            lat = tryparse(Float64, tok[18])
            lon = tryparse(Float64, tok[20])
            elev = tryparse(Float64, tok[17])
            if lat !== nothing && lon !== nothing && !(lat == 0 && lon == 0)
                pos = (_lemi_position(lat) * _lemi_hemisphere_sign(tok[19]),
                       _lemi_position(lon) * _lemi_hemisphere_sign(tok[21]),
                       something(elev, NaN))
            end
        end
    end
    if !(isfinite(fs) && fs > 0)
        dt = length(head) == 2 ? Dates.value(_text_line_time(head[2], fmt) - t0) / 1000 : 1.0
        fs = dt > 0 ? 1 / dt : 1.0
    end
    n = round(Int, Dates.value(t1 - t0) / 1000 * fs) + 1
    stop = t1 + Millisecond(round(Int, 1000 / fs))
    comps = _text_components(vcat(head, tail), fmt)
    return SurveyRun(_norm_path(path), fs, t0, stop, n, comps), pos, fmt
end

#---------- time overlap -----

const _Interval = Tuple{Float64, Float64}     # unix seconds, [start, stop)

_unix(t::DateTime) = Dates.datetime2unix(t)
_datetime(s::Real) = Dates.unix2datetime(s)

_same_rate(a::Real, b::Real) = isapprox(a, b; rtol = 1.0e-6)

"""
    site_rates(site) -> Vector{Float64}

The different sampling rates of a [`SurveySite`](@ref), in increasing order.
"""
site_rates(site::SurveySite) = sort!(unique(_rate_key.(r.sample_rate for r in site.runs)))

"""
    survey_rates(survey) -> Vector{Float64}

The different sampling rates of all the sites of a [`Survey`](@ref), in
increasing order.
"""
survey_rates(s::Survey) = sort!(unique(reduce(vcat, (site_rates(x) for x in s.sites); init = Float64[])))

function _merge_intervals!(iv::Vector{_Interval})
    isempty(iv) && return iv
    sort!(iv)
    out = _Interval[iv[1]]
    for (a, b) in @view iv[2:end]
        if a <= out[end][2]
            out[end] = (out[end][1], max(out[end][2], b))
        else
            push!(out, (a, b))
        end
    end
    return out
end

function _intersect_intervals(a::Vector{_Interval}, b::Vector{_Interval})
    out = _Interval[]
    i = j = 1
    while i <= length(a) && j <= length(b)
        lo = max(a[i][1], b[j][1])
        hi = min(a[i][2], b[j][2])
        hi > lo && push!(out, (lo, hi))
        a[i][2] < b[j][2] ? (i += 1) : (j += 1)
    end
    return out
end

_total(iv) = sum((b - a for (a, b) in iv); init = 0.0)

# Tells if a run counts at `rate`. A number selects that rate. `nothing` (a
# shared rate) and `:all` (all rates) select all the runs
_at_rate(r::SurveyRun, rate) = rate === nothing || rate === :all || _same_rate(r.sample_rate, rate)

function _site_intervals(site::SurveySite, rate)
    iv = [(_unix(r.start), _unix(r.stop)) for r in site.runs if _at_rate(r, rate)]
    return _merge_intervals!(iv)
end

# Where `a` and `b` both record: at one rate. With `rate = nothing`, at a rate
# that they share. With `rate = :all`, at all rates
function _overlap(a::SurveySite, b::SurveySite, rate)
    rate === :all && return _intersect_intervals(_site_intervals(a, :all), _site_intervals(b, :all))
    rates = rate === nothing ? intersect(site_rates(a), site_rates(b)) : [rate]
    acc = _Interval[]
    for fs in rates
        append!(acc, _intersect_intervals(_site_intervals(a, fs), _site_intervals(b, fs)))
    end
    return _merge_intervals!(acc)
end

"""
    recording_intervals(site; rate = nothing) -> Vector{Tuple{DateTime, DateTime}}

The spans that `site` recorded, at `rate`, or at all rates if `rate` is
`nothing` or `:all`. Runs that touch become one span.
"""
recording_intervals(site::SurveySite; rate = nothing) =
    [(_datetime(a), _datetime(b)) for (a, b) in _site_intervals(site, rate)]

"""
    recording_seconds(site; rate = nothing) -> Float64

The total time that `site` recorded, at `rate`, or at all rates if `rate` is
`nothing` or `:all`.
"""
recording_seconds(site::SurveySite; rate = nothing) = _total(_site_intervals(site, rate))

"""
    overlap_intervals(a, b; rate = nothing) -> Vector{Tuple{DateTime, DateTime}}

The spans in which the sites `a` and `b` both recorded. `rate` sets which
recordings make a pair:

- a number -- the two sites at that sampling rate;
- `nothing` (the default) -- the two sites at the same rate, at a rate that
  they share;
- `:all` -- at all rates. Use it for a survey with different instruments (for
  example, a LEMI-424 at 1 Hz and a Metronix at 128 Hz). The processing will
  resample these records.
"""
overlap_intervals(a::SurveySite, b::SurveySite; rate = nothing) =
    [(_datetime(x), _datetime(y)) for (x, y) in _overlap(a, b, rate)]

"""
    overlap_seconds(a, b; rate = nothing) -> Float64

The total length of [`overlap_intervals`](@ref)`(a, b; rate)`.
"""
overlap_seconds(a::SurveySite, b::SurveySite; rate = nothing) = _total(_overlap(a, b, rate))

"""
    overlap_matrix(survey; rate = nothing) -> Matrix{Float64}

The hours that each pair of sites recorded together, in the order of
`survey.sites`. The diagonal holds the recording hours of each site.
"""
function overlap_matrix(s::Survey; rate = nothing)
    n = length(s.sites)
    m = zeros(n, n)
    for i in 1:n
        m[i, i] = recording_seconds(s.sites[i]; rate) / 3600
        for j in (i + 1):n
            m[i, j] = m[j, i] = overlap_seconds(s.sites[i], s.sites[j]; rate) / 3600
        end
    end
    return m
end

#---------- distance -----

const _EARTH_RADIUS_KM = 6371.0088

"""
    site_distance(a, b) -> Float64

The great-circle distance between two sites, in kilometres. `NaN` if one of
the sites has no position.
"""
function site_distance(a::SurveySite, b::SurveySite)
    φ1, φ2 = deg2rad(a.latitude), deg2rad(b.latitude)
    Δφ = φ2 - φ1
    Δλ = deg2rad(b.longitude - a.longitude)
    h = sin(Δφ / 2)^2 + cos(φ1) * cos(φ2) * sin(Δλ / 2)^2
    return 2 * _EARTH_RADIUS_KM * asin(min(1.0, sqrt(h)))
end

#---------- base and remote sites -----

"""
    site_references(survey, site; rate = nothing, base_km = 5.0, remote_km = 20.0,
                    min_overlap_hours = 1.0, exclude = ()) -> (base, remote)

The sites that recorded together with `site` (a [`SurveySite`](@ref), its
name or its index) for `min_overlap_hours` or more, at a common rate (refer
to [`overlap_intervals`](@ref)):

- **base** sites: within `base_km` of the site;
- **remote** sites: `remote_km` or more from the site.

`rate` has the same meaning as for [`overlap_intervals`](@ref). Base and
remote sites supply their magnetic field. Thus, only sites that recorded Hx
and Hy count (refer to [`has_magnetic`](@ref)). These sites are in neither
list: sites between the two distances, sites without a position in their
headers, and sites with only electric channels. With `remote_km = base_km`,
there is no gap.

Each list is in order of the time recorded together, longest first. If two
sites have the same time, the nearer site comes first. An entry is a
`NamedTuple`:

- `site` -- the name of the site
- `distance_km` -- the great-circle distance from `site`
- `overlap_hours` -- the time recorded together
- `overlap_fraction` -- that time as a fraction of the recording of `site`
- `windows` -- the overlap spans, as `(start, stop)` pairs
- `excluded` -- tells if its name is in `exclude`. Excluded sites stay in the
  lists with this flag. Thus, a window can show them.
  [`reference_plan`](@ref) does not include them
"""
function site_references(s::Survey, site; rate = nothing, base_km::Real = 5.0,
                         remote_km::Real = 20.0, min_overlap_hours::Real = 1.0, exclude = ())
    t = _target_site(s, site)
    own = recording_seconds(t; rate)
    base, remote = NamedTuple[], NamedTuple[]
    for c in s.sites
        (c.path == t.path || !has_magnetic(c)) && continue
        iv = _overlap(t, c, rate)
        h = _total(iv) / 3600
        (h > 0 && h >= min_overlap_hours) || continue
        d = site_distance(t, c)
        isfinite(d) || continue
        list = d <= base_km ? base : d >= remote_km ? remote : nothing
        list === nothing && continue
        push!(list,
              (site = c.name, distance_km = d, overlap_hours = h,
               overlap_fraction = own > 0 ? h * 3600 / own : 0.0,
               windows = [(_datetime(a), _datetime(b)) for (a, b) in iv],
               excluded = c.name in exclude))
    end
    order(c) = (-round(c.overlap_hours; digits = 6), c.distance_km)
    return (base = sort!(base; by = order), remote = sort!(remote; by = order))
end

_target_site(s::Survey, t::SurveySite) = t
_target_site(s::Survey, t::AbstractString) = s[t]
_target_site(s::Survey, t::Integer) = s.sites[t]

"""
    common_window(survey, site, others; rate = nothing) -> Vector{Tuple{DateTime, DateTime}}

The spans in which `site` and all the sites in `others` (sites or names)
recorded together: at `rate`, at one rate that they all share (`nothing`), or
at all rates (`:all`). Processing can use this part of the data when it uses
`site` with all of `others` at the same time. The result is empty if `others`
is empty or if the sites never recorded together.
"""
function common_window(s::Survey, site, others; rate = nothing)
    t = _target_site(s, site)
    os = [_target_site(s, o) for o in others]
    isempty(os) && return Tuple{DateTime, DateTime}[]
    rates = rate === :all ? [:all] : rate === nothing ? site_rates(t) : [rate]
    acc = _Interval[]
    for fs in rates
        iv = _site_intervals(t, fs)
        for o in os
            isempty(iv) && break
            iv = _intersect_intervals(iv, _site_intervals(o, fs))
        end
        append!(acc, iv)
    end
    return [(_datetime(a), _datetime(b)) for (a, b) in _merge_intervals!(acc)]
end

_window_hours(w) = sum((Dates.value(b - a) / 3.6e6 for (a, b) in w); init = 0.0)

"""
    reference_plan(survey; rate = nothing, base_km = 5.0, remote_km = 20.0,
                   min_overlap_hours = 1.0, exclude = Dict()) -> Vector{NamedTuple}

All the sites of `survey` with their base and remote sites from
[`site_references`](@ref). The result has one `NamedTuple` for each site:
`(site, base, base_hours, remote, remote_hours, common)`.

- `base` and `remote` are vectors of site names, with the longest overlap
  first.
- `base_hours` and `remote_hours` are the hours that each of these sites
  recorded with the site, in the same order.
- `common` is the [`common_window`](@ref) of the site with all of them.

`exclude[site]` is a set of names that the lists of that site must not
include.

See also [`write_reference_plan`](@ref), [`read_reference_plan`](@ref).
"""
function reference_plan(s::Survey; rate = nothing, base_km::Real = 5.0, remote_km::Real = 20.0,
                        min_overlap_hours::Real = 1.0, exclude = Dict{String, Set{String}}())
    return map(s.sites) do site
        refs = site_references(s, site; rate, base_km, remote_km, min_overlap_hours,
                               exclude = get(exclude, site.name, ()))
        kept(list) = [c for c in list if !c.excluded]
        base, remote = kept(refs.base), kept(refs.remote)
        names = String[c.site for c in vcat(base, remote)]
        (site = site.name,
         base = String[c.site for c in base], base_hours = Float64[c.overlap_hours for c in base],
         remote = String[c.site for c in remote], remote_hours = Float64[c.overlap_hours for c in remote],
         common = common_window(s, site, names; rate))
    end
end

const _PLAN_COLUMNS = ("site", "base", "overlap (h)", "remote", "overlap (h)")

"""
    write_reference_plan(path, survey; kwargs...) -> String
    write_reference_plan(path, plan) -> String

Write a [`reference_plan`](@ref) as a plain text table. The table has one
header line and one row for each site. The columns are aligned:

```
site      base                       overlap (h)          remote              overlap (h)
site002   site004, site006           11.77, 9.00          site099, site100    11.77, 11.77
site099   -                          -                    site100, site004    22.98, 12.00
```

The base and remote sites are in order of longest overlap first. The overlaps
below `overlap (h)` are in the same order. A `-` shows an empty list. The
keywords are the same as for [`reference_plan`](@ref). Use
[`read_reference_plan`](@ref) to read the table. The function returns `path`.
"""
write_reference_plan(path::AbstractString, s::Survey; kwargs...) =
    write_reference_plan(path, reference_plan(s; kwargs...))

function write_reference_plan(path::AbstractString, plan::AbstractVector)
    list(v) = isempty(v) ? "-" : join(v, ", ")
    hours(v) = isempty(v) ? "-" : join((@sprintf("%.2f", h) for h in v), ", ")
    rows = [[r.site, list(r.base), hours(r.base_hours), list(r.remote), hours(r.remote_hours)]
            for r in plan]
    widths = [maximum(length, (row[k] for row in vcat([collect(_PLAN_COLUMNS)], rows))) for k in 1:5]
    line(cells) = rstrip(join((rpad(c, w) for (c, w) in zip(cells, widths)), "   "))
    mkpath(dirname(abspath(path)))
    open(path, "w") do io
        println(io, line(_PLAN_COLUMNS))
        foreach(row -> println(io, line(row)), rows)
    end
    return path
end

"""
    read_reference_plan(path) -> Vector{NamedTuple}

Read a table that [`write_reference_plan`](@ref) (or **Export** in TKDash)
wrote. The result has one `(site, base, base_hours, remote, remote_hours)` for
each row. The sites are vectors of names, and the hours are vectors of
numbers, in the order of the table.
"""
function read_reference_plan(path::AbstractString)
    lines = filter(!isempty ∘ strip, readlines(path))
    isempty(lines) && error("Empty reference plan: $path")
    head = lines[1]
    # each column starts at the position of its header
    starts = Int[]
    from = 1
    for name in _PLAN_COLUMNS
        r = findnext(name, head, from)
        r === nothing && error("$path is not a reference plan: no '$name' column")
        push!(starts, first(r))
        from = last(r) + 1
    end
    # cut by character, not by byte. Thus, names such as Sarıçam keep their
    # columns
    function cell(chars, k)
        stop = k < 5 ? min(starts[k + 1] - 1, length(chars)) : length(chars)
        return starts[k] > length(chars) ? "" : strip(String(chars[starts[k]:stop]))
    end
    list(x) = x == "-" || isempty(x) ? String[] : String[strip(n) for n in split(x, ',')]
    hours(x) = Float64[parse(Float64, h) for h in list(x)]
    return map(lines[2:end]) do line
        c = collect(line)
        (site = cell(c, 1), base = list(cell(c, 2)), base_hours = hours(cell(c, 3)),
         remote = list(cell(c, 4)), remote_hours = hours(cell(c, 5)))
    end
end
