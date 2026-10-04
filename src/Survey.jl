# Survey.jl - survey scanning and reference-site selection.
# Author: @pankajkmishra
#
# A survey is a directory of sites recorded over the same weeks. For each site
# the sites that recorded at the same time and rate are its base sites when
# close and its remote sites when far. Finding them needs only when,
# how fast and where each site recorded, so this file reads headers alone: the
# .ats headers of a Metronix site, and the first and last lines of LEMI-424
# and GEOMAG files.
#
# On top of that index it measures the time two sites recorded together at a
# common rate and the distance between them, sorts each site's base and remote
# sites, and writes them as a reference plan that TKApp can read back. The
# TKDash window (Dashboard.jl) draws all of it

"""
    SurveyRun

One continuous recording found by [`scan_survey`](@ref), read from headers
only.

# Fields
- `path::String` -- the run's XML (Metronix) or data file.
- `sample_rate::Float64` -- samples per second.
- `start::DateTime` -- time of the first sample.
- `stop::DateTime` -- time just after the last sample.
- `n_samples::Int` -- samples per channel.
- `components::Vector{Symbol}` -- the channels recorded (`:e1`, `:e2`, `:bx`,
  `:by`, `:bz`): every `.ats` file of a Metronix run; for LEMI-424 and GEOMAG
  the columns that are not zero or `NaN` at both ends of the file.
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

A site found by [`scan_survey`](@ref): its runs and where it stood.

# Fields
- `name::String` -- the site directory's name.
- `path::String` -- the site directory.
- `format::Symbol` -- `:metronix`, `:lemi424` or `:geomag`.
- `latitude::Float64`, `longitude::Float64` -- decimal degrees, `NaN` if the
  headers carry no position.
- `elevation::Float64` -- metres, `NaN` if unknown.
- `runs::Vector{SurveyRun}` -- ordered by start time.
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

Every site under a survey directory, as returned by [`scan_survey`](@ref).
Index it by position or by site name: `survey["site004"]`.
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

Every channel any run of a [`SurveySite`](@ref) recorded, electric first.
"""
function site_components(site::SurveySite)
    all = unique(c for r in site.runs for c in r.components)
    order = Dict(c => i for (i, c) in enumerate((:e1, :e2, :e3, :e4, :bx, :by, :bz)))
    return sort!(all; by = c -> (get(order, c, 99), string(c)))
end

"""
    has_magnetic(site) -> Bool

Whether a [`SurveySite`](@ref) recorded both horizontal magnetic channels,
Hx and Hy - what a base or remote site has to offer. A site with electric
channels alone (a telluric site) can be paired with base and remote sites but
never be one.
"""
has_magnetic(site::SurveySite) = (c = site_components(site); :bx in c && :by in c)

# "Ex Ey Hx Hy Hz", as the instruments name them
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
    scan_survey(root; include_split = false, maxdepth = 4) -> Survey

Find every site under `root` and index its runs from headers alone, without
reading samples. A site is a directory holding Metronix `meas_*` directories
(or `.ats` files), or LEMI-424 / GEOMAG files; directories that are not sites
are searched down to `maxdepth` levels. Each site's position comes from its
headers: the `.ats` header of a Metronix run, the GPS columns of a LEMI-424
record, the header block of a GEOMAG file.

The rate directories that splitting a Metronix site makes (`site002.128`
beside `site002`) are skipped, so a site is not counted twice; pass
`include_split = true` to index them as sites of their own.

A LEMI-424 or GEOMAG file is taken as recording from its first line to its
last; gaps inside one file are not seen.

See also [`site_references`](@ref), [`reference_plan`](@ref), [`run_tkdash`](@ref).
"""
function scan_survey(root::AbstractString; include_split::Bool = false, maxdepth::Integer = 4)
    root = _norm_path(root)
    isdir(root) || error("Not a directory: $root")
    sites = SurveySite[]
    _scan_dir!(sites, root, 0, maxdepth, include_split)
    sort!(sites; by = s -> s.name)
    # two sites of one name in different subdirectories keep their relative path
    names = [s.name for s in sites]
    for (i, s) in enumerate(sites)
        count(==(s.name), names) > 1 || continue
        sites[i] = SurveySite(relpath(s.path, root), s.path, s.format, s.latitude,
                              s.longitude, s.elevation, s.runs)
    end
    return Survey(root, sites)
end

function _scan_dir!(sites, dir, depth, maxdepth, include_split)
    site = try
        _scan_site(dir)
    catch err
        @warn "Skipping $(basename(dir)): its headers could not be read" exception = err
        nothing
    end
    if site !== nothing
        push!(sites, site)
        return sites
    end
    depth >= maxdepth && return sites
    names = readdir(dir; sort = true)
    for name in names
        startswith(name, ".") && continue
        full = joinpath(dir, name)
        isdir(full) || continue
        !include_split && _is_split_dir(name, names) && continue
        _scan_dir!(sites, full, depth + 1, maxdepth, include_split)
    end
    return sites
end

# "site002.128" beside "site002": a rate directory split off a site
function _is_split_dir(name::AbstractString, siblings)
    m = match(r"^(.+)\.\d+(?:\.\d+)?$", name)
    return m !== nothing && m.captures[1] in siblings
end

# The site `dir` is, or nothing when it holds no recordings of its own
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

# The LEMI-424 and GEOMAG files of a site directory. A cleaned copy
# (`*_clean.txt`) or a combined `<site>.txt` written beside raw files is the
# same data again and left out; alone, either is the site's data
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

# (latitude, longitude, elevation) from an .ats header; NaN when unset
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
    # a single meas_ directory is named after the site above it
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

# The last data line, read from the file's tail so the body is never scanned
function _last_data_line(path::AbstractString)
    sz = filesize(path)
    return open(path, "r") do io
        chunk = min(sz, 8192)
        while true
            seek(io, sz - chunk)
            lines = split(read(io, String), '\n')
            chunk < sz && popfirst!(lines)               # cut mid-line
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

# The channels of a LEMI-424 or GEOMAG record that hold a value - finite and
# not zero - on any of `lines`, the first and last of the file; an unconnected
# input is logged as zeros
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

# A data line opens with year, month, day, hour, minute, second
const _RECORD_START = r"^\s*\d{4}\s+\d{1,2}\s+\d{1,2}\s+\d{1,2}\s+\d{1,2}\s+\d"

# (run, (latitude, longitude, elevation), format) of one LEMI-424 or GEOMAG
# file, or nothing when the file holds no timestamped records
function _scan_text_run(path::AbstractString)
    fmt = endswith(lowercase(path), ".txt") ? _detect_format(path) : :lemi424
    head = _first_data_lines(path, 2)
    # calibration tables and other text that is not a timestamped record
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

Sorted unique sampling rates of a [`SurveySite`](@ref).
"""
site_rates(site::SurveySite) = sort!(unique(_rate_key.(r.sample_rate for r in site.runs)))

"""
    survey_rates(survey) -> Vector{Float64}

Sorted unique sampling rates across every site of a [`Survey`](@ref).
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

# Whether a run counts at `rate`: a number picks that rate; `nothing` (any
# rate shared) and `:all` (every rate) take every run
_at_rate(r::SurveyRun, rate) = rate === nothing || rate === :all || _same_rate(r.sample_rate, rate)

function _site_intervals(site::SurveySite, rate)
    iv = [(_unix(r.start), _unix(r.stop)) for r in site.runs if _at_rate(r, rate)]
    return _merge_intervals!(iv)
end

# Where `a` and `b` both record: at one rate; with `rate = nothing`, at any
# rate they share; with `rate = :all`, at any rates at all
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

The spans `site` recorded, runs that touch merged into one, at `rate`, or at
any rate when `rate` is `nothing` or `:all`.
"""
recording_intervals(site::SurveySite; rate = nothing) =
    [(_datetime(a), _datetime(b)) for (a, b) in _site_intervals(site, rate)]

"""
    recording_seconds(site; rate = nothing) -> Float64

Total time `site` recorded, at `rate`, or at any rate when `rate` is `nothing`
or `:all`.
"""
recording_seconds(site::SurveySite; rate = nothing) = _total(_site_intervals(site, rate))

"""
    overlap_intervals(a, b; rate = nothing) -> Vector{Tuple{DateTime, DateTime}}

The spans in which sites `a` and `b` were both recording. `rate` sets which
recordings pair:

- a number -- both at that sampling rate;
- `nothing` (the default) -- both at the same rate, any rate the two share;
- `:all` -- at any rates at all, for a survey of mixed instruments (a LEMI-424
  at 1 Hz with a Metronix at 128 Hz) whose records processing will resample.
"""
overlap_intervals(a::SurveySite, b::SurveySite; rate = nothing) =
    [(_datetime(x), _datetime(y)) for (x, y) in _overlap(a, b, rate)]

"""
    overlap_seconds(a, b; rate = nothing) -> Float64

Total length of [`overlap_intervals`](@ref)`(a, b; rate)`.
"""
overlap_seconds(a::SurveySite, b::SurveySite; rate = nothing) = _total(_overlap(a, b, rate))

"""
    overlap_matrix(survey; rate = nothing) -> Matrix{Float64}

Hours each pair of sites recorded together, in the order of `survey.sites`;
the diagonal holds each site's own recording hours.
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

Great-circle distance between two sites in kilometres, `NaN` when either has
no position.
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
name or its index) for at least `min_overlap_hours` at a common rate (see
[`overlap_intervals`](@ref)): **base** sites within `base_km` of it, **remote**
sites `remote_km` or more away; `rate` is as for [`overlap_intervals`](@ref).
Both lend their magnetic field, so only sites
that recorded Hx and Hy count (see [`has_magnetic`](@ref)). Sites in between,
sites without a position in their headers and electric-only sites are in
neither; `remote_km = base_km` leaves no gap.

Each list is ordered by the time recorded together, longest first, a tie
going to the nearer site. An entry is a `NamedTuple`:

- `site` -- the site's name
- `distance_km` -- great-circle distance from `site`
- `overlap_hours` -- time recorded together
- `overlap_fraction` -- that time as a share of `site`'s own recording
- `windows` -- the overlap spans, as `(start, stop)` pairs
- `excluded` -- whether its name is in `exclude`; excluded sites stay in the
  lists, flagged, so a window can show them, and [`reference_plan`](@ref)
  leaves them out
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

The spans in which `site` and every site in `others` (sites or names) were all
recording: at `rate`, at any one rate they all share (`nothing`), or at any
rates at all (`:all`). This
is the stretch of data processing can use when `site` is paired with all of
`others` at once. Empty when `others` is empty or the sites never recorded
together.
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

Every site of `survey` with its base and remote sites from
[`site_references`](@ref), one `NamedTuple` per site:
`(site, base, base_hours, remote, remote_hours, common)`. `base` and `remote`
are vectors of site names, longest overlap first; `base_hours` and
`remote_hours` the hours each recorded with the site, in the same order; and
`common` the [`common_window`](@ref) of the site with all of them.
`exclude[site]` is a collection of names to leave out of that site's lists.

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

Write a [`reference_plan`](@ref) as a plain-text table, one header line and
one row per site, columns aligned:

```
site      base                       overlap (h)          remote              overlap (h)
site002   site004, site006           11.77, 9.00          site099, site100    11.77, 11.77
site099   -                          -                    site100, site004    22.98, 12.00
```

Base and remote sites are listed longest overlap first, each overlap below
`overlap (h)` in the same order; `-` marks an empty list. Keywords are those
of [`reference_plan`](@ref). [`read_reference_plan`](@ref) reads the table
back. Returns `path`.
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

Read a table written by [`write_reference_plan`](@ref) (or by TKDash's
**Export**): one `(site, base, base_hours, remote, remote_hours)` per row,
the sites as vectors of names and the hours as vectors of numbers, in the
table's order.
"""
function read_reference_plan(path::AbstractString)
    lines = filter(!isempty ∘ strip, readlines(path))
    isempty(lines) && error("Empty reference plan: $path")
    head = lines[1]
    # each column starts where its header does
    starts = Int[]
    from = 1
    for name in _PLAN_COLUMNS
        r = findnext(name, head, from)
        r === nothing && error("$path is not a reference plan: no '$name' column")
        push!(starts, first(r))
        from = last(r) + 1
    end
    # cut by character, not byte, so names such as Sarıçam keep their columns
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
