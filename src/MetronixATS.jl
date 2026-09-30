# MetronixATS.jl - Metronix ADU binary ATS format reader and writer.
# Author: @pankajkmishra
#
# A Metronix recording is a directory of per-channel .ats files (a binary
# header giving sample count, rate, start time and LSB scaling, followed by
# Int32 samples) plus an XML run descriptor. A site holds many meas_*
# directories, and one meas_* directory can hold several runs - each a run
# number (R000, R001, ...) at one sampling rate, with its own XML.
#
# Beyond reading and writing single runs, this file separates a site by
# sampling rate (DF002 -> DF002.TK/128, DF002.TK/4096, ...) and implements the
# amputation workflow: applying mask intervals to a whole site, splitting each
# run into its surviving contiguous segments, writing them as new meas_*
# directories under DF002.TK<date>_<time>/<rate> with headers and XML rewritten
# to match, and appending a README history of every write session.

"""
    METRONIX_CHANNEL_MAP

Maps Metronix ATS channel-type strings to Timekeepers component names:
`"Ex" => :e1`, `"Ey" => :e2`, `"Hx" => :bx`, `"Hy" => :by`, `"Hz" => :bz`.
Channel types outside this table are carried through as symbols of their own
name.
"""
const METRONIX_CHANNEL_MAP = Dict("Ex" => :e1, "Ey" => :e2, "Hx" => :bx, "Hy" => :by, "Hz" => :bz)
const METRONIX_DEFAULT_COMPONENTS = [:e1, :e2, :bx, :by, :bz]

const _ATS_OFF_SAMPLE_LENGTH = 4
const _ATS_OFF_SAMPLING_RATE = 8
const _ATS_OFF_START = 12
const _ATS_OFF_LSBVAL = 16
const _ATS_OFF_CHANNEL_TYPE = 38

_ats_get(::Type{T}, bytes::Vector{UInt8}, off::Int) where {T} =
    reinterpret(T, @view bytes[(off + 1):(off + sizeof(T))])[1]

function _ats_put!(bytes::Vector{UInt8}, off::Int, v::T) where {T}
    bytes[(off + 1):(off + sizeof(T))] = reinterpret(UInt8, [v])
    return bytes
end

function _parse_ats_header_bytes(hbytes::Vector{UInt8})
    ct_raw = hbytes[(_ATS_OFF_CHANNEL_TYPE + 1):(_ATS_OFF_CHANNEL_TYPE + 2)]
    return Dict{String, Any}(
        "header_bytes" => hbytes,
        "header_length" => length(hbytes),
        "sample_length" => Int(_ats_get(Int32, hbytes, _ATS_OFF_SAMPLE_LENGTH)),
        "sampling_rate" => Float64(_ats_get(Float32, hbytes, _ATS_OFF_SAMPLING_RATE)),
        "start_unix" => Int(_ats_get(Int32, hbytes, _ATS_OFF_START)),
        "lsbval" => _ats_get(Float64, hbytes, _ATS_OFF_LSBVAL),
        "channel_type" => String(filter(!=(0x00), ct_raw)),
    )
end

function _read_ats_header(path::AbstractString)
    open(path, "r") do f
        header_length = read(f, UInt16)
        seekstart(f)
        hbytes = read(f, Int(header_length))
        length(hbytes) == header_length || error("Truncated ATS header in $path")
        info = _parse_ats_header_bytes(hbytes)
        info["source_file"] = abspath(path)
        return info
    end
end

function _read_ats(path::AbstractString)
    open(path, "r") do f
        header_length = read(f, UInt16)
        seekstart(f)
        hbytes = read(f, Int(header_length))
        length(hbytes) == header_length || error("Truncated ATS header in $path")
        info = _parse_ats_header_bytes(hbytes)
        info["source_file"] = abspath(path)
        seek(f, Int(header_length))
        raw = Vector{Int32}(undef, info["sample_length"])
        read!(f, raw)
        data = Float64.(raw) .* info["lsbval"]
        return data, info
    end
end

function _write_ats(path::AbstractString, data::AbstractVector{<:Real}, header_bytes::Vector{UInt8},
                    lsbval::Real, start_unix::Integer)
    bytes = copy(header_bytes)
    _ats_put!(bytes, _ATS_OFF_SAMPLE_LENGTH, Int32(length(data)))
    _ats_put!(bytes, _ATS_OFF_START, Int32(start_unix))
    raw = round.(Int32, data ./ lsbval)
    open(path, "w") do f
        write(f, bytes)
        write(f, raw)
    end
    return path
end

"""
    _MetronixRunFiles

One recorded run inside a `meas_*` directory: the `.ats` files sharing a run
number (`R000`, `R001`, ...) and a sampling rate, plus the XML that describes
them. A `meas_*` directory can hold several of these - an ADU scheduling a
long 128 Hz run and 4096 Hz bursts writes them all into one directory.
"""
struct _MetronixRunFiles
    meas_dir::String
    run_token::String
    rate::Float64
    xml::Union{Nothing, String}
    ats::Vector{String}
    start_unix::Int
    n_samples::Int
end

# A run is named by its XML, the file a user picks to load it; the first .ats
# stands in when the XML is missing.
_run_id(r::_MetronixRunFiles) = something(r.xml, first(r.ats))

_rate_key(fs::Real) = round(Float64(fs); digits = 6)

_norm_path(p::AbstractString) = rstrip(abspath(String(p)), ['/', '\\'])

function _ats_run_token(path::AbstractString)
    for part in split(splitext(basename(path))[1], '_')
        occursin(r"^R\d+$", part) && return String(part)
    end
    return "R000"
end

# "4096H" is 4096 Hz; "8S" is one sample every 8 s.
function _freq_token_rate(token::AbstractString)
    m = match(r"^(\d+(?:\.\d+)?)([HhSs])$", token)
    m === nothing && return nothing
    v = parse(Float64, m.captures[1])
    return lowercase(m.captures[2]) == "h" ? v : 1 / v
end

function _xml_rate(xml_path::AbstractString)
    _, _, freq_token = _parse_xml_filename_tokens(xml_path)
    rate = _freq_token_rate(freq_token)
    rate === nothing || return rate
    try
        node = findfirst("//sample_freq", EzXML.readxml(xml_path))
        node === nothing || return parse(Float64, strip(node.content))
    catch
    end
    return nothing
end

"""
    _metronix_meas_runs(meas_dir) -> (runs, empty_xmls)

Group the `.ats` files of one `meas_*` directory into runs by run number and
the sampling rate in their headers, and pair each run with the XML whose
filename carries the same run number and rate. XMLs left without data - a job
the ADU scheduled but never recorded - come back in `empty_xmls`. Runs are
ordered by start time.
"""
function _metronix_meas_runs(meas_dir::AbstractString)
    meas_dir = _norm_path(meas_dir)
    isdir(meas_dir) || error("Not a directory: $meas_dir")
    groups = Dict{Tuple{String, Float64}, Vector{Tuple{String, Dict{String, Any}}}}()
    xmls = String[]
    for name in readdir(meas_dir; sort = true)
        full = joinpath(meas_dir, name)
        isfile(full) || continue
        ext = lowercase(splitext(name)[2])
        if ext == ".ats"
            info = _read_ats_header(full)
            key = (_ats_run_token(full), _rate_key(info["sampling_rate"]))
            push!(get!(groups, key, Tuple{String, Dict{String, Any}}[]), (full, info))
        elseif ext == ".xml"
            push!(xmls, full)
        end
    end

    xml_by_key = Dict{Tuple{String, Float64}, String}()
    loose_xmls = String[]
    for x in xmls
        _, run_token, _ = _parse_xml_filename_tokens(x)
        rate = _xml_rate(x)
        key = rate === nothing ? nothing : (run_token, _rate_key(rate))
        if key === nothing || haskey(xml_by_key, key)
            push!(loose_xmls, x)
        else
            xml_by_key[key] = x
        end
    end

    runs = _MetronixRunFiles[]
    unpaired = Tuple{String, Float64}[]
    for (key, files) in groups
        xml = pop!(xml_by_key, key, nothing)
        xml === nothing && push!(unpaired, key)
        info = files[1][2]
        push!(runs, _MetronixRunFiles(meas_dir, key[1], info["sampling_rate"], xml,
                                      first.(files), info["start_unix"],
                                      minimum(f[2]["sample_length"] for f in files)))
    end
    empty_xmls = vcat(collect(values(xml_by_key)), loose_xmls)
    # Files that break the naming convention: one run and one XML left over
    # belong together, as they did before runs were told apart.
    if length(unpaired) == 1 && length(empty_xmls) == 1
        i = findfirst(r -> r.xml === nothing, runs)
        r = runs[i]
        runs[i] = _MetronixRunFiles(r.meas_dir, r.run_token, r.rate, only(empty_xmls),
                                    r.ats, r.start_unix, r.n_samples)
        empty!(empty_xmls)
    end
    sort!(runs; by = r -> (r.start_unix, r.rate, r.run_token))
    return runs, sort!(empty_xmls)
end

# The directories to index: `dir` itself when it holds .ats files (a single
# meas_ directory), otherwise its meas_* children (a site), including those
# one level down in rate directories (`DF002.TK/128/meas_*`).
function _metronix_meas_dirs(dir::AbstractString)
    isdir(dir) || return String[]
    _has_ats(dir) && return [_norm_path(dir)]
    root = _norm_path(dir)
    dirs = String[]
    for name in readdir(root; sort = true)
        full = joinpath(root, name)
        if startswith(name, "meas_")
            _has_ats(full) && push!(dirs, full)
        elseif _is_rate_dirname(name) && isdir(full)
            append!(dirs, (joinpath(full, m) for m in readdir(full; sort = true)
                           if startswith(m, "meas_") && _has_ats(joinpath(full, m))))
        end
    end
    return dirs
end

_is_rate_dirname(name::AbstractString) = tryparse(Float64, name) !== nothing

# "128", "4096"; a rate below 1 Hz keeps its fraction: "0.125".
_rate_dirname(rate::Real) = isinteger(rate) ? string(Int(rate)) : string(Float64(rate))

# A raw ADU site: meas_* directories directly inside it, not yet separated by
# rate, and not itself a single meas_ directory.
_is_raw_metronix_site(dir::AbstractString) =
    isdir(dir) && !_has_ats(dir) &&
    any(n -> startswith(n, "meas_") && _has_ats(joinpath(dir, n)), readdir(dir))

_has_ats(dir) = isdir(dir) && any(n -> lowercase(splitext(n)[2]) == ".ats", readdir(dir))

"""
    _metronix_site_index(dir) -> (runs, empty_xmls)

Every run in a Metronix site, or in a single `meas_*` directory, ordered by
start time, plus the XMLs that describe no recorded data.
"""
function _metronix_site_index(dir::AbstractString)
    runs = _MetronixRunFiles[]
    empty_xmls = String[]
    for meas in _metronix_meas_dirs(dir)
        r, e = _metronix_meas_runs(meas)
        append!(runs, r)
        append!(empty_xmls, e)
    end
    sort!(runs; by = r -> (r.start_unix, r.rate, r.run_token))
    return runs, empty_xmls
end

"""
    _metronix_resolve_run(path) -> _MetronixRunFiles

The run `path` names: an `.xml` names the run it describes, an `.ats` the run
it belongs to, and a directory the one run inside it. A directory holding
several runs is ambiguous, and the error lists them.
"""
function _metronix_resolve_run(path::AbstractString)
    p = _norm_path(path)
    if isdir(p)
        runs, _ = _metronix_meas_runs(p)
        isempty(runs) && error("No .ats files found in: $p")
        length(runs) == 1 && return only(runs)
        error("$p holds $(length(runs)) runs; pick the .xml of the one to read:\n" *
              join(("  " * basename(_run_id(r)) for r in runs), "\n"))
    end
    isfile(p) || error("No such file or directory: $p")
    ext = lowercase(splitext(p)[2])
    ext in (".ats", ".xml") || error("Not a Metronix .ats or .xml file: $p")
    runs, empty_xmls = _metronix_meas_runs(dirname(p))
    for r in runs
        (r.xml == p || p in r.ats) && return r
    end
    p in empty_xmls && error("$(basename(p)) describes a run with no recorded .ats data")
    error("No Metronix run found for: $p")
end

function _parse_xml_filename_tokens(xml_path::AbstractString)
    base = splitext(basename(xml_path))[1]
    parts = split(base, '_')
    prefix = isempty(parts) ? "000" : parts[1]
    run_token = length(parts) >= 2 ? parts[end - 1] : "R000"
    freq_token = length(parts) >= 1 ? parts[end] : "128H"
    return String(prefix), String(run_token), String(freq_token)
end

"""
    read_metronix(path; site, components = nothing, include_aux = true) -> TimekeeperRun

Read one Metronix ADU run -- its `.ats` binaries plus their `.xml` sidecar --
into a [`TimekeeperRun`](@ref).

A `meas_*` directory can hold several runs (a 128 Hz run and 4096 Hz bursts,
say), told apart by the run number and rate in their filenames. `path` picks
one: the run's `.xml`, any of its `.ats` files, or a `meas_*` directory that
holds a single run. [`metronix_site_runs`](@ref) lists the runs of a site.

Samples are scaled by each file's LSB value and channel types are mapped
through [`METRONIX_CHANNEL_MAP`](@ref). Pass `components` to read a subset.
`site` defaults to the name of the directory above the `meas_*` directory.

The XML path and filename tokens are kept in metadata so
[`write_metronix`](@ref) can reproduce the original naming.
"""
function read_metronix(path::AbstractString;
                       site = nothing, components = nothing, include_aux = true)
    files = _metronix_resolve_run(path)
    meas_dir = files.meas_dir
    site === nothing && (site = basename(dirname(meas_dir)))
    ats_files = files.ats
    xml_path = files.xml
    xml_path === nothing && error("No .xml file describes run $(files.run_token) at " *
                                  "$(files.rate) Hz in: $meas_dir")
    prefix, run_token, freq_token = _parse_xml_filename_tokens(xml_path)

    channels = Dict{Symbol, TimekeeperChannel}()
    fs = NaN
    start_unix = nothing
    for path in ats_files
        data, info = _read_ats(path)
        ct = info["channel_type"]
        comp = get(METRONIX_CHANNEL_MAP, ct, nothing)
        comp === nothing && (comp = _symbolize(ct))
        if components !== nothing && !(comp in _symbolize.(collect(components)))
            continue
        end
        fs = info["sampling_rate"]
        start_unix = info["start_unix"]
        start_dt = Dates.unix2datetime(info["start_unix"])
        header = Dict{String, Any}(
            "format" => "Metronix-ATS",
            "channel_type" => ct,
            "ats_data_file" => basename(path),
            "ats_header_bytes" => info["header_bytes"],
            "lsbval" => info["lsbval"],
            "start_unix" => info["start_unix"],
            "sample_rate" => fs,
        )
        channels[comp] = TimekeeperChannel(comp, data, fs, start_dt, component_units(comp), path, header)
    end
    isempty(channels) && error("No requested Metronix channels found in: $meas_dir")

    site_dir = dirname(meas_dir)
    metadata = Dict{Symbol, Any}(
        :source_format => :metronix,
        :meas_dir => meas_dir,
        :site_dir => site_dir,
        :metronix_run_id => _run_id(files),
        :metronix_xml_path => xml_path,
        :metronix_prefix => prefix,
        :metronix_run_token => run_token,
        :metronix_freq_token => freq_token,
        :sample_rate => fs,
        :start_time => Dates.unix2datetime(start_unix),
        :n_samples => minimum(length(ch.data) for ch in values(channels)),
        :instrument_model => "Metronix ADU",
        :data_logger_manufacturer => "Metronix",
    )
    return TimekeeperRun(String(site), "Metronix ADU", :metronix, channels, metadata)
end

"""
    load_metronix(meas_dir; components = nothing, kwargs...) -> TimeArray

[`read_metronix`](@ref) followed by [`to_timearray`](@ref) -- the one-step route
when you want a `TimeArray` rather than a run.
"""
function load_metronix(meas_dir::AbstractString; components = nothing, kwargs...)
    run = read_metronix(meas_dir; components = components, kwargs...)
    comps = components === nothing ? default_components(run) : _symbolize.(collect(components))
    return to_timearray(run; components = comps)
end

function _segment_ranges(n::Integer, mask::Union{Nothing, TimekeeperMask}, min_samples::Integer)
    mask === nothing && return UnitRange{Int}[1:n]
    length(mask.masked) == n ||
        error("Mask length $(length(mask.masked)) does not match run length $n")
    ranges = UnitRange{Int}[]
    active = false
    start_index = 1
    for i in 1:n
        if !mask.masked[i] && !active
            active = true
            start_index = i
        elseif mask.masked[i] && active
            (i - 1) - start_index + 1 >= min_samples && push!(ranges, start_index:(i - 1))
            active = false
        end
    end
    active && (n - start_index + 1 >= min_samples) && push!(ranges, start_index:n)
    return ranges
end

_metronix_date_str(dt::DateTime) = Dates.format(dt, "yyyy-mm-dd")
_metronix_time_str(dt::DateTime) = Dates.format(dt, "HH:MM:SS")
_meas_dir_name(dt::DateTime) = "meas_" * Dates.format(dt, "yyyy-mm-dd_HH-MM-SS")

function _metronix_xml_filename(prefix, start_dt::DateTime, stop_dt::DateTime, run_token, freq_token)
    fmt(dt) = Dates.format(dt, "yyyy-mm-dd_HH-MM-SS")
    return "$(prefix)_$(fmt(start_dt))_$(fmt(stop_dt))_$(run_token)_$(freq_token).xml"
end

_set_node_text!(::Nothing, _s) = nothing
function _set_node_text!(node, s::AbstractString)
    node.content = s
    return node
end

# The fields of a run's XML that describe one segment of it: the recording's
# start and stop, each ATSWriter channel's start and sample count, and the
# size of each .ats file. Everything else in the XML stays as the ADU wrote it.
function _set_segment_fields!(doc, start_dt::DateTime, stop_dt::DateTime, n_samples::Integer, file_size::Integer)
    date_s = _metronix_date_str(start_dt)
    time_s = _metronix_time_str(start_dt)
    _set_node_text!(findfirst("//recording/start_date", doc), date_s)
    _set_node_text!(findfirst("//recording/start_time", doc), time_s)
    _set_node_text!(findfirst("//recording/stop_date", doc), _metronix_date_str(stop_dt))
    _set_node_text!(findfirst("//recording/stop_time", doc), _metronix_time_str(stop_dt))
    for ch in findall("//ATSWriter/configuration/channel", doc)
        _set_node_text!(findfirst("./start_date", ch), date_s)
        _set_node_text!(findfirst("./start_time", ch), time_s)
        _set_node_text!(findfirst("./num_samples", ch), string(n_samples))
    end
    _set_node_text!(findfirst("//ATSWriter/output_file/ats_file_size", doc), string(file_size))
    return doc
end

# Set the content of every <tag>…</tag> (or empty <tag/>) inside text[span].
function _xml_text_set(text::String, span::UnitRange{Int}, tag::AbstractString, value::AbstractString)
    inner = replace(SubString(text, first(span), last(span)),
                    Regex("<$(tag)>[^<]*</$(tag)>|<$(tag)/>") => "<$(tag)>$(value)</$(tag)>")
    return SubString(text, 1, prevind(text, first(span))) * inner * SubString(text, nextind(text, last(span)))
end

# The span from an opening <tag> to the first of `ends` after it (exclusive),
# or to the closing </tag>.
function _xml_span(text::String, tag::AbstractString, ends = ())
    open_r = findfirst("<$(tag)>", text)
    open_r === nothing && return nothing
    stop = something(findnext("</$(tag)>", text, last(open_r)), (lastindex(text) + 1):0)
    stop_i = first(stop)
    for e in ends
        r = findnext(e, text, last(open_r))
        r === nothing || (stop_i = min(stop_i, first(r)))
    end
    return first(open_r):prevind(text, stop_i)
end

"""
    _write_segment_xml(out_path, template_path, start_dt, stop_dt, n_samples, header_length)

Write the XML of one segment of a run: the run's own XML with only the fields
that describe the segment changed (see `_set_segment_fields!`). The template's
text is edited in place rather than re-serialised, so every other byte - the
declaration, whitespace, comments, character escapes - is kept as the ADU
wrote it. The edit is checked against the same change made through the XML
tree; if they disagree, or the result does not parse, nothing is written.
"""
function _write_segment_xml(out_path::AbstractString, template_path::AbstractString,
                            start_dt::DateTime, stop_dt::DateTime, n_samples::Integer, header_length::Integer)
    file_size = header_length + n_samples * 4
    template = read(template_path, String)
    date_s = _metronix_date_str(start_dt)
    time_s = _metronix_time_str(start_dt)

    text = template
    # the recording's own fields come before its first child block
    rec = _xml_span(text, "recording", ("<input", "<output", "<ATSWriter"))
    if rec !== nothing
        for (tag, v) in (("start_date", date_s), ("start_time", time_s),
                         ("stop_date", _metronix_date_str(stop_dt)), ("stop_time", _metronix_time_str(stop_dt)))
            text = _xml_text_set(text, rec, tag, v)
            rec = _xml_span(text, "recording", ("<input", "<output", "<ATSWriter"))
        end
    end
    for (tag, v) in (("start_date", date_s), ("start_time", time_s),
                     ("num_samples", string(n_samples)), ("ats_file_size", string(file_size)))
        atsw = _xml_span(text, "ATSWriter")
        atsw === nothing || (text = _xml_text_set(text, atsw, tag, v))
    end

    expected = _set_segment_fields!(EzXML.parsexml(template), start_dt, stop_dt, n_samples, file_size)
    got = try
        EzXML.parsexml(text)
    catch err
        error("Segment XML from $(basename(template_path)) does not parse; not written: " *
              sprint(showerror, err))
    end
    string(got) == string(expected) ||
        error("Segment XML from $(basename(template_path)) would change more than the segment's " *
              "start, stop, sample count and file size; not written")
    write(out_path, text)
    return out_path
end

function _metronix_output_channels(run::TimekeeperRun)
    comps = [c for c in METRONIX_DEFAULT_COMPONENTS if haskey(run.channels, c)]
    isempty(comps) && error("Run has no Metronix channels (e1/e2/bx/by/bz) to write")
    for c in comps
        haskey(run.channels[c].header, "ats_header_bytes") ||
            error("Channel $c is missing the original ATS header; write_metronix requires a run loaded via read_metronix")
    end
    return comps
end

function _write_meas_dir(dest_meas_dir::AbstractString, run::TimekeeperRun, comps,
                         range::UnitRange{Int}, template_path, prefix, run_token, freq_token)
    mkpath(dest_meas_dir)
    fs = sampling_rate(run)
    base_unix = run.channels[first(comps)].header["start_unix"]
    seg_start_unix = _segment_start_unix(base_unix, first(range), fs)
    n_samples = length(range)
    start_dt = Dates.unix2datetime(seg_start_unix)
    # The ADU's stop time is the start plus the whole seconds recorded: two
    # hours at 4096 Hz, 29491200 samples, runs 00:00:00 to 02:00:00.
    stop_dt = Dates.unix2datetime(seg_start_unix + floor(Int, n_samples / fs))

    header_length = 0
    for c in comps
        ch = run.channels[c]
        hb = ch.header["ats_header_bytes"]::Vector{UInt8}
        header_length = length(hb)
        out = joinpath(dest_meas_dir, ch.header["ats_data_file"])
        _write_ats(out, view(ch.data, range), hb, ch.header["lsbval"], seg_start_unix)
    end

    xml_name = _metronix_xml_filename(prefix, start_dt, stop_dt, run_token, freq_token)
    _write_segment_xml(joinpath(dest_meas_dir, xml_name), template_path,
                       start_dt, stop_dt, n_samples, header_length)
    return dest_meas_dir
end

# Sample i of a run starting at base_unix, in whole seconds; segment starts
# are snapped to whole seconds, so this is exact for them.
_segment_start_unix(base_unix::Integer, i::Integer, fs::Real) = base_unix + floor(Int, (i - 1) / fs)

function _snap_range_to_second(range::UnitRange{Int}, sps::Int)
    sps <= 1 && return range
    a = first(range)
    rem = (a - 1) % sps
    a2 = rem == 0 ? a : a + (sps - rem)
    a2 > last(range) && return nothing
    a2 == a || @info "Trimmed $(a2 - a) sample(s) to align segment start to an integer second" segment_start = a
    return a2:last(range)
end

"""
    write_metronix(dest_meas_dir, run::TimekeeperRun) -> String

Write `run` as a Metronix measurement directory: one `.ats` file per channel
plus an `.xml` sidecar derived from the template the run was read with.
Requires `run.metadata[:metronix_xml_path]`, so the run must have come from
[`read_metronix`](@ref). Returns `dest_meas_dir`.
"""
function write_metronix(dest_meas_dir::AbstractString, run::TimekeeperRun)
    comps = _metronix_output_channels(run)
    template = run.metadata[:metronix_xml_path]
    prefix = get(run.metadata, :metronix_prefix, "000")
    run_token = get(run.metadata, :metronix_run_token, "R000")
    freq_token = get(run.metadata, :metronix_freq_token, "128H")
    n = minimum(length(run.channels[c].data) for c in comps)
    return _write_meas_dir(dest_meas_dir, run, comps, 1:n, template, prefix, run_token, freq_token)
end

# Timekeepers' directories beside a site DF002: DF002.TK holds the site
# separated by rate, DF002.TK20260930_141205 one write of it. Each holds one
# directory per rate - 128/, 4096/, 131072/ - of meas_* directories.
const _TK_SUFFIX = ".TK"
const _TK_DIR_RE = r"\.TK(\d{8}_\d{6}(_\d+)?)?$"

# The site directory above a rate directory of a .TK layout, else `dir`.
function _metronix_layout_root(dir::AbstractString)
    d = _norm_path(dir)
    _is_rate_dirname(basename(d)) && occursin(_TK_DIR_RE, basename(dirname(d))) && return dirname(d)
    return d
end

# The source site a layout came from: DF002, DF002.TK, DF002.TK/128 and
# DF002.TK20260930_141205 all give .../DF002.
_metronix_site_stem(dir::AbstractString) = replace(_metronix_layout_root(dir), _TK_DIR_RE => "")

"""
    _tk_write_dir(site_dir) -> String

A fresh destination for one write of `site_dir`: `<site>.TK<date>_<time>`,
e.g. `DF002.TK20260930_141205`, next to the source site. Every write gets its
own directory, so an earlier write is never overwritten; two writes within
the same second get `_2`, `_3`, ... appended.
"""
function _tk_write_dir(site_dir::AbstractString)
    base = _metronix_site_stem(site_dir) * _TK_SUFFIX * Dates.format(Dates.now(), "yyyymmdd_HHMMSS")
    dest = base
    k = 1
    while ispath(dest)
        k += 1
        dest = base * "_$(k)"
    end
    return dest
end

function _tk_site_dir(run::TimekeeperRun)
    site_dir = get(run.metadata, :site_dir, nothing)
    site_dir === nothing && error("Run has no :site_dir metadata; pass dest= explicitly")
    return joinpath(_tk_write_dir(site_dir), _rate_dirname(_rate_key(sampling_rate(run))))
end

"""
    write_metronix_site(run; mask = nothing, dest = nothing, min_samples = 1) -> (String, Vector{String})

Write `run` as a Metronix *site* directory, splitting it at the masked
intervals: each contiguous good stretch becomes its own `meas_<start>`
directory, so downstream tools see clean continuous runs instead of one record
with holes. Segments shorter than `min_samples` are dropped, and segment starts
are trimmed to whole seconds where the rate requires it.

`dest` defaults to the run's rate directory in a new write of its site,
`<site>.TK<date>_<time>/<rate>` (e.g. `DF002.TK20260930_141205/128`). With
`mask = nothing` the run is written whole. Returns the destination directory
and the list of `meas_*` directories created.

For a whole site on disk rather than a single loaded run, see
[`write_metronix_site_masked`](@ref).
"""
function write_metronix_site(run::TimekeeperRun; mask::Union{Nothing, TimekeeperMask} = nothing,
                             dest::Union{Nothing, AbstractString} = nothing, min_samples::Integer = 1)
    comps = _metronix_output_channels(run)
    template = run.metadata[:metronix_xml_path]
    prefix = get(run.metadata, :metronix_prefix, "000")
    run_token = get(run.metadata, :metronix_run_token, "R000")
    freq_token = get(run.metadata, :metronix_freq_token, "128H")
    fs = sampling_rate(run)
    sps = round(Int, fs)

    n = minimum(length(run.channels[c].data) for c in comps)
    ranges = _segment_ranges(n, mask, min_samples)
    isempty(ranges) && error("No unmasked segments to write")

    dest_root = dest === nothing ? _tk_site_dir(run) : String(dest)
    mkpath(dest_root)

    written = String[]
    for range in ranges
        snapped = _snap_range_to_second(range, sps)
        snapped === nothing && continue
        base_unix = run.channels[first(comps)].header["start_unix"]
        seg_start_unix = _segment_start_unix(base_unix, first(snapped), fs)
        meas_name = _meas_dir_name(Dates.unix2datetime(seg_start_unix))
        dest_meas = joinpath(dest_root, meas_name)
        _write_meas_dir(dest_meas, run, comps, snapped, template, prefix, run_token, freq_token)
        push!(written, dest_meas)
    end
    isempty(written) && error("No segments long enough to write (min_samples=$min_samples)")
    @info "Wrote Metronix site" dest = dest_root meas_dirs = length(written)
    return dest_root, written
end

"""
    metronix_site_runs(dir) -> Dict{Float64, Vector{String}}

Every run in a Metronix site, grouped by sampling rate and ordered by start
time. Each run is named by the path of its `.xml`, which
[`read_metronix`](@ref) accepts.

A run is one run number at one rate, not one `meas_*` directory: an ADU can
write a long 128 Hz run and several 4096 Hz bursts into the same directory.
XMLs with no recorded `.ats` data are left out. `dir` is a site of `meas_*`
directories or a single `meas_*` directory. Only the `.ats` headers are read.
"""
function metronix_site_runs(dir::AbstractString)
    isdir(dir) || error("Not a directory: $dir")
    runs = Dict{Float64, Vector{String}}()
    for r in first(_metronix_site_index(dir))
        push!(get!(runs, _rate_key(r.rate), String[]), _run_id(r))
    end
    return runs
end

"""
    metronix_site_rates(dir) -> Vector{Float64}

Sorted unique sampling rates present in a Metronix site (empty if not a site).
"""
function metronix_site_rates(dir::AbstractString)
    isdir(dir) || return Float64[]
    return sort(collect(keys(metronix_site_runs(dir))))
end

"""
    is_metronix_site(dir) -> Bool

Whether `dir` holds Metronix data: `.ats` files of its own, or `meas_*`
directories that do.
"""
is_metronix_site(dir::AbstractString) = !isempty(_metronix_meas_dirs(dir))

_run_files(r::_MetronixRunFiles) = r.xml === nothing ? copy(r.ats) : vcat(r.ats, r.xml)

# `dst` already holds `src`: same size. A file is written whole or not at all
# by cp, so a size match means an earlier split finished it.
_split_done(src::AbstractString, dst::AbstractString) = isfile(dst) && filesize(dst) == filesize(src)

"""
    _metronix_split_plan(site_dir, dest_root) -> Vector{Pair{String, Vector{Pair{String, String}}}}

The copies that separate `site_dir` by rate into `dest_root`, grouped as
`run name => [src => dst, ...]`, in the order they are made. The files beside
the runs of a `meas_*` directory - its `.kml`, XMLs of jobs that never
recorded - are grouped with its last run. Only the `.ats` headers are read.
"""
function _metronix_split_plan(site_dir::AbstractString, dest_root::AbstractString)
    plan = Pair{String, Vector{Pair{String, String}}}[]
    for meas in _metronix_meas_dirs(site_dir)
        runs, _ = _metronix_meas_runs(meas)
        isempty(runs) && continue
        targets = Dict(_rate_key(r.rate) => joinpath(dest_root, _rate_dirname(_rate_key(r.rate)), basename(meas))
                       for r in runs)
        run_files = Set{String}()
        for r in runs
            target = targets[_rate_key(r.rate)]
            push!(plan, _run_id(r) => [f => joinpath(target, basename(f)) for f in _run_files(r)])
            union!(run_files, _run_files(r))
        end
        for name in readdir(meas; sort = true)
            src = joinpath(meas, name)
            (isfile(src) && !(src in run_files)) || continue
            homes = collect(values(targets))
            if lowercase(splitext(name)[2]) == ".xml"
                rate = _xml_rate(src)
                rate !== nothing && haskey(targets, _rate_key(rate)) && (homes = [targets[_rate_key(rate)]])
            end
            append!(last(plan[end]), [src => joinpath(t, name) for t in sort(homes)])
        end
    end
    return plan
end

"""
    metronix_site_is_split(site_dir; dest = nothing) -> Bool

Whether `site_dir` has already been separated by rate into `dest` (default
`<site_dir>.TK`) - every file [`split_metronix_site`](@ref) would copy is
there - so a split would copy nothing.
"""
function metronix_site_is_split(site_dir::AbstractString; dest::Union{Nothing, AbstractString} = nothing)
    site_dir = _norm_path(site_dir)
    dest_root = dest === nothing ? site_dir * _TK_SUFFIX : _norm_path(dest)
    isdir(dest_root) || return false
    return all(_split_done(src, dst) for (_, files) in _metronix_split_plan(site_dir, dest_root)
               for (src, dst) in files)
end

"""
    split_metronix_site(site_dir; dest = nothing, on_run = nothing) -> String

Separate a raw Metronix site by sampling rate into `dest` (default
`<site_dir>.TK`, beside it): one directory per rate, each holding the
`meas_*` directories with runs at that rate.

    DF002/meas_2021-09-25_14-02-01/   ->   DF002.TK/128/meas_2021-09-25_14-02-01/
                                           DF002.TK/4096/meas_2021-09-25_14-02-01/

Each run's `.ats` files and `.xml` go to its rate. An XML of a job that never
recorded goes to the rate in its filename; any other file of a `meas_*`
directory - the `.kml` - goes to every rate the directory holds. Files are
copied byte for byte and the source is not changed. A file already present
at the same size is not copied again, so a site already split is left as it
is ([`metronix_site_is_split`](@ref) tells beforehand). `on_run(i, n,
run_name)` is called before each run that still has files to copy, counting
only those. Returns the destination.
"""
function split_metronix_site(site_dir::AbstractString; dest::Union{Nothing, AbstractString} = nothing,
                             on_run = nothing)
    site_dir = _norm_path(site_dir)
    _is_raw_metronix_site(site_dir) ||
        error("Not a Metronix site of meas_* directories (already separated by rate?): $site_dir")
    dest_root = dest === nothing ? site_dir * _TK_SUFFIX : _norm_path(dest)
    dest_root == site_dir && error("Destination is the site itself: $site_dir")
    todo = [id => [p for p in files if !_split_done(p...)] for (id, files) in _metronix_split_plan(site_dir, dest_root)]
    filter!(t -> !isempty(last(t)), todo)
    for (i, (id, files)) in enumerate(todo)
        on_run === nothing || on_run(i, length(todo), id)
        for (src, dst) in files
            mkpath(dirname(dst))
            cp(src, dst; force = true)
        end
    end
    @info "Separated Metronix site by sampling rate" site = site_dir dest = dest_root runs_copied = length(todo)
    return dest_root
end

function _mask_from_intervals(run::TimekeeperRun, intervals)
    n = run.metadata[:n_samples]
    fs = sampling_rate(run)
    t0 = start_time(run)
    masked = falses(n)
    for iv in intervals
        lo, hi = iv[1] <= iv[2] ? (iv[1], iv[2]) : (iv[2], iv[1])
        lo_s = Dates.value(lo - t0) / 1000.0
        hi_s = Dates.value(hi - t0) / 1000.0
        a = max(1, ceil(Int, lo_s * fs) + 1)
        b = min(n, floor(Int, hi_s * fs) + 1)
        a <= b && (masked[a:b] .= true)
    end
    return TimekeeperMask(DateTime[], masked, Tuple{DateTime, DateTime}[])
end

_rate_label(r::Real) = (isinteger(r) ? string(Int(r)) : string(r)) * " Hz"

# `cuts` pairs each sampling rate with the intervals cut out of its runs.
function _append_mask_history(dest_dir::AbstractString, site_dir::AbstractString,
                              cuts, run_segments)
    path = joinpath(dest_dir, "README.md")
    new_file = !isfile(path)
    open(path, "a") do io
        if new_file
            println(io, "# Timekeepers — manipulated Metronix site")
            println(io)
            println(io, "Source site: `", site_dir, "`")
            println(io)
            println(io, "Each section below records one write (mask/unmask) session,")
            println(io, "with the intervals that were amputated and the runs written.")
            println(io)
        end
        println(io, "## Write session ", Dates.format(Dates.now(), "yyyy-mm-dd HH:MM:SS"))
        println(io)
        cuts = [(r, ivs) for (r, ivs) in cuts if !isempty(ivs)]
        if isempty(cuts)
            println(io, "No intervals masked — runs written unchanged.")
        end
        for (rate, intervals) in cuts
            println(io, "Amputated (masked) intervals", rate === nothing ? "" : " at " * _rate_label(rate), ":")
            println(io)
            println(io, "| # | Start | End | Duration |")
            println(io, "|---|-------|-----|----------|")
            for (i, iv) in enumerate(intervals)
                lo, hi = iv[1] <= iv[2] ? (iv[1], iv[2]) : (iv[2], iv[1])
                dur = Dates.canonicalize(Dates.CompoundPeriod(hi - lo))
                println(io, "| ", i, " | ", lo, " | ", hi, " | ", dur, " |")
            end
            println(io)
        end
        println(io)
        println(io, "Runs written:")
        println(io)
        for (src, segs) in run_segments
            seglist = isempty(segs) ? "_(none — fully masked)_" :
                      join(("`" * replace(relpath(s, dest_dir), '\\' => '/') * "`" for s in segs), ", ")
            println(io, "- `", basename(src), "` → ", length(segs), " segment(s): ", seglist)
        end
        println(io)
    end
    return path
end

"""
    write_metronix_site_masked(site_dir; intervals=[], rate_intervals=Dict(), dest=nothing, min_samples=1, only=nothing)

Write a copy of a Metronix site - every run at every sampling rate - with
masked intervals amputated, into `dest`: by default a new
`<site>.TK<date>_<time>` directory beside the site, e.g.
`DF002.TK20260930_141205`, so no earlier write is overwritten. `site_dir` is a
raw site, the `<site>.TK` directory [`split_metronix_site`](@ref) makes, or an
earlier write; each gives the same destination name.

Intervals are `DateTime` tuples, typically the mask intervals from the app.
`rate_intervals` maps a sampling rate to the intervals cut from its runs;
runs at any other rate are cut by `intervals`. Keep them per rate when the
masks were drawn at one rate: a 128 Hz run recorded through the same hours as
4096 Hz bursts must not lose what was masked on the bursts.

The destination holds one directory per rate, laid out as
[`split_metronix_site`](@ref) lays out the source:

    DF002.TK20260930_141205/128/meas_2021-09-25_14-02-01/
    DF002.TK20260930_141205/4096/meas_2021-09-25_14-02-01/     R001, untouched
    DF002.TK20260930_141205/4096/meas_2021-09-26_00-00-00/     R000 before the cut
    DF002.TK20260930_141205/4096/meas_2021-09-26_00-35-01/     R000 after it

A run no interval touches is copied byte for byte into the `meas_*` directory
it came from, with the `.kml` and the XMLs of jobs at its rate that never
recorded, so a write with no intervals reproduces the site file for file. A
run an interval does touch is split: each unmasked stretch becomes its own
`meas_<start>` directory, with `.ats` files cut from the run, the `.kml`, and
a copy of the run's XML in which only the start, stop, sample counts and file
size are changed. The masked time is simply not recorded - nothing is filled
in. A `README.md` in the destination records the source, the intervals cut at
each rate, and every directory written. Returns the destination directory.

Pass `only` - run names from [`metronix_site_runs`](@ref), or `meas_*`
directories - to write just those runs instead of the whole site.
"""
function write_metronix_site_masked(site_dir::AbstractString;
                                    intervals = Tuple{DateTime, DateTime}[],
                                    rate_intervals::AbstractDict = Dict{Float64, Vector{Tuple{DateTime, DateTime}}}(),
                                    dest::Union{Nothing, AbstractString} = nothing,
                                    min_samples::Integer = 1,
                                    only = nothing)
    site_dir = _norm_path(site_dir)
    runs, _ = _metronix_site_index(site_dir)
    isempty(runs) && error("No Metronix runs found in: $site_dir")
    if only !== nothing
        keep = Set(_norm_path(d) for d in only)
        filter!(r -> _run_id(r) in keep || r.meas_dir in keep, runs)
        isempty(runs) && error("None of the $(length(keep)) requested run(s) are under: $site_dir")
    end

    dest_dir = dest === nothing ? _tk_write_dir(site_dir) : String(dest)
    mkpath(dest_dir)
    run_segments = Tuple{String, Vector{String}}[]
    by_rate = Dict(_rate_key(k) => v for (k, v) in rate_intervals)
    cuts_for(r) = get(by_rate, _rate_key(r.rate), intervals)
    for r in runs
        rate_dir = joinpath(dest_dir, _rate_dirname(_rate_key(r.rate)))
        ivs = cuts_for(r)
        written = _run_touched(r, ivs) ? _write_run_cut(r, ivs, rate_dir, min_samples) : nothing
        if written === nothing
            # untouched: back into the meas_ directory it came from, as it was
            same = joinpath(rate_dir, basename(r.meas_dir))
            _copy_run_verbatim(r, same)
            _copy_meas_extras(r, same; scheduled_xmls = true)
            written = [same]
        else
            foreach(d -> _copy_meas_extras(r, d), written)
        end
        push!(run_segments, (_run_id(r), written))
    end
    cuts = [(rate, cuts_for(first(filter(r -> _rate_key(r.rate) == rate, runs))))
            for rate in sort(unique(_rate_key(r.rate) for r in runs))]
    _append_mask_history(dest_dir, site_dir, cuts, run_segments)
    @info "Wrote manipulated Metronix site" dest = dest_dir
    return dest_dir
end

# Cut run `r` at `intervals` into meas_ directories under rate_dir. `nothing`
# when the intervals mask none of its samples, so the run is copied whole;
# empty when they mask all of it.
function _write_run_cut(r::_MetronixRunFiles, intervals, rate_dir::AbstractString, min_samples::Integer)
    run = read_metronix(_run_id(r))
    mask = _mask_from_intervals(run, intervals)
    any(mask.masked) || return nothing
    try
        return last(write_metronix_site(run; mask = mask, dest = rate_dir, min_samples = min_samples))
    catch err
        msg = sprint(showerror, err)
        occursin("No unmasked segments", msg) || occursin("No segments long enough", msg) || rethrow()
        return String[]
    end
end

# Whether any interval overlaps the run's span, from its header alone.
function _run_touched(r::_MetronixRunFiles, intervals)
    isempty(intervals) && return false
    t0 = Dates.unix2datetime(r.start_unix)
    t1 = t0 + Millisecond(round(Int, 1000 * (r.n_samples - 1) / r.rate))
    for iv in intervals
        lo, hi = iv[1] <= iv[2] ? (iv[1], iv[2]) : (iv[2], iv[1])
        (lo <= t1 && hi >= t0) && return true
    end
    return false
end

function _copy_run_verbatim(r::_MetronixRunFiles, dest_meas::AbstractString)
    mkpath(dest_meas)
    for f in _run_files(r)
        cp(f, joinpath(dest_meas, basename(f)); force = true)
    end
    return dest_meas
end

# Copy the files beside run `r` that belong to no run - the .kml - into
# `target`, a directory written from it. With `scheduled_xmls`, also the XMLs
# of jobs at r's rate that never recorded (any rate, if the XML names none), as
# the rate split places them.
function _copy_meas_extras(r::_MetronixRunFiles, target::AbstractString; scheduled_xmls::Bool = false)
    runs, empty_xmls = _metronix_meas_runs(r.meas_dir)
    run_files = Set(f for x in runs for f in _run_files(x))
    for name in readdir(r.meas_dir; sort = true)
        src = joinpath(r.meas_dir, name)
        (isfile(src) && !(src in run_files)) || continue
        if lowercase(splitext(name)[2]) == ".xml"
            scheduled_xmls && src in empty_xmls || continue
            rate = _xml_rate(src)
            rate === nothing || _rate_key(rate) == _rate_key(r.rate) ||
                !any(x -> _rate_key(x.rate) == _rate_key(rate), runs) || continue
        end
        cp(src, joinpath(target, name); force = true)
    end
    return nothing
end
