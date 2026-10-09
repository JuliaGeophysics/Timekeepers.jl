# MetronixATS.jl - Metronix ADU binary ATS format reader and writer.
# Author: @pankajkmishra
#
# A Metronix recording is a directory with one .ats file for each channel and
# an XML file that describes the run. An .ats file has a binary header (number
# of samples, rate, start time and LSB scale) and then Int32 samples. A site
# holds many meas_* directories. One meas_* directory can hold more than one
# run. Each run is a run number (R000, R001, ...) at one sampling rate, with
# its own XML
#
# This file reads and writes single runs. It also copies a site into one
# directory for each sampling rate (DF002 -> DF002.128, DF002.4096, ...). In
# these rate directories, it cuts the masked intervals out of the runs:
# - It cuts each run into its good contiguous segments, in the same location.
# - Each segment after a cut becomes a new run number in its meas_* directory.
# - It writes the headers and the XML again to agree with each segment.
# The site itself never changes. Each write adds its cuts to the mask.csv and
# writes the README in the rate directory

"""
    METRONIX_CHANNEL_MAP

The map from Metronix ATS channel types to Timekeepers component names:
`"Ex" => :e1`, `"Ey" => :e2`, `"Hx" => :bx`, `"Hy" => :by`, `"Hz" => :bz`.
A channel type that is not in this table stays as a symbol with the same
name.
"""
const METRONIX_CHANNEL_MAP = Dict("Ex" => :e1, "Ey" => :e2, "Hx" => :bx, "Hy" => :by, "Hz" => :bz)
const METRONIX_DEFAULT_COMPONENTS = [:e1, :e2, :bx, :by, :bz]

const _ATS_OFF_SAMPLE_LENGTH = 4
const _ATS_OFF_SAMPLING_RATE = 8
const _ATS_OFF_START = 12
const _ATS_OFF_LSB = 16
const _ATS_OFF_CHANNEL_TYPE = 38
# ADU-07 and newer headers give the name of the XML of the run here, with NUL
# bytes at the end
const _ATS_OFF_XML_NAME = 448
const _ATS_XML_NAME_LEN = 64

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
        "lsb_mv" => _ats_get(Float64, hbytes, _ATS_OFF_LSB),
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
        data = Float64.(raw) .* info["lsb_mv"]
        return data, info
    end
end

# The XML name in a header, or nothing if the header has no XML name
function _ats_xml_name(bytes::Vector{UInt8})
    length(bytes) >= _ATS_OFF_XML_NAME + _ATS_XML_NAME_LEN || return nothing
    field = bytes[(_ATS_OFF_XML_NAME + 1):(_ATS_OFF_XML_NAME + _ATS_XML_NAME_LEN)]
    stop = something(findfirst(==(0x00), field), length(field) + 1)
    name = String(field[1:(stop - 1)])
    return endswith(lowercase(name), ".xml") ? name : nothing
end

# Change a header to its new XML. If the header has no XML name, do not change
# it
function _set_ats_xml_name!(bytes::Vector{UInt8}, name::AbstractString)
    _ats_xml_name(bytes) === nothing && return bytes
    raw = codeunits(name)
    length(raw) < _ATS_XML_NAME_LEN || error("XML name too long for the ATS header: $name")
    field = zeros(UInt8, _ATS_XML_NAME_LEN)
    field[1:length(raw)] = raw
    bytes[(_ATS_OFF_XML_NAME + 1):(_ATS_OFF_XML_NAME + _ATS_XML_NAME_LEN)] = field
    return bytes
end

function _write_ats(path::AbstractString, data::AbstractVector{<:Real}, header_bytes::Vector{UInt8},
                    lsb_mv::Real, start_unix::Integer; xml_name = nothing)
    bytes = copy(header_bytes)
    _ats_put!(bytes, _ATS_OFF_SAMPLE_LENGTH, Int32(length(data)))
    _ats_put!(bytes, _ATS_OFF_START, Int32(start_unix))
    xml_name === nothing || _set_ats_xml_name!(bytes, xml_name)
    raw = round.(Int32, data ./ lsb_mv)
    open(path, "w") do f
        write(f, bytes)
        write(f, raw)
    end
    return path
end

"""
    _MetronixRunFiles

One recorded run in a `meas_*` directory. It is the `.ats` files with the same
run number (`R000`, `R001`, ...) and sampling rate, and the XML that describes
them. A `meas_*` directory can hold more than one run. An ADU that schedules a
long 128 Hz run and 4096 Hz bursts writes all of them into one directory.
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

# The name of a run is its XML, which is the file that a user selects to load
# it. If the XML is missing, the first .ats file is the name
_run_id(r::_MetronixRunFiles) = something(r.xml, first(r.ats))

_rate_key(fs::Real) = round(Float64(fs); digits = 6)

_norm_path(p::AbstractString) = rstrip(abspath(String(p)), ['/', '\\'])

# "999_V01_C00_R000_TEx_BL_8H.ats" with run R003 is "999_V01_C00_R003_TEx_BL_8H.ats"
function _with_run_token(name::AbstractString, run_token::AbstractString)
    m = match(r"_R\d+_", name)
    if m === nothing
        _ats_run_token(name) == run_token && return String(name)
        error("Cannot renumber $name: its name carries no run number")
    end
    return string(SubString(name, 1, m.offset - 1), "_", run_token, "_",
                  SubString(name, m.offset + length(m.match)))
end

_run_number(run_token::AbstractString) = (m = match(r"^R(\d+)$", run_token); m === nothing ? nothing : parse(Int, m.captures[1]))

function _ats_run_token(path::AbstractString)
    for part in split(splitext(basename(path))[1], '_')
        occursin(r"^R\d+$", part) && return String(part)
    end
    return "R000"
end

# "4096H" is 4096 Hz. "8S" is one sample each 8 s
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

Put the `.ats` files of one `meas_*` directory into runs, by the run number and
the sampling rate in their headers. Then find the XML of each run: its
filename has the same run number and rate. The XMLs without data are jobs
that the ADU scheduled but did not record. The function returns them in
`empty_xmls`. The runs are in order of start time.
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
    # Files with names that do not obey the convention: if one run and one XML
    # stay without a pair, they are one run, as before the function could tell
    # runs apart
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

# The directories for the index: `dir` itself if it holds .ats files (one
# meas_ directory). If not, the meas_* directories in it (a site)
function _metronix_meas_dirs(dir::AbstractString)
    isdir(dir) || return String[]
    _has_ats(dir) && return [_norm_path(dir)]
    root = _norm_path(dir)
    return [joinpath(root, n) for n in readdir(root; sort = true)
            if startswith(n, "meas_") && _has_ats(joinpath(root, n))]
end

# "128", "4096". A rate below 1 Hz keeps its fraction: "0.125"
_rate_dirname(rate::Real) = isinteger(rate) ? string(Int(rate)) : string(Float64(rate))

# A site of meas_* directories, which is not one meas_ directory itself
_is_raw_metronix_site(dir::AbstractString) =
    isdir(dir) && !_has_ats(dir) &&
    any(n -> startswith(n, "meas_") && _has_ats(joinpath(dir, n)), readdir(dir))

_has_ats(dir) = isdir(dir) && any(n -> lowercase(splitext(n)[2]) == ".ats", readdir(dir))

"""
    _metronix_site_index(dir) -> (runs, empty_xmls)

All the runs in a Metronix site, or in one `meas_*` directory, in order of
start time. The function also returns the XMLs that describe no recorded
data.
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

The run that `path` identifies:
- an `.xml` identifies the run that it describes;
- an `.ats` identifies the run that it is part of;
- a directory identifies the one run in it.

A directory that holds more than one run is ambiguous. The error then gives a
list of the runs.
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

# The rate token of the ADU: 4096 Hz is "4096H", one sample each 8 s is "8S"
function _freq_token(rate::Real)
    rate >= 1 && return isinteger(rate) ? "$(Int(rate))H" : "$(rate)H"
    period = 1 / rate
    return isinteger(round(period; digits = 6)) ? "$(round(Int, period))S" : "$(period)S"
end

# The name tokens of a run that has no XML. The prefix and the run number come
# from the .ats filename ("406_V01_C00_R000_TEx_BL_128H.ats"). The rate token
# comes from the header. Thus, a filename that is not usual cannot give an
# incorrect rate
function _ats_filename_tokens(files::_MetronixRunFiles)
    prefix = split(splitext(basename(first(files.ats)))[1], '_')[1]
    return String(prefix), files.run_token, _freq_token(files.rate)
end

"""
    read_metronix(path; site, components = nothing, include_aux = true) -> TimekeeperRun

Read one Metronix ADU run (its `.ats` binary files and their `.xml` sidecar)
into a [`TimekeeperRun`](@ref). The XML is optional. All the data that the
function reads comes from the `.ats` headers. Thus, the function can read a
run without its XML. It then gives a warning, and [`write_metronix`](@ref)
writes only the `.ats` files of the run.

A `meas_*` directory can hold more than one run, for example a 128 Hz run and
4096 Hz bursts. The run number and the rate in the filenames identify each
run. `path` selects one run: the `.xml` of the run, one of its `.ats` files, or
a `meas_*` directory that holds one run. [`metronix_site_runs`](@ref) gives
the runs of a site.

The function multiplies the samples by the LSB value of each file. It changes
the channel types with [`METRONIX_CHANNEL_MAP`](@ref). To read only some
components, give `components`. The default `site` is the name of the
directory above the `meas_*` directory.

The metadata keeps the XML path and the filename tokens. Thus,
[`write_metronix`](@ref) can use the original names again.
"""
function read_metronix(path::AbstractString;
                       site = nothing, components = nothing, include_aux = true)
    files = _metronix_resolve_run(path)
    meas_dir = files.meas_dir
    site === nothing && (site = basename(dirname(meas_dir)))
    ats_files = files.ats
    xml_path = files.xml
    if xml_path === nothing
        @warn "No .xml file describes run $(files.run_token) at $(files.rate) Hz; reading the .ats files alone" meas_dir maxlog = 1
        prefix, run_token, freq_token = _ats_filename_tokens(files)
    else
        prefix, run_token, freq_token = _parse_xml_filename_tokens(xml_path)
    end

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
            "lsb_mv" => info["lsb_mv"],
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

[`read_metronix`](@ref) and then [`to_timearray`](@ref). Use this function to
get a `TimeArray` in one step, not a run.
"""
function load_metronix(meas_dir::AbstractString; components = nothing, kwargs...)
    run = read_metronix(meas_dir; components = components, kwargs...)
    comps = components === nothing ? default_components(run) : _symbolize.(collect(components))
    return to_timearray(run; components = comps)
end

function _segment_ranges(n::Integer, mask::Union{Nothing, TimekeeperMask})
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
            push!(ranges, start_index:(i - 1))
            active = false
        end
    end
    active && push!(ranges, start_index:n)
    return ranges
end

_metronix_date_str(dt::DateTime) = Dates.format(dt, "yyyy-mm-dd")
_metronix_time_str(dt::DateTime) = Dates.format(dt, "HH:MM:SS")

function _metronix_xml_filename(prefix, start_dt::DateTime, stop_dt::DateTime, run_token, freq_token)
    fmt(dt) = Dates.format(dt, "yyyy-mm-dd_HH-MM-SS")
    return "$(prefix)_$(fmt(start_dt))_$(fmt(stop_dt))_$(run_token)_$(freq_token).xml"
end

_set_node_text!(::Nothing, _s) = nothing
function _set_node_text!(node, s::AbstractString)
    node.content = s
    return node
end

# The fields of the XML of a run that describe one segment of it:
# - the start and the stop of the recording
# - the start, the number of samples and the .ats filename of each ATSWriter
#   channel
# - the size of each .ats file
# All the other data in the XML stays as the ADU wrote it
function _set_segment_fields!(doc, start_dt::DateTime, stop_dt::DateTime, n_samples::Integer, file_size::Integer,
                              run_token::AbstractString)
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
        f = findfirst("./ats_data_file", ch)
        f === nothing || _set_node_text!(f, _with_run_token(strip(f.content), run_token))
    end
    _set_node_text!(findfirst("//ATSWriter/output_file/ats_file_size", doc), string(file_size))
    return doc
end

# Set the content of each <tag>…</tag> (or empty <tag/>) in text[span]
function _xml_text_set(text::String, span::UnitRange{Int}, tag::AbstractString, value::AbstractString)
    inner = replace(SubString(text, first(span), last(span)),
                    Regex("<$(tag)>[^<]*</$(tag)>|<$(tag)/>") => "<$(tag)>$(value)</$(tag)>")
    return SubString(text, 1, prevind(text, first(span))) * inner * SubString(text, nextind(text, last(span)))
end

# The span from an opening <tag> to the first of `ends` after it (not
# included), or to the closing </tag>
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
    _write_segment_xml(out_path, template_path, start_dt, stop_dt, n_samples, header_length, run_token)

Write the XML of one segment of a run. It is the XML of the run, with only the
fields that describe the segment changed (refer to `_set_segment_fields!`).
Its `.ats` names have `run_token`. The function edits the text of the
template in its location. It does not serialize the XML again. Thus, all the
other bytes stay as the ADU wrote them: the declaration, the whitespace, the
comments and the character escapes. The function compares the edit with the
same change made through the XML tree. If the two are different, or if the
result does not parse, the function writes nothing.
"""
function _write_segment_xml(out_path::AbstractString, template_path::AbstractString,
                            start_dt::DateTime, stop_dt::DateTime, n_samples::Integer, header_length::Integer,
                            run_token::AbstractString)
    file_size = header_length + n_samples * 4
    template = read(template_path, String)
    date_s = _metronix_date_str(start_dt)
    time_s = _metronix_time_str(start_dt)

    text = template
    # the fields of the recording come before its first child block
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
    atsw = _xml_span(text, "ATSWriter")
    if atsw !== nothing
        inner = replace(SubString(text, first(atsw), last(atsw)),
                        r"<ats_data_file>[^<]*</ats_data_file>" =>
                        t -> "<ats_data_file>" * _with_run_token(t[16:(end - 16)], run_token) * "</ats_data_file>")
        text = SubString(text, 1, prevind(text, first(atsw))) * inner * SubString(text, nextind(text, last(atsw)))
    end

    expected = _set_segment_fields!(EzXML.parsexml(template), start_dt, stop_dt, n_samples, file_size, run_token)
    got = try
        EzXML.parsexml(text)
    catch err
        error("Segment XML from $(basename(template_path)) does not parse; not written: " *
              sprint(showerror, err))
    end
    string(got) == string(expected) ||
        error("Segment XML from $(basename(template_path)) would change more than the segment's " *
              "start, stop, sample count, file size and .ats names; not written")
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
    # The stop time of the ADU is the start plus the whole seconds recorded. For
    # example, two hours at 4096 Hz (29491200 samples) go from 00:00:00 to
    # 02:00:00
    stop_dt = Dates.unix2datetime(seg_start_unix + floor(Int, n_samples / fs))

    # write a run that has no XML as it came: only its .ats files
    xml_name = template_path === nothing ? nothing :
               _metronix_xml_filename(prefix, start_dt, stop_dt, run_token, freq_token)
    header_length = 0
    for c in comps
        ch = run.channels[c]
        hb = ch.header["ats_header_bytes"]::Vector{UInt8}
        header_length = length(hb)
        out = joinpath(dest_meas_dir, _with_run_token(ch.header["ats_data_file"], run_token))
        ispath(out) && error("Not overwriting $out")
        _write_ats(out, view(ch.data, range), hb, ch.header["lsb_mv"], seg_start_unix; xml_name = xml_name)
    end
    xml_name === nothing && return dest_meas_dir
    _write_segment_xml(joinpath(dest_meas_dir, xml_name), template_path,
                       start_dt, stop_dt, n_samples, header_length, run_token)
    return dest_meas_dir
end

# The time of sample i of a run that starts at base_unix, in whole seconds. The
# starts of segments are moved to whole seconds. Thus, this time is exact for
# them
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

Write `run` as a Metronix measurement directory. The directory gets one `.ats`
file for each channel and an `.xml` sidecar. The sidecar comes from the
template that the reader used. The run must come from
[`read_metronix`](@ref). If the run has no XML, the function writes it without
an XML and gives a warning. The function returns `dest_meas_dir`.
"""
function write_metronix(dest_meas_dir::AbstractString, run::TimekeeperRun)
    comps = _metronix_output_channels(run)
    template = _metronix_template(run)
    prefix = get(run.metadata, :metronix_prefix, "000")
    run_token = get(run.metadata, :metronix_run_token, "R000")
    freq_token = get(run.metadata, :metronix_freq_token, "128H")
    n = minimum(length(run.channels[c].data) for c in comps)
    return _write_meas_dir(dest_meas_dir, run, comps, 1:n, template, prefix, run_token, freq_token)
end

function _metronix_template(run::TimekeeperRun)
    template = get(run.metadata, :metronix_xml_path, nothing)
    template === nothing &&
        @warn "Run $(run.site) was read without an XML; writing its .ats files without one" maxlog = 1
    return template
end

# The work on a site DF002 occurs in its rate directories DF002.128,
# DF002.4096, ... next to it. Each rate directory is a site that holds the
# runs at one rate

# The rate directory of a site: DF002 at 128 Hz is DF002.128
_metronix_rate_dir(site_dir::AbstractString, rate::Real) =
    _norm_path(site_dir) * "." * _rate_dirname(_rate_key(rate))

# The site that `dir` came from when the site was split by rate (DF002 for
# DF002.128). If `dir` is not a rate directory of a site next to it, the
# result is nothing
function _metronix_split_source(dir::AbstractString)
    d = _norm_path(dir)
    name = basename(d)
    for i in findall(==('.'), name)
        suffix = name[(i + 1):end]
        rate = tryparse(Float64, suffix)
        (rate === nothing || rate <= 0 || _rate_dirname(_rate_key(rate)) != suffix) && continue
        stem = joinpath(dirname(d), name[1:(i - 1)])
        _is_raw_metronix_site(stem) && return stem
    end
    return nothing
end

"""
    write_metronix_site(run; mask = nothing, dest = nothing, min_samples = 1) -> (String, Vector{String})

Cut `run` at the masked intervals. Write the parts into `dest`, a `meas_*`
directory, in place of the files of the run there.

- The default `dest` is the `meas_*` directory of the run in its rate
  directory. For example, it is `DF002.128/meas_2021-09-25_14-02-01` for a
  128 Hz run of `DF002`. A run that the reader read from there is replaced in
  the same location.
- The first good part keeps the run number. Each subsequent part gets the
  next free run number at the rate of the run in `dest`.
- If the rate needs it, the function moves the start of a segment to a whole
  second.
- Some segments are too short to keep as a run: shorter than one whole second
  (the ATS header and the XML cannot describe them), or shorter than
  `min_samples`. The function ignores such a segment and gives a warning.
- With `mask = nothing`, the function writes the full run.

The function returns `dest` and the runs that it wrote. The name of each run
is a name that [`read_metronix`](@ref) accepts. For a full site, use
[`write_metronix_site_masked`](@ref).
"""
function write_metronix_site(run::TimekeeperRun; mask::Union{Nothing, TimekeeperMask} = nothing,
                             dest::Union{Nothing, AbstractString} = nothing, min_samples::Integer = 1)
    if dest === nothing
        meas_dir = get(run.metadata, :meas_dir, nothing)
        meas_dir === nothing && error("Run has no :meas_dir metadata; pass dest= explicitly")
        site = dirname(_norm_path(meas_dir))
        if _metronix_split_source(site) === nothing
            # a run of the site itself: first make its full rate directory
            _is_raw_metronix_site(site) && split_metronix_site(site)
            site = _metronix_rate_dir(site, sampling_rate(run))
        end
        dest = joinpath(site, basename(meas_dir))
    end
    dest = _norm_path(dest)
    written, kept = _replace_run!(dest, run, mask, min_samples)
    rate_dir = dirname(dest)
    site = _metronix_split_source(rate_dir)
    site === nothing || _record_cuts!(rate_dir, site, [(run, kept)])
    @info "Wrote Metronix run" dest runs = length(written)
    return dest, written
end

# The good parts of `run` that can be stored as runs, with their starts moved
# to whole seconds. The ADU stores the start and the stop in whole seconds.
# Thus, the function ignores each part shorter than one second (a run with its
# stop equal to its start) or shorter than `min_samples`, and gives a warning
function _metronix_pieces(run::TimekeeperRun, comps, mask, min_samples::Integer)
    fs = sampling_rate(run)
    sps = round(Int, fs)
    n = minimum(length(run.channels[c].data) for c in comps)
    ranges = _segment_ranges(n, mask)
    isempty(ranges) && error("No unmasked segments to write")
    shortest = max(min_samples, ceil(Int, fs))
    pieces = UnitRange{Int}[]
    for range in ranges
        snapped = _snap_range_to_second(range, sps)
        if snapped === nothing || length(snapped) < shortest
            seg_start = run.channels[first(comps)].start + Millisecond(round(Int, 1000 * (first(range) - 1) / fs))
            @warn "Skipped a segment too short to store as a Metronix run" site = run.site rate = fs start = seg_start samples = length(range) seconds = round(length(range) / fs; digits = 3) minimum_samples = shortest
            continue
        end
        push!(pieces, snapped)
    end
    isempty(pieces) && error("No segments long enough to write (minimum $shortest samples)")
    return pieces
end

"""
    metronix_site_runs(dir) -> Dict{Float64, Vector{String}}

All the runs in a Metronix site, in groups by sampling rate, in order of start
time. The name of each run is the path of its `.xml`, which
[`read_metronix`](@ref) accepts.

A run is one run number at one rate. It is not one `meas_*` directory: an ADU
can write a long 128 Hz run and some 4096 Hz bursts into the same directory.
The result does not include XMLs without recorded `.ats` data. `dir` is a site
of `meas_*` directories or one `meas_*` directory. The function reads only the
`.ats` headers.
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

The different sampling rates in a Metronix site, in increasing order. The list
is empty if `dir` is not a site.
"""
function metronix_site_rates(dir::AbstractString)
    isdir(dir) || return Float64[]
    return sort(collect(keys(metronix_site_runs(dir))))
end

"""
    is_metronix_site(dir) -> Bool

Tells if `dir` holds Metronix data: its own `.ats` files, or `meas_*`
directories with `.ats` files.
"""
is_metronix_site(dir::AbstractString) = !isempty(_metronix_meas_dirs(dir))

_run_files(r::_MetronixRunFiles) = r.xml === nothing ? copy(r.ats) : vcat(r.ats, r.xml)

"""
    _metronix_split_plan(site_dir) -> Dict{Float64, Vector{Pair{String, Vector{Pair{String, String}}}}}

The copies that put the runs of `site_dir` into its rate directories, for each
rate. The groups have the form `run name => [src => dst, ...]`, in the order
of the copies. `dst` is relative to the rate directory. The other files of a
`meas_*` directory (its `.kml`, and the XMLs of jobs that did not record) go
in the group of its last run at each rate that they go to. The function reads
only the `.ats` headers.
"""
function _metronix_split_plan(site_dir::AbstractString)
    plan = Dict{Float64, Vector{Pair{String, Vector{Pair{String, String}}}}}()
    for meas in _metronix_meas_dirs(site_dir)
        runs, _ = _metronix_meas_runs(meas)
        isempty(runs) && continue
        m = basename(meas)
        last_of = Dict{Float64, Vector{Pair{String, String}}}()
        run_files = Set{String}()
        for r in runs
            k = _rate_key(r.rate)
            files = [f => joinpath(m, basename(f)) for f in _run_files(r)]
            push!(get!(plan, k, Pair{String, Vector{Pair{String, String}}}[]), _run_id(r) => files)
            last_of[k] = files
            union!(run_files, _run_files(r))
        end
        for name in readdir(meas; sort = true)
            src = joinpath(meas, name)
            (isfile(src) && !(src in run_files)) || continue
            homes = collect(keys(last_of))
            if lowercase(splitext(name)[2]) == ".xml"
                rate = _xml_rate(src)
                rate !== nothing && haskey(last_of, _rate_key(rate)) && (homes = [_rate_key(rate)])
            end
            foreach(k -> push!(last_of[k], src => joinpath(m, name)), homes)
        end
    end
    return plan
end

"""
    metronix_site_is_split(site_dir; dest = nothing) -> Bool

Tells if `site_dir` is already split by rate: each rate of the site has its
rate directory in `dest` (default: next to the site). The function makes a
rate directory only as a full directory. Thus, a rate directory that exists
is complete.
"""
function metronix_site_is_split(site_dir::AbstractString; dest::Union{Nothing, AbstractString} = nothing)
    site_dir = _norm_path(site_dir)
    parent = dest === nothing ? dirname(site_dir) : _norm_path(dest)
    return all(r -> isdir(joinpath(parent, basename(_metronix_rate_dir(site_dir, r)))), metronix_site_rates(site_dir))
end

"""
    split_metronix_site(site_dir; dest = nothing, on_run = nothing) -> Vector{String}

Copy a Metronix site into one directory for each sampling rate. The name of
each directory is `<site>.<rate>`, and it is in `dest` (default: next to the
site). Each directory is a usual site of `meas_*` directories with the runs at
that rate:

    DF002/meas_2021-09-25_14-02-01/   ->   DF002.128/meas_2021-09-25_14-02-01/
                                           DF002.4096/meas_2021-09-25_14-02-01/

The `.ats` files and the `.xml` of each run go to its rate. The XML of a job
that did not record goes to the rate in its filename. Each other file of a
`meas_*` directory (the `.kml`) goes to each rate that the directory holds.
The function copies the files byte for byte and does not change the source.

[`write_metronix_site_masked`](@ref) cuts the masked intervals in the rate
directories. Thus, the function never copies over a rate directory that
exists. It makes only the missing rates. It copies each rate directory fully
under a temporary name and then renames it. Thus, a split that stops before
the end leaves no partial rate directory. To start a rate again from the site,
delete its rate directory. The function calls `on_run(i, n, run_name)` before
it copies each run. It returns the rate directories, the lowest rate first.
"""
function split_metronix_site(site_dir::AbstractString; dest::Union{Nothing, AbstractString} = nothing,
                             on_run = nothing)
    site_dir = _norm_path(site_dir)
    _is_raw_metronix_site(site_dir) ||
        error("Not a Metronix site of meas_* directories: $site_dir")
    _metronix_split_source(site_dir) === nothing ||
        error("Already one rate of a site split by rate: $site_dir")
    parent = dest === nothing ? dirname(site_dir) : _norm_path(dest)
    plan = _metronix_split_plan(site_dir)
    target(k) = joinpath(parent, basename(_metronix_rate_dir(site_dir, k)))
    todo = sort([k for k in keys(plan) if !isdir(target(k))])
    n = sum(k -> length(plan[k]), todo; init = 0)
    i = 0
    for k in todo
        partial = target(k) * ".partial"
        rm(partial; force = true, recursive = true)
        for (id, files) in plan[k]
            i += 1
            on_run === nothing || on_run(i, n, id)
            for (src, rel) in files
                mkpath(dirname(joinpath(partial, rel)))
                cp(src, joinpath(partial, rel); force = true)
            end
        end
        mv(partial, target(k))
    end
    rate_dirs = [target(k) for k in sort(collect(keys(plan)))]
    @info "Separated Metronix site by sampling rate" site = site_dir rate_dirs rates_made = length(todo)
    return rate_dirs
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

const _MASK_CSV_HEADER = "start_sample,end_sample,start_time,end_time"

_iso_ms(dt::DateTime) = Dates.format(dt, dateformat"yyyy-mm-ddTHH:MM:SS.sss")

# The start of `run` in the run that the site recorded, in samples: zero for a
# run without a cut, more for a part after a previous cut
function _recorded_offset(site::AbstractString, run::TimekeeperRun)
    meas = joinpath(site, basename(get(run.metadata, :meas_dir, "")))
    isdir(meas) || return 0
    fs = sampling_rate(run)
    start_unix = first(values(run.channels)).header["start_unix"]
    for r in first(_metronix_meas_runs(meas))
        _rate_key(r.rate) == _rate_key(fs) || continue
        offset = round(Int, (start_unix - r.start_unix) * fs)
        0 <= offset < r.n_samples && return offset
    end
    return 0
end

# The sample ranges of 1:n that are not in `kept`
function _removed_ranges(n::Integer, kept)
    out = UnitRange{Int}[]
    next = 1
    for k in sort(kept; by = first)
        first(k) > next && push!(out, next:(first(k) - 1))
        next = last(k) + 1
    end
    next <= n && push!(out, next:n)
    return out
end

# Add the parts removed from each cut run to the mask.csv of the rate
# directory, and write its README
function _record_cuts!(rate_dir::AbstractString, site::AbstractString, cut_runs)
    path = joinpath(rate_dir, "mask.csv")
    rows = isfile(path) ? [l for l in readlines(path)[2:end] if !isempty(strip(l))] : String[]
    for (run, kept) in cut_runs
        n = run.metadata[:n_samples]
        fs = sampling_rate(run)
        t0 = start_time(run)
        offset = _recorded_offset(site, run)
        at(i) = t0 + Millisecond(floor(Int, 1000 * (i - 1) / fs))
        for r in _removed_ranges(n, kept)
            push!(rows, join((offset + first(r), offset + last(r), _iso_ms(at(first(r))), _iso_ms(at(last(r)))), ","))
        end
    end
    sort!(unique!(rows); by = l -> (split(l, ',')[3], parse(Int, split(l, ',')[1])))
    open(path, "w") do io
        println(io, _MASK_CSV_HEADER)
        foreach(l -> println(io, l), rows)
    end
    label = _rate_label(sampling_rate(first(first(cut_runs))))
    write(joinpath(rate_dir, "README.md"), """
        # $(basename(rate_dir))

        The $(label) runs of `$(basename(site))` with masked stretches cut out by Timekeepers; `mask.csv` lists them.
        `$(basename(site))` itself is never changed: delete this directory to start again from it.
        """)
    return path
end

_rate_label(r::Real) = (isinteger(r) ? string(Int(r)) : string(r)) * " Hz"

"""
    write_metronix_site_masked(site_dir; intervals=[], rate_intervals=Dict(), min_samples=1, only=nothing, format=:default) -> Vector{String}

Cut the masked intervals out of a Metronix site. The cuts occur in the rate
directories of the site, `<site>.128`, `<site>.4096`, ... next to it. If these
directories are missing, [`split_metronix_site`](@ref) first makes them. The
site itself never changes. `site_dir` is the site or one of its rate
directories. A rate directory stands for the full site.

The intervals are `DateTime` tuples, usually the mask intervals from the app.
`rate_intervals` is a map from a sampling rate to the intervals for its runs.
`intervals` cuts the runs at all the other rates. Keep the intervals for each
rate if you drew the masks at one rate. For example, a 128 Hz run that
recorded in the same hours as 4096 Hz bursts must not lose the parts that you
masked on the bursts.

If no interval touches a run, the run stays as it is. If an interval touches a
run, the good parts of the run replace it in its `meas_*` directory. The first
part keeps the run number. Each subsequent part gets the next free run number
at that rate in that directory, after all the runs and scheduled jobs there.
Thus, the function does not rename other files:

    DF002.4096/meas_2021-09-26_00-00-00/   R000 up to the cut, R002 after it,
                                           R001 untouched, the .kml

The function cuts the `.ats` files of a part from the run. The XML of the part
is the XML of the run. Only the start, the stop, the numbers of samples and
the file size change, and the `.ats` names when the part gets a new number.
Each `.ats` header that gives the name of its XML changes to the new XML. The
`.kml` stays. The masked time is not in the data, and nothing fills it.

Some parts are too short to keep as a run: shorter than one whole second, or
shorter than `min_samples`. The function ignores such a part and gives a
warning. If no run stays in a `meas_*` directory, the function removes the
directory. If a run had no XML, the function writes it without an XML. The
function returns the rate directories that it cut.

Each rate directory with a cut gets a `mask.csv`. The file gives each removed
part, one row for each part. The rows of all the writes stay in the file, in
time order:

    start_sample,end_sample,start_time,end_time
    17,40,2025-04-01T07:00:08.000,2025-04-01T07:00:10.875

The samples have numbers from 1 at the start of the run, as the site recorded
it. The two ends are included. Thus, the numbers stay the same after more
cuts. The times are the times of the first and the last removed sample. The
file includes the parts moved to a whole second and the parts that were too
short. Thus, the file describes the data as written. A short `README.md` next
to the file tells where the directory came from.

`format = :default` gives the layout above. `:MTH5` is reserved for MTH5
output, which is not available yet. To cut only some runs, give `only`: run
names from [`metronix_site_runs`](@ref), or `meas_*` directories, in the rate
directories.
"""
function write_metronix_site_masked(site_dir::AbstractString;
                                    intervals = Tuple{DateTime, DateTime}[],
                                    rate_intervals::AbstractDict = Dict{Float64, Vector{Tuple{DateTime, DateTime}}}(),
                                    min_samples::Integer = 1,
                                    only = nothing,
                                    format::Symbol = :default)
    format === :MTH5 && error("MTH5 output is not available yet")
    format === :default || error("format must be :default or :MTH5, not :$format")
    site_dir = _norm_path(site_dir)
    site = something(_metronix_split_source(site_dir), site_dir)
    _is_raw_metronix_site(site) || error("Not a Metronix site of meas_* directories: $site")
    split_metronix_site(site)
    by_rate = Dict(_rate_key(k) => v for (k, v) in rate_intervals)
    cuts_for(rate) = get(by_rate, _rate_key(rate), intervals)
    keep = only === nothing ? nothing : Set(_norm_path(d) for d in only)

    written_dirs = String[]
    for rate in metronix_site_rates(site)
        ivs = cuts_for(rate)
        isempty(ivs) && continue
        rate_dir = _metronix_rate_dir(site, rate)
        runs, _ = _metronix_site_index(rate_dir)
        keep === nothing || filter!(r -> _run_id(r) in keep || r.meas_dir in keep, runs)
        filter!(r -> _run_touched(r, ivs), runs)
        isempty(runs) && continue
        cut_runs = Tuple{TimekeeperRun, Vector{UnitRange{Int}}}[]
        for r in runs
            run = read_metronix(_run_id(r))
            mask = _mask_from_intervals(run, ivs)
            any(mask.masked) || continue
            _, kept = _replace_run!(r.meas_dir, run, mask, min_samples)
            push!(cut_runs, (run, kept))
        end
        isempty(cut_runs) && continue
        for meas in unique(r.meas_dir for r in runs)
            _has_ats(meas) && continue
            @warn "Every run in $(basename(meas)) was masked; removing it" meas
            rm(meas; recursive = true)
        end
        _record_cuts!(rate_dir, site, cut_runs)
        push!(written_dirs, rate_dir)
    end
    @info "Cut masked intervals from Metronix site" site rate_dirs = written_dirs
    return written_dirs
end

# Tells if an interval overlaps the span of the run. Uses only the header
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

# For each meas_ directory and rate, the highest run number of a run or a
# scheduled XML there. New runs get numbers after it
function _next_run_numbers(runs, empty_xmls)
    used = Dict{Tuple{String, Float64}, Int}()
    any_rate = Dict{String, Int}()
    for r in runs
        k = (basename(r.meas_dir), _rate_key(r.rate))
        used[k] = max(get(used, k, -1), something(_run_number(r.run_token), -1))
    end
    for x in empty_xmls
        _, token, _ = _parse_xml_filename_tokens(x)
        num = something(_run_number(token), -1)
        rate = _xml_rate(x)
        meas = basename(dirname(x))
        if rate === nothing
            any_rate[meas] = max(get(any_rate, meas, -1), num)
        else
            k = (meas, _rate_key(rate))
            used[k] = max(get(used, k, -1), num)
        end
    end
    return (; used, any_rate)
end

function _take_run_number!(next_run, meas::AbstractString, rate::Real)
    k = (meas, _rate_key(rate))
    num = max(get(next_run.used, k, -1), get(next_run.any_rate, meas, -1)) + 1
    next_run.used[k] = num
    return "R" * lpad(num, 3, '0')
end

# Replace the files of `run` in the meas_ directory `dest` with its good parts.
# The first part gets the number of the run. Each subsequent part gets the
# next free number at its rate there. The function first writes the parts to
# a staging directory. Thus, it removes the run only when the replacement is
# complete. It returns the runs that it wrote (empty if the mask leaves
# nothing to store) and the sample ranges of `run` in them
function _replace_run!(dest::AbstractString, run::TimekeeperRun, mask, min_samples::Integer)
    comps = _metronix_output_channels(run)
    pieces = try
        _metronix_pieces(run, comps, mask, min_samples)
    catch err
        msg = sprint(showerror, err)
        occursin("No unmasked segments", msg) || occursin("No segments long enough", msg) || rethrow()
        UnitRange{Int}[]
    end
    template = _metronix_template(run)
    prefix = get(run.metadata, :metronix_prefix, "000")
    run_token = get(run.metadata, :metronix_run_token, "R000")
    freq_token = get(run.metadata, :metronix_freq_token, "128H")
    fs = sampling_rate(run)
    mkpath(dest)
    here, empty_xmls = _metronix_meas_runs(dest)
    old = filter(x -> x.run_token == run_token && _rate_key(x.rate) == _rate_key(fs), here)
    next_run = _next_run_numbers(here, empty_xmls)
    if template !== nothing
        # the XML of the run will be replaced: keep its text for the edits
        staged_template = tempname()
        cp(template, staged_template)
        template = staged_template
    end
    staging = mktempdir(dirname(dest); prefix = ".tk_staging_")
    try
        tokens = String[]
        for (i, range) in enumerate(pieces)
            token = i == 1 ? run_token : _take_run_number!(next_run, basename(dest), fs)
            _write_meas_dir(staging, run, comps, range, template, prefix, token, freq_token)
            push!(tokens, token)
        end
        foreach(x -> foreach(f -> rm(f), _run_files(x)), old)
        for name in readdir(staging)
            mv(joinpath(staging, name), joinpath(dest, name); force = false)
        end
        here, _ = _metronix_meas_runs(dest)
        ids = [_run_id(only(filter(x -> x.run_token == t && _rate_key(x.rate) == _rate_key(fs), here))) for t in tokens]
        return ids, pieces
    finally
        rm(staging; force = true, recursive = true)
        template === nothing || rm(template; force = true)
    end
end
