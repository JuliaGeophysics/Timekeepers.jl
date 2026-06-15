const METRONIX_CHANNEL_MAP = Dict("Ex" => :e1, "Ey" => :e2, "Hx" => :bx, "Hy" => :by, "Hz" => :bz)
const METRONIX_COMPONENT_TO_CHTYPE = Dict(v => k for (k, v) in METRONIX_CHANNEL_MAP)
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

function _metronix_files(meas_dir::AbstractString)
    isdir(meas_dir) || error("Not a directory: $meas_dir")
    ats = String[]
    xml = String[]
    for name in readdir(meas_dir; sort = true)
        full = joinpath(meas_dir, name)
        isfile(full) || continue
        ext = lowercase(splitext(name)[2])
        ext == ".ats" && push!(ats, full)
        ext == ".xml" && push!(xml, full)
    end
    isempty(ats) && error("No .ats files found in: $meas_dir")
    isempty(xml) && error("No .xml file found in: $meas_dir")
    return ats, first(xml)
end

function _parse_xml_filename_tokens(xml_path::AbstractString)
    base = splitext(basename(xml_path))[1]
    parts = split(base, '_')
    prefix = isempty(parts) ? "000" : parts[1]
    run_token = length(parts) >= 2 ? parts[end - 1] : "R000"
    freq_token = length(parts) >= 1 ? parts[end] : "128H"
    return String(prefix), String(run_token), String(freq_token)
end

function read_metronix(meas_dir::AbstractString;
                       site = basename(rstrip(dirname(abspath(meas_dir)), ['/', '\\'])),
                       components = nothing, include_aux = true)
    meas_dir = abspath(meas_dir)
    ats_files, xml_path = _metronix_files(meas_dir)
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

function _write_segment_xml(out_path::AbstractString, template_path::AbstractString,
                            start_dt::DateTime, stop_dt::DateTime, n_samples::Integer, header_length::Integer)
    doc = EzXML.readxml(template_path)
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
    _set_node_text!(findfirst("//ATSWriter/output_file/ats_file_size", doc),
                    string(header_length + n_samples * 4))

    EzXML.write(out_path, doc)
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
    sps = round(Int, fs)
    base_unix = run.channels[first(comps)].header["start_unix"]
    seg_start_unix = base_unix + (first(range) - 1) ÷ max(sps, 1)
    n_samples = length(range)
    start_dt = Dates.unix2datetime(seg_start_unix)
    stop_dt = Dates.unix2datetime(seg_start_unix + (n_samples - 1) ÷ max(sps, 1))

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

function _snap_range_to_second(range::UnitRange{Int}, sps::Int)
    sps <= 1 && return range
    a = first(range)
    rem = (a - 1) % sps
    a2 = rem == 0 ? a : a + (sps - rem)
    a2 > last(range) && return nothing
    a2 == a || @info "Trimmed $(a2 - a) sample(s) to align segment start to an integer second" segment_start = a
    return a2:last(range)
end

function write_metronix(dest_meas_dir::AbstractString, run::TimekeeperRun)
    comps = _metronix_output_channels(run)
    template = run.metadata[:metronix_xml_path]
    prefix = get(run.metadata, :metronix_prefix, "000")
    run_token = get(run.metadata, :metronix_run_token, "R000")
    freq_token = get(run.metadata, :metronix_freq_token, "128H")
    n = minimum(length(run.channels[c].data) for c in comps)
    return _write_meas_dir(dest_meas_dir, run, comps, 1:n, template, prefix, run_token, freq_token)
end

const _TK_WRITE_SUFFIX = ".W"

_tk_write_dir(site_dir::AbstractString) = rstrip(String(site_dir), ['/', '\\']) * _TK_WRITE_SUFFIX

function _tk_site_dir(run::TimekeeperRun)
    site_dir = get(run.metadata, :site_dir, nothing)
    site_dir === nothing && error("Run has no :site_dir metadata; pass dest= explicitly")
    return _tk_write_dir(site_dir)
end

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
        seg_start_unix = base_unix + (first(snapped) - 1) ÷ max(sps, 1)
        meas_name = _meas_dir_name(Dates.unix2datetime(seg_start_unix))
        dest_meas = joinpath(dest_root, meas_name)
        _write_meas_dir(dest_meas, run, comps, snapped, template, prefix, run_token, freq_token)
        push!(written, dest_meas)
    end
    isempty(written) && error("No segments long enough to write (min_samples=$min_samples)")
    @info "Wrote Metronix site" dest = dest_root meas_dirs = length(written)
    return dest_root, written
end

_rate_label(fs::Real) = string(round(Int, fs))

function _meas_dir_rate(meas_dir::AbstractString)
    for name in readdir(meas_dir; sort = true)
        lowercase(splitext(name)[2]) == ".ats" || continue
        return _read_ats_header(joinpath(meas_dir, name))["sampling_rate"]
    end
    return nothing
end

"""
    metronix_site_runs(site_dir) -> Dict{Float64, Vector{String}}

Scan a Metronix site directory and group its `meas_*` run directories by the
sampling rate detected from each run's ATS header (one rate per `meas_*` dir).
"""
function metronix_site_runs(site_dir::AbstractString)
    isdir(site_dir) || error("Not a directory: $site_dir")
    runs = Dict{Float64, Vector{String}}()
    for name in readdir(site_dir; sort = true)
        d = joinpath(site_dir, name)
        (isdir(d) && startswith(name, "meas_")) || continue
        rate = _meas_dir_rate(d)
        rate === nothing && continue
        push!(get!(runs, round(rate; digits = 6), String[]), d)
    end
    return runs
end

"""
    metronix_site_rates(site_dir) -> Vector{Float64}

Sorted unique sampling rates present in a Metronix site (empty if not a site).
"""
function metronix_site_rates(site_dir::AbstractString)
    isdir(site_dir) || return Float64[]
    return sort(collect(keys(metronix_site_runs(site_dir))))
end

is_metronix_site(dir::AbstractString) = !isempty(metronix_site_rates(dir))

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

function _append_mask_history(dest_dir::AbstractString, site_dir::AbstractString,
                              intervals, run_segments)
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
        if isempty(intervals)
            println(io, "No intervals masked — runs written unchanged.")
        else
            println(io, "Amputated (masked) intervals:")
            println(io)
            println(io, "| # | Start | End | Duration |")
            println(io, "|---|-------|-----|----------|")
            for (i, iv) in enumerate(intervals)
                lo, hi = iv[1] <= iv[2] ? (iv[1], iv[2]) : (iv[2], iv[1])
                dur = Dates.canonicalize(Dates.CompoundPeriod(hi - lo))
                println(io, "| ", i, " | ", lo, " | ", hi, " | ", dur, " |")
            end
        end
        println(io)
        println(io, "Runs written:")
        println(io)
        for (src, segs) in run_segments
            seglist = isempty(segs) ? "_(none — fully masked)_" :
                      join(("`" * basename(s) * "`" for s in segs), ", ")
            println(io, "- `", basename(src), "` → ", length(segs), " segment(s): ", seglist)
        end
        println(io)
    end
    return path
end

"""
    write_metronix_site_masked(site_dir; intervals=[], dest=nothing, min_samples=1)

Apply amputation `intervals` (DateTime tuples, typically the mask intervals
from the app) to every `meas_*` run of a single-rate Metronix site and write
the remaining contiguous segments as `meas_*` directories under `dest`
(default `<site_dir>.W`, e.g. `RK137.TK128` → `RK137.TK128.W`). Runs untouched
by any interval are written whole. A `README.md` mask/unmask history (one
datetime-stamped section per write) is appended in the destination. Assumes the
site contains a single sampling rate (use the splitter script to separate rates
first). Returns the destination directory.
"""
function write_metronix_site_masked(site_dir::AbstractString;
                                    intervals = Tuple{DateTime, DateTime}[],
                                    dest::Union{Nothing, AbstractString} = nothing,
                                    min_samples::Integer = 1)
    site_dir = rstrip(abspath(site_dir), ['/', '\\'])
    runs_by_rate = metronix_site_runs(site_dir)
    isempty(runs_by_rate) && error("No Metronix meas_ directories found in: $site_dir")
    length(runs_by_rate) > 1 &&
        @warn "Site has multiple sampling rates; writing all into one dir. Split by rate first." rates = sort(collect(keys(runs_by_rate)))

    dest_dir = dest === nothing ? _tk_write_dir(site_dir) : String(dest)
    mkpath(dest_dir)
    run_segments = Tuple{String, Vector{String}}[]
    for (_, meas_dirs) in sort(collect(runs_by_rate); by = first)
        for meas in sort(meas_dirs)
            run = read_metronix(meas)
            mask = isempty(intervals) ? nothing : _mask_from_intervals(run, intervals)
            _, written = write_metronix_site(run; mask = mask, dest = dest_dir, min_samples = min_samples)
            push!(run_segments, (meas, written))
        end
    end
    _append_mask_history(dest_dir, site_dir, intervals, run_segments)
    @info "Wrote manipulated Metronix site" dest = dest_dir
    return dest_dir
end
