const METRONIX_COMPONENT_MAP = Dict(
    "TEx" => :Ex,
    "TEy" => :Ey,
    "THx" => :Hx,
    "THy" => :Hy,
    "THz" => :Hz,
    "Ex" => :Ex,
    "Ey" => :Ey,
    "Hx" => :Hx,
    "Hy" => :Hy,
    "Hz" => :Hz,
)

function _read_ats_header_fields(io::IO)
    header = Dict{String, Any}()
    header["header_length"] = read(io, UInt16)
    header["header_version"] = read(io, Int16)
    header["sample_length"] = read(io, Int32)
    header["sampling_rate"] = read(io, Float32)
    header["start"] = read(io, Int32)
    header["lsbval"] = read(io, Float64)
    header["GMToffset"] = read(io, Int32)
    header["Res1"] = read(io, Int32)
    header["serial_number_ADU06"] = read(io, UInt16)
    header["serial_number_ADC_board"] = read(io, Int16)
    header["channel_number"] = read(io, Int8)
    header["Res2"] = read(io, Int8)
    header["channel_type"] = _read_cstring(io, 2)
    header["sensor_type"] = _read_cstring(io, 6)
    header["sensor_serial_number"] = read(io, Int16)
    for key in ["x1", "y1", "z1", "x2", "y2", "z2", "E_field_dipole_length", "angle", "rho_probe", "DC_offset_voltage", "internal_gain_ampli"]
        header[key] = read(io, Float32)
    end
    header["Res3"] = read(io, Int32)
    header["ADU_Lat_raw"] = read(io, Int32)
    header["ADU_Long_raw"] = read(io, Int32)
    header["ADU_Elev_raw"] = read(io, Int32)
    header["Lat_Long_TYPE"] = _read_cstring(io, 1)
    header["add_coordinates"] = _read_cstring(io, 1)
    header["ref_meridian"] = read(io, Int16)
    header["xcoor"] = read(io, Float64)
    header["ycoor"] = read(io, Float64)
    header["gps_clock_status"] = _read_cstring(io, 1)
    header["accuracy_GPS"] = read(io, Int8)
    header["offset_UTC"] = read(io, Int16)
    header["Res4a"] = read(io, Int32)
    header["Res4b"] = read(io, Int32)
    header["Res4c"] = read(io, Int32)
    header["survey_header_filename"] = _read_cstring(io, 12)
    header["type_of_meas"] = _read_cstring(io, 4)
    header["logfile"] = _read_cstring(io, 12)
    header["result_selftest"] = _read_cstring(io, 2)
    header["Res5"] = read(io, Int16)
    header["number_of_calib_freq"] = read(io, Int16)
    header["length_of_freq_entry"] = read(io, Int16)
    header["version_calib"] = read(io, Int16)
    header["start_addres"] = read(io, Int16)
    header["Res6"] = read(io, Int64)
    header["cal_filename_ADU06"] = _read_cstring(io, 12)
    header["datetime_calib"] = read(io, Int32)
    header["cal_sensor_filename"] = _read_cstring(io, 12)
    header["datetime_calib_sens"] = read(io, Int32)
    header["powerline1"] = read(io, Float32)
    header["powerline2"] = read(io, Float32)
    header["Res7"] = read(io, Int64)
    header["CSAMT_Tx_freq"] = read(io, Float32)
    header["CSAMT_TS_blocks"] = read(io, Int16)
    header["CSAMT_stacks"] = read(io, Int16)
    header["CSAMT_blk_length"] = read(io, Int32)
    header["Res8"] = read(io, Int32)
    header["Client"] = _read_cstring(io, 16)
    header["Contractor"] = _read_cstring(io, 16)
    header["Area"] = _read_cstring(io, 16)
    header["SurveyID"] = _read_cstring(io, 16)
    header["Operator"] = _read_cstring(io, 16)
    header["Res9"] = read(io, 112)
    header["Weather"] = _read_cstring(io, 64)
    header["Comments"] = _read_cstring(io, 512)

    header["ADU_Lat"] = header["ADU_Lat_raw"] * 1e-3 / 3600
    header["ADU_Long"] = header["ADU_Long_raw"] * 1e-3 / 3600
    header["ADU_Elev"] = header["ADU_Elev_raw"] * 1e-2
    return header
end

function read_ats_header(filename::AbstractString)
    open(filename, "r") do io
        header_length = read(io, UInt16)
        seekstart(io)
        raw = read(io, header_length)
        seekstart(io)
        header = _read_ats_header_fields(io)
        header["_raw_header_bytes"] = raw
        header["_source_file"] = _as_path(filename)
        return header
    end
end

function read_ats_sample(filename::AbstractString)
    header = read_ats_header(filename)
    n = Int(header["sample_length"])
    data = Vector{Int32}(undef, n)
    open(filename, "r") do io
        seek(io, Int(header["header_length"]))
        read!(io, data)
    end
    scaled = Float64.(data) .* Float64(header["lsbval"])
    return scaled, Float64(header["sampling_rate"]), Int(header["start"]), header
end

function _metronix_component(path::AbstractString, header::Dict{String, Any})
    base = basename(path)
    for (token, comp) in METRONIX_COMPONENT_MAP
        occursin(token, base) && return comp
    end
    ch_type = get(header, "channel_type", "")
    ch_type == "Ex" && return :Ex
    ch_type == "Ey" && return :Ey
    ch_type == "Hx" && return :Hx
    ch_type == "Hy" && return :Hy
    ch_type == "Hz" && return :Hz
    channel_number = get(header, "channel_number", 0)
    return Symbol("ch$(channel_number)")
end

function _discover_ats_files(root::AbstractString; frequency = nothing)
    if isfile(root)
        return endswith(lowercase(root), ".ats") ? [_as_path(root)] : String[]
    end
    files = String[]
    for (dir, _, names) in walkdir(root)
        for name in names
            lower = lowercase(name)
            endswith(lower, ".ats") || continue
            if frequency !== nothing
                tag1 = lowercase("$(frequency)H")
                tag2 = lowercase("$(frequency)Hz")
                occursin(tag1, lower) || occursin(tag2, lower) || continue
            end
            push!(files, joinpath(dir, name))
        end
    end
    return sort(files)
end

function _choose_metronix_files(files; measurement = nothing)
    isempty(files) && return files
    measurement === nothing && return files
    selected = filter(path -> occursin(String(measurement), path), files)
    isempty(selected) && error("No ATS files matched measurement=$(measurement)")
    return selected
end

function read_metronix(
    path::AbstractString;
    frequency = nothing,
    site = _site_from_path(path),
    measurement = nothing,
)
    root = _as_path(path)
    files = _choose_metronix_files(_discover_ats_files(root; frequency = frequency); measurement = measurement)
    isempty(files) && error("No Metronix ATS files found in $path")

    channels = Dict{Symbol, TimekeeperChannel}()
    source_paths = Dict{Symbol, String}()
    for file in files
        data, fs, start_unix, header = read_ats_sample(file)
        comp = _metronix_component(file, header)
        haskey(channels, comp) && continue
        header["_relative_path"] = isdir(root) ? relpath(file, root) : basename(file)
        header["_root"] = root
        start = _unix_datetime(start_unix)
        channels[comp] = TimekeeperChannel(comp, data, fs, start, component_units(comp), file, header)
        source_paths[comp] = file
    end

    first_header = first(values(channels)).header
    metadata = Dict{Symbol, Any}(
        :source_root => root,
        :source_files => source_paths,
        :sample_rate => sampling_rate(TimekeeperRun(String(site), "Metronix ADU", :metronix_ats, channels, Dict{Symbol, Any}())),
        :adu_number => extract_adu_number(basename(first(values(source_paths)))),
        :latitude => get(first_header, "ADU_Lat", NaN),
        :longitude => get(first_header, "ADU_Long", NaN),
        :elevation => get(first_header, "ADU_Elev", NaN),
        :instrument_model => "ADU",
        :data_logger_manufacturer => "Metronix",
    )
    return TimekeeperRun(String(site), "Metronix ADU", :metronix_ats, channels, metadata)
end

function load_metronix(path::AbstractString; components = nothing, kwargs...)
    run = read_metronix(path; kwargs...)
    comps = components === nothing ? default_components(run) : components
    return to_timearray(run; components = comps)
end

function _patch_raw_ats_header(raw::Vector{UInt8}, sample_length::Integer)
    patched = copy(raw)
    io = IOBuffer(patched; write = true, read = true)
    seek(io, 4)
    write(io, Int32(sample_length))
    return take!(io)
end

function write_ats_file(filename::AbstractString, data::AbstractVector, header::Dict{String, Any})
    haskey(header, "_raw_header_bytes") ||
        error("ATS header does not contain raw header bytes; read the source ATS first before writing native ATS.")
    lsbval = Float64(header["lsbval"])
    lsbval != 0.0 || error("Cannot write ATS data with lsbval=0")
    raw_data = round.(Int32, Float64.(data) ./ lsbval)
    raw_header = _patch_raw_ats_header(header["_raw_header_bytes"], length(raw_data))
    mkpath(dirname(filename))
    open(filename, "w") do io
        write(io, raw_header)
        write(io, raw_data)
    end
    return filename
end

function write_metronix(output_dir::AbstractString, run::TimekeeperRun; prefix = "")
    mkpath(output_dir)
    written = String[]
    for comp in components(run)
        ch = run.channels[comp]
        rel = get(ch.header, "_relative_path", isempty(ch.source_file) ? "$(comp).ats" : basename(ch.source_file))
        dir = dirname(rel)
        base = prefix * basename(rel)
        output_path = isempty(dir) || dir == "." ? joinpath(output_dir, base) : joinpath(output_dir, dir, base)
        write_ats_file(output_path, ch.data, ch.header)
        push!(written, output_path)
    end
    return written
end

function _metronix_run_from_template(data_run::TimekeeperRun, template_run::TimekeeperRun)
    channels = Dict{Symbol, TimekeeperChannel}()
    for comp in components(data_run)
        haskey(template_run.channels, comp) ||
            error("Template Metronix run does not contain component $comp")
        data_ch = data_run.channels[comp]
        template_ch = template_run.channels[comp]
        isapprox(data_ch.sample_rate, template_ch.sample_rate; atol = 1e-9) ||
            error("Sample rate mismatch for $comp: TimeArray has $(data_ch.sample_rate), template has $(template_ch.sample_rate)")
        channels[comp] = TimekeeperChannel(
            comp,
            data_ch.data,
            template_ch.sample_rate,
            template_ch.start,
            template_ch.units,
            template_ch.source_file,
            copy(template_ch.header),
        )
    end
    metadata = merge(
        copy(template_run.metadata),
        Dict{Symbol, Any}(:template_source => get(template_run.metadata, :source_root, "")),
    )
    return TimekeeperRun(template_run.site, template_run.instrument, :metronix_ats, channels, metadata)
end

function write_metronix(
    output_dir::AbstractString,
    tv::TimeArray,
    ;
    template_run = nothing,
    template_dir = nothing,
    frequency = nothing,
    measurement = nothing,
    site = "unknown",
    kwargs...,
)
    data_run = from_timearray(tv; site = site, instrument = "Metronix ADU", source_format = :timearray, kwargs...)
    if template_run === nothing && template_dir !== nothing
        template_run = read_metronix(template_dir; frequency = frequency, measurement = measurement, site = site)
    end
    template_run === nothing &&
        error("write_metronix requires template_run or template_dir so original ATS headers can be reused.")
    return write_metronix(output_dir, _metronix_run_from_template(data_run, template_run))
end
