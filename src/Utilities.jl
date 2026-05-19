const MT_COMPONENT_UNITS = Dict(
    :Ex => "mV/km",
    :Ey => "mV/km",
    :Hx => "nT",
    :Hy => "nT",
    :Hz => "nT",
    :bx => "nT",
    :by => "nT",
    :bz => "nT",
    :e1 => "mV/km",
    :e2 => "mV/km",
    :e3 => "mV/km",
    :e4 => "mV/km",
    :temperature_e => "C",
    :temperature_h => "C",
    :battery => "V",
    :elevation => "m",
    :latitude => "degrees",
    :longitude => "degrees",
    :n_satellites => "count",
    :gps_fix => "flag",
    :time_diff => "s",
)

component_units(component::Symbol) = get(MT_COMPONENT_UNITS, component, "")

function default_data_dir()
    root_data = normpath(joinpath(@__DIR__, "..", "data"))
    isdir(root_data) && return root_data
    examples_data = normpath(joinpath(@__DIR__, "..", "examples", "data"))
    isdir(examples_data) && return examples_data
    return root_data
end

function _as_path(path)
    return abspath(String(path))
end

function _symbolize(x)
    x isa Symbol && return x
    return Symbol(String(x))
end

function _stringify_metadata(metadata::Dict{Symbol, Any})
    return Dict(String(k) => v for (k, v) in metadata)
end

function _symbol_metadata(metadata::AbstractDict)
    return Dict(Symbol(k) => v for (k, v) in metadata)
end

function _read_cstring(io::IO, n::Integer)
    bytes = Vector{UInt8}(undef, n)
    read!(io, bytes)
    stop = findfirst(==(0x00), bytes)
    last_index = stop === nothing ? n : stop - 1
    return String(bytes[1:last_index])
end

function _padded_bytes(value, n::Integer)
    out = zeros(UInt8, n)
    bytes = Vector{UInt8}(String(value))
    count = min(length(bytes), n)
    count > 0 && copyto!(out, 1, bytes, 1, count)
    return out
end

function _unix_datetime(seconds)
    return Dates.unix2datetime(round(Int, seconds))
end

function format_duration(seconds)
    hours = Int(floor(seconds / 3600))
    minutes = Int(floor((seconds - hours * 3600) / 60))
    secs = Int(round(seconds - hours * 3600 - minutes * 60))
    return @sprintf("%02d:%02d:%02d", hours, minutes, secs)
end

function extract_adu_number(filename::AbstractString)
    m = match(r"^(\d+)_", basename(filename))
    return m === nothing ? "Unknown" : m.captures[1]
end

function _site_from_path(path::AbstractString)
    base = basename(path)
    isempty(base) && return basename(dirname(path))
    return splitext(base)[1]
end

function _component_header(header::Dict{String, Any}, key::String, default)
    return get(header, key, default)
end

function _copy_with_new_data(ch::TimekeeperChannel, data::AbstractVector)
    return TimekeeperChannel(
        ch.component,
        Float64.(data),
        ch.sample_rate,
        ch.start,
        ch.units,
        ch.source_file,
        copy(ch.header),
    )
end
