function read_timekeeper(path::AbstractString; format = :auto, kwargs...)
    fmt = format == :auto ? _detect_format(path) : Symbol(format)
    fmt == :lemi424 || error("Unsupported Timekeepers format: $fmt")
    return read_lemi424(path; kwargs...)
end

function write_timekeeper(path::AbstractString, run::TimekeeperRun; format = :auto)
    fmt = format == :auto ? _detect_output_format(path, run) : Symbol(format)
    fmt == :lemi424 || error("Unsupported Timekeepers output format: $fmt")
    return write_lemi424(path, run)
end

function _detect_format(path::AbstractString)
    lower = lowercase(path)
    endswith(lower, ".txt") && return :lemi424
    error("Could not infer input format for $path")
end

function _detect_output_format(path::AbstractString, run::TimekeeperRun)
    lower = lowercase(path)
    endswith(lower, ".txt") && return :lemi424
    run.source_format == :lemi424 && return :lemi424
    error("Could not infer output format for $path")
end
