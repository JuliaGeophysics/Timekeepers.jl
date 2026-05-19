function read_timekeeper(path::AbstractString; format = :auto, kwargs...)
    fmt = format == :auto ? _detect_format(path) : Symbol(format)
    if fmt == :lemi424
        return read_lemi424(path; kwargs...)
    elseif fmt in (:metronix, :metronix_ats, :ats)
        return read_metronix(path; kwargs...)
    else
        error("Unsupported Timekeepers format: $fmt")
    end
end

function write_timekeeper(path::AbstractString, run::TimekeeperRun; format = :auto)
    fmt = format == :auto ? _detect_output_format(path, run) : Symbol(format)
    if fmt == :lemi424
        return write_lemi424(path, run)
    elseif fmt in (:metronix, :metronix_ats, :ats)
        return write_metronix(path, run)
    else
        error("Unsupported Timekeepers output format: $fmt")
    end
end

function _detect_format(path::AbstractString)
    isdir(path) && return :metronix_ats
    lower = lowercase(path)
    endswith(lower, ".txt") && return :lemi424
    endswith(lower, ".ats") && return :metronix_ats
    error("Could not infer input format for $path")
end

function _detect_output_format(path::AbstractString, run::TimekeeperRun)
    lower = lowercase(path)
    endswith(lower, ".txt") && return :lemi424
    run.source_format == :lemi424 && return :lemi424
    run.source_format == :metronix_ats && return :metronix_ats
    error("Could not infer output format for $path")
end
