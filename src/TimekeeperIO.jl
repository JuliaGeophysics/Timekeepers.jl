# TimekeeperIO.jl - format-agnostic read/write front end.
# Author: @pankajkmishra
#
# Dispatches read_timekeeper and write_timekeeper to the LEMI-424, GEOMAG or
# Metronix implementation. When the format is not given it is inferred: from
# the path for Metronix (.ats file or a directory containing one), and
# otherwise by sniffing the first non-blank line to tell GEOMAG from LEMI-424.

function read_timekeeper(path::AbstractString; format = :auto, kwargs...)
    fmt = format == :auto ? _detect_format(path) : Symbol(format)
    fmt == :lemi424 && return read_lemi424(path; kwargs...)
    fmt == :geomag && return read_geomag(path; kwargs...)
    fmt == :metronix && return read_metronix(path; kwargs...)
    error("Unsupported Timekeepers format: $fmt")
end

function write_timekeeper(path::AbstractString, run::TimekeeperRun; format = :auto)
    fmt = format == :auto ? _detect_output_format(path, run) : Symbol(format)
    fmt == :lemi424 && return write_lemi424(path, run)
    fmt == :geomag && return write_geomag(path, run)
    fmt == :metronix && return write_metronix(path, run)
    error("Unsupported Timekeepers output format: $fmt")
end

function _is_metronix_dir(path::AbstractString)
    isdir(path) || return false
    for name in readdir(path)
        lowercase(splitext(name)[2]) == ".ats" && return true
    end
    return false
end

function _detect_format(path::AbstractString)
    lower = lowercase(path)
    endswith(lower, ".ats") && return :metronix
    _is_metronix_dir(path) && return :metronix
    if endswith(lower, ".txt")
        detected = open(path, "r") do io
            for line in eachline(io)
                isempty(strip(line)) && continue
                occursin("GEOMAG", uppercase(line)) && return :geomag
                startswith(strip(line), ";") && continue
                break
            end
            return :lemi424
        end
        return detected
    end
    error("Could not infer input format for $path")
end

function _detect_output_format(path::AbstractString, run::TimekeeperRun)
    lower = lowercase(path)
    run.source_format == :geomag && return :geomag
    run.source_format == :lemi424 && return :lemi424
    run.source_format == :metronix && return :metronix
    endswith(lower, ".txt") && return :lemi424
    error("Could not infer output format for $path")
end
