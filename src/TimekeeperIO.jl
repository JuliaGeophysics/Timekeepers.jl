# TimekeeperIO.jl - read and write functions for all formats.
# Author: @pankajkmishra
#
# This file sends read_timekeeper and write_timekeeper to the LEMI-424, GEOMAG
# or Metronix code. If you do not give the format, the file finds it:
# - Metronix: from the path (an .ats file, or a directory that contains one)
# - GEOMAG or LEMI-424: from the first line that is not blank

"""
    read_timekeeper(path; format = :auto, kwargs...) -> TimekeeperRun

Read a supported format into a [`TimekeeperRun`](@ref). The function calls
[`read_lemi424`](@ref), [`read_geomag`](@ref) or [`read_metronix`](@ref).

With `format = :auto`, the function finds the format:
- A `.ats` file, or a directory that contains one, is Metronix.
- A `.txt` file with a `GEOMAG` header is GEOMAG. Other `.txt` files are
  LEMI-424.

To skip the detection, give `format = :lemi424`, `:geomag` or `:metronix`.
The function gives the other keywords to the selected reader.
"""
function read_timekeeper(path::AbstractString; format = :auto, kwargs...)
    fmt = format == :auto ? _detect_format(path) : Symbol(format)
    fmt == :lemi424 && return read_lemi424(path; kwargs...)
    fmt == :geomag && return read_geomag(path; kwargs...)
    fmt == :metronix && return read_metronix(path; kwargs...)
    error("Unsupported Timekeepers format: $fmt")
end

"""
    write_timekeeper(path, run::TimekeeperRun; format = :auto) -> String

Write `run` in its native format. The function calls [`write_lemi424`](@ref),
[`write_geomag`](@ref) or [`write_metronix`](@ref).

With `format = :auto`, the `source_format` of the run selects the writer.
Thus, you can write a file that you read with [`read_timekeeper`](@ref)
without more arguments.
"""
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
