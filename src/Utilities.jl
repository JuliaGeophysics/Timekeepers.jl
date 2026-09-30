# Utilities.jl - small shared helpers.
# Author: @pankajkmishra
#
# Component-to-unit lookup for the magnetotelluric channel names used across
# the package, the default data directory search, and the path and symbol
# normalisation helpers the readers share.

const MT_COMPONENT_UNITS = Dict(
    :Ex => "mV/km",
    :Ey => "mV/km",
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

"""
    default_data_dir() -> String

Path the app and examples use when no data directory is given: `data/` next to
the package root if it exists, otherwise `examples/data/`. The returned path is
not guaranteed to exist.
"""
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

function _site_from_path(path::AbstractString)
    base = basename(path)
    isempty(base) && return basename(dirname(path))
    return splitext(base)[1]
end
