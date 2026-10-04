# Utilities.jl - small shared helpers.
# Author: @pankajkmishra
#
# This file has:
# - the table of units for the magnetotelluric channel names of the package
# - the search for the default data directory
# - the helpers that the readers use to normalize paths and symbols

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

The path that the app and the examples use if you do not give a data
directory. It is `data/` next to the package root if that directory exists.
If not, it is `examples/data/`. The path that the function returns possibly
does not exist.
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
