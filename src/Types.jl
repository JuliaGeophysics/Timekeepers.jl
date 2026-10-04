# Types.jl - core data model.
# Author: @pankajkmishra
#
# This file defines the two containers. Each reader makes them, and each
# writer uses them:
# - TimekeeperChannel: the samples of one component, with its rate, start and
#   header
# - TimekeeperRun: a named set of channels, with the run metadata
# It also defines the accessors for the components, sample rate, time span and
# duration of a run

const MetadataMap = Dict{Symbol, Any}

"""
    TimekeeperChannel

One component of a recording. It holds the samples and all the data that you
need to put them in time and to interpret them.

# Fields
- `component::Symbol` -- the standard component name (`:bx`, `:by`, `:bz`, `:e1`, `:e2`, ...).
- `data::Vector{Float64}` -- the samples, in `units`.
- `sample_rate::Float64` -- samples per second.
- `start::DateTime` -- timestamp of `data[1]`.
- `units::String` -- physical units, e.g. `"nT"` or `"mV/km"`.
- `source_file::String` -- the file that the samples came from (empty for synthetic data).
- `header::Dict{String, Any}` -- the header fields of the format, kept for the writer.

See also [`TimekeeperRun`](@ref), [`end_time`](@ref).
"""
struct TimekeeperChannel
    component::Symbol
    data::Vector{Float64}
    sample_rate::Float64
    start::DateTime
    units::String
    source_file::String
    header::Dict{String, Any}
end

"""
    TimekeeperRun

One continuous recording. It is a set of [`TimekeeperChannel`](@ref)s with the
same sample rate and start time. It also holds the run metadata that a writer
needs to make the original file again.

# Fields
- `site::String` -- site name, usually from the file or directory name.
- `instrument::String` -- instrument description, e.g. `"Metronix ADU"`.
- `source_format::Symbol` -- `:lemi424`, `:geomag` or `:metronix`.
- `channels::Dict{Symbol, TimekeeperChannel}` -- channels, with the component as key.
- `metadata::Dict{Symbol, Any}` -- run metadata: position, sample rate, header
  values and, for Metronix, the paths of the XML templates for the writer.

The readers ([`read_timekeeper`](@ref), [`read_lemi424`](@ref),
[`read_geomag`](@ref), [`read_metronix`](@ref)) return a run. The writers use
it. Use [`to_timearray`](@ref) to change it to a `TimeSeries.TimeArray`.
"""
struct TimekeeperRun
    site::String
    instrument::String
    source_format::Symbol
    channels::Dict{Symbol, TimekeeperChannel}
    metadata::MetadataMap
end

function Base.show(io::IO, run::TimekeeperRun)
    comps = join(string.(components(run)), ", ")
    print(
        io,
        "TimekeeperRun(site=\"$(run.site)\", instrument=\"$(run.instrument)\", ",
        "format=:$(run.source_format), components=[$comps])",
    )
end

"""
    components(run::TimekeeperRun) -> Vector{Symbol}

All the component names in `run`, in alphabetical order.

See also [`default_components`](@ref).
"""
components(run::TimekeeperRun) = sort(collect(keys(run.channels)); by = string)

"""
    default_components(run::TimekeeperRun) -> Vector{Symbol}

The components of `run` in the usual magnetotelluric plot order
(`bx, by, bz, e1, e2`). Components that are not in the run are not in the
list. If the run has none of these names, the function returns
[`components`](@ref).
"""
function default_components(run::TimekeeperRun)
    preferred = [:bx, :by, :bz, :e1, :e2, :Bx, :By, :Bz, :Ex, :Ey]
    present = components(run)
    ordered = [c for c in preferred if c in present]
    return isempty(ordered) ? present : ordered
end

"""
    sampling_rate(run::TimekeeperRun) -> Float64

The sample rate of `run` in Hz. If the channels have different rates, the
function gives an error, because a `TimekeeperRun` holds one rate. To
separate a Metronix site with more than one rate, use
[`metronix_site_rates`](@ref) and the procedure that splits a site by rate.
"""
function sampling_rate(run::TimekeeperRun)
    isempty(run.channels) && return NaN
    rates = unique(round.(getfield.(collect(values(run.channels)), :sample_rate); digits = 9))
    length(rates) == 1 || error("Run has multiple sample rates: $(join(rates, ", "))")
    return first(rates)
end

"""
    start_time(run::TimekeeperRun) -> Union{DateTime, Nothing}

The earliest start of the channels in `run`, or `nothing` if the run has no
channels.
"""
function start_time(run::TimekeeperRun)
    isempty(run.channels) && return nothing
    return minimum(ch.start for ch in values(run.channels))
end

"""
    end_time(ch::TimekeeperChannel) -> DateTime
    end_time(run::TimekeeperRun) -> Union{DateTime, Nothing}

The timestamp of the last sample. The function calculates it from the start,
the number of samples and the sample rate. For a run, it is the latest end of
its channels.
"""
function end_time(ch::TimekeeperChannel)
    isempty(ch.data) && return ch.start
    dt_seconds = (length(ch.data) - 1) / ch.sample_rate
    return ch.start + Millisecond(round(Int, dt_seconds * 1000))
end

function end_time(run::TimekeeperRun)
    isempty(run.channels) && return nothing
    return maximum(end_time(ch) for ch in values(run.channels))
end

"""
    duration_seconds(run::TimekeeperRun) -> Float64

The length of `run` in seconds (`nsamples / sample_rate`), or `0.0` if the run
is empty.
"""
function duration_seconds(run::TimekeeperRun)
    isempty(run.channels) && return 0.0
    ch = first(values(run.channels))
    return length(ch.data) / ch.sample_rate
end
