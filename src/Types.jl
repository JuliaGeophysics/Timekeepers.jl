# Types.jl - core data model.
# Author: @pankajkmishra
#
# Defines the two containers every reader produces and every writer consumes:
# TimekeeperChannel (one component's samples plus its rate, start and header)
# and TimekeeperRun (a named set of channels plus run metadata), along with the
# accessors for a run's components, sample rate, time span and duration.

const ChannelMap = Dict{Symbol, Any}
const MetadataMap = Dict{Symbol, Any}

"""
    TimekeeperChannel

One component of a recording: the samples plus everything needed to place them
in time and interpret them.

# Fields
- `component::Symbol` -- canonical component name (`:bx`, `:by`, `:bz`, `:e1`, `:e2`, ...).
- `data::Vector{Float64}` -- the samples, in `units`.
- `sample_rate::Float64` -- samples per second.
- `start::DateTime` -- timestamp of `data[1]`.
- `units::String` -- physical units, e.g. `"nT"` or `"mV/km"`.
- `source_file::String` -- file the samples were read from (empty if synthesised).
- `header::Dict{String, Any}` -- format-specific header fields kept for round-tripping.

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

A single continuous recording: a set of [`TimekeeperChannel`](@ref)s sharing a
sample rate and start time, plus the run-level metadata a writer needs to
reproduce the original file.

# Fields
- `site::String` -- site name, usually derived from the file or directory name.
- `instrument::String` -- instrument description, e.g. `"Metronix ADU"`.
- `source_format::Symbol` -- `:lemi424`, `:geomag` or `:metronix`.
- `channels::Dict{Symbol, TimekeeperChannel}` -- channels keyed by component.
- `metadata::Dict{Symbol, Any}` -- run metadata (position, sample rate, header
  values, and for Metronix the XML template paths used on write).

Readers ([`read_timekeeper`](@ref), [`read_lemi424`](@ref),
[`read_geomag`](@ref), [`read_metronix`](@ref)) return one of these; writers
consume it. Use [`to_timearray`](@ref) to move to a `TimeSeries.TimeArray`.
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

All component names present in `run`, sorted alphabetically.

See also [`default_components`](@ref).
"""
components(run::TimekeeperRun) = sort(collect(keys(run.channels)); by = string)

"""
    default_components(run::TimekeeperRun) -> Vector{Symbol}

The components of `run` in conventional magnetotelluric plotting order
(`bx, by, bz, e1, e2`), skipping any that are absent. Falls back to
[`components`](@ref) when none of the preferred names are present.
"""
function default_components(run::TimekeeperRun)
    preferred = [:bx, :by, :bz, :e1, :e2, :Bx, :By, :Bz, :Ex, :Ey]
    present = components(run)
    ordered = [c for c in preferred if c in present]
    return isempty(ordered) ? present : ordered
end

"""
    sampling_rate(run::TimekeeperRun) -> Float64

Sample rate of `run` in Hz. Errors if the channels disagree, since a
`TimekeeperRun` is meant to hold one rate; use [`metronix_site_rates`](@ref)
and the split-by-rate workflow to separate mixed-rate Metronix sites.
"""
function sampling_rate(run::TimekeeperRun)
    isempty(run.channels) && return NaN
    rates = unique(round.(getfield.(collect(values(run.channels)), :sample_rate); digits = 9))
    length(rates) == 1 || error("Run has multiple sample rates: $(join(rates, ", "))")
    return first(rates)
end

"""
    start_time(run::TimekeeperRun) -> Union{DateTime, Nothing}

Earliest channel start in `run`, or `nothing` when the run has no channels.
"""
function start_time(run::TimekeeperRun)
    isempty(run.channels) && return nothing
    return minimum(ch.start for ch in values(run.channels))
end

"""
    end_time(ch::TimekeeperChannel) -> DateTime
    end_time(run::TimekeeperRun) -> Union{DateTime, Nothing}

Timestamp of the last sample, computed from the start, the sample count and the
sample rate. For a run this is the latest end across its channels.
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

Length of `run` in seconds (`nsamples / sample_rate`), or `0.0` when empty.
"""
function duration_seconds(run::TimekeeperRun)
    isempty(run.channels) && return 0.0
    ch = first(values(run.channels))
    return length(ch.data) / ch.sample_rate
end
