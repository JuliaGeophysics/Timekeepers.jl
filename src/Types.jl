const ChannelMap = Dict{Symbol, Any}
const MetadataMap = Dict{Symbol, Any}

struct TimekeeperChannel
    component::Symbol
    data::Vector{Float64}
    sample_rate::Float64
    start::DateTime
    units::String
    source_file::String
    header::Dict{String, Any}
end

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

components(run::TimekeeperRun) = sort(collect(keys(run.channels)); by = string)

function default_components(run::TimekeeperRun)
    preferred = [:Ex, :Ey, :Hx, :Hy, :Hz, :bx, :by, :bz, :e1, :e2]
    present = components(run)
    ordered = [c for c in preferred if c in present]
    return isempty(ordered) ? present : ordered
end

function sampling_rate(run::TimekeeperRun)
    isempty(run.channels) && return NaN
    rates = unique(round.(getfield.(collect(values(run.channels)), :sample_rate); digits = 9))
    length(rates) == 1 || error("Run has multiple sample rates: $(join(rates, ", "))")
    return first(rates)
end

function start_time(run::TimekeeperRun)
    isempty(run.channels) && return nothing
    return minimum(ch.start for ch in values(run.channels))
end

function end_time(ch::TimekeeperChannel)
    isempty(ch.data) && return ch.start
    dt_seconds = (length(ch.data) - 1) / ch.sample_rate
    return ch.start + Millisecond(round(Int, dt_seconds * 1000))
end

function end_time(run::TimekeeperRun)
    isempty(run.channels) && return nothing
    return maximum(end_time(ch) for ch in values(run.channels))
end

function duration_seconds(run::TimekeeperRun)
    isempty(run.channels) && return 0.0
    ch = first(values(run.channels))
    return length(ch.data) / ch.sample_rate
end
