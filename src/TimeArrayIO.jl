function _channel_lengths(run::TimekeeperRun, comps)
    return [length(run.channels[c].data) for c in comps]
end

function _time_axis(start::DateTime, n::Integer, fs::Real; axis = :auto)
    n <= 0 && return DateTime[]
    fs > 0 || error("Sample rate must be positive, got $fs")
    if axis == :auto
        axis = fs > 1000 ? :time : :datetime
    end

    if axis == :time
        t0 = Time(start)
        step = Nanosecond(round(Int, 1_000_000_000 / fs))
        return [t0 + (i - 1) * step for i in 1:n]
    elseif axis == :datetime
        return [start + Millisecond(round(Int, (i - 1) * 1000 / fs)) for i in 1:n]
    else
        error("Unsupported time axis $(axis). Use :auto, :datetime, or :time.")
    end
end

function to_timearray(run::TimekeeperRun; components = default_components(run), axis = :auto)
    comps = _symbolize.(collect(components))
    isempty(comps) && error("No components selected")
    missing = [c for c in comps if !haskey(run.channels, c)]
    isempty(missing) || error("Run does not contain components: $(join(missing, ", "))")

    lengths = _channel_lengths(run, comps)
    n = minimum(lengths)
    all(==(n), lengths) || @warn "Components have different lengths; truncating TimeArray to $n samples"

    first_ch = run.channels[first(comps)]
    times = _time_axis(first_ch.start, n, first_ch.sample_rate; axis = axis)
    values_matrix = hcat([run.channels[c].data[1:n] for c in comps]...)
    meta = merge(
        Dict{Symbol, Any}(
            :site => run.site,
            :instrument => run.instrument,
            :source_format => run.source_format,
            :sample_rate => first_ch.sample_rate,
            :start_time => first_ch.start,
            :units => Dict(c => run.channels[c].units for c in comps),
        ),
        run.metadata,
    )
    return TimeArray(times, values_matrix, comps, meta)
end

function _ta_timestamps(ta::TimeArray)
    try
        return timestamp(ta)
    catch
        return getfield(ta, :timestamp)
    end
end

function _ta_values(ta::TimeArray)
    try
        return values(ta)
    catch
        return getfield(ta, :values)
    end
end

function _ta_colnames(ta::TimeArray)
    try
        return colnames(ta)
    catch
        return getfield(ta, :colnames)
    end
end

function _ta_meta(ta::TimeArray)
    try
        return meta(ta)
    catch
        return getfield(ta, :meta)
    end
end

function _sample_rate_from_timearray(ta::TimeArray)
    md = _ta_meta(ta)
    if md isa AbstractDict && haskey(md, :sample_rate)
        return Float64(md[:sample_rate])
    end
    times = _ta_timestamps(ta)
    length(times) < 2 && return 1.0
    dt = times[2] - times[1]
    if dt isa Millisecond
        return 1000.0 / Dates.value(dt)
    elseif dt isa Nanosecond
        return 1_000_000_000.0 / Dates.value(dt)
    elseif dt isa Second
        return 1.0 / Dates.value(dt)
    else
        return 1.0 / (Dates.value(dt) / 1000.0)
    end
end

function _start_from_timearray(ta::TimeArray)
    md = _ta_meta(ta)
    if md isa AbstractDict && haskey(md, :start_time) && md[:start_time] isa DateTime
        return md[:start_time]
    end
    first_time = first(_ta_timestamps(ta))
    if first_time isa DateTime
        return first_time
    elseif first_time isa Time
        date = md isa AbstractDict && haskey(md, :start_date) ? md[:start_date] : Date(1970, 1, 1)
        return DateTime(date) + (first_time - Time(0))
    else
        return DateTime(1970, 1, 1)
    end
end

function from_timearray(
    ta::TimeArray;
    site = "unknown",
    instrument = "unknown",
    source_format = :timearray,
    units = Dict{Symbol, String}(),
    metadata = Dict{Symbol, Any}(),
)
    names = _symbolize.(_ta_colnames(ta))
    vals = _ta_values(ta)
    fs = _sample_rate_from_timearray(ta)
    start = _start_from_timearray(ta)
    channels = Dict{Symbol, TimekeeperChannel}()
    for (i, comp) in enumerate(names)
        ch_units = get(units, comp, component_units(comp))
        channels[comp] = TimekeeperChannel(comp, Float64.(vals[:, i]), fs, start, ch_units, "", Dict{String, Any}())
    end
    md = merge(Dict{Symbol, Any}(:from_timearray => true), metadata)
    return TimekeeperRun(String(site), String(instrument), Symbol(source_format), channels, md)
end
