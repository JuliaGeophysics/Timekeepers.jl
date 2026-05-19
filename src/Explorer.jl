const TK_BLACK = RGBf(0.05, 0.05, 0.07)
const TK_BLUE = RGBf(0.114, 0.306, 0.847)
const TK_GREY = RGBf(0.42, 0.45, 0.50)
const TK_MUTED = RGBf(0.62, 0.66, 0.71)
const TK_FRAME = RGBf(0.18, 0.20, 0.24)
const TK_PANEL_BG = RGBf(0.985, 0.987, 0.990)
const TK_MASK_FILL = RGBAf(0.40, 0.43, 0.50, 0.18)
const TK_SEL_FILL = RGBAf(0.114, 0.306, 0.847, 0.18)
const TK_SEL_EDGE = RGBAf(0.114, 0.306, 0.847, 0.70)

const WINDOW_OPTIONS = [
    ("1 minute",  60.0),
    ("10 minutes", 600.0),
    ("30 minutes", 1800.0),
    ("1 hour",    3600.0),
    ("6 hours",   21600.0),
    ("12 hours",  43200.0),
    ("1 day",     86400.0),
    ("3 days",    259200.0),
    ("7 days",    604800.0),
    ("All",       Inf),
]

mutable struct TKApp
    data::TimeArray
    mask::TimekeeperMask
    figure::Figure
    plot_layout::GridLayout
    summary_label::Label
    status_label::Label
    axes::Vector{Axis}
    origin::DateTime
    time_seconds::Vector{Float64}
    span_seconds::Float64
    raw_values::Matrix{Float64}
    line_clean::Vector{Observable{Vector{Float32}}}
    line_masked::Vector{Observable{Vector{Float32}}}
    line_x::Observable{Vector{Float64}}
    window_seconds::Observable{Float64}
    window_start::Observable{Float64}
    slider::Any
    window_menu::Any
    selection::Observable{Tuple{Float64, Float64}}
    selection_visible::Observable{Bool}
    mask_lows::Observable{Vector{Float64}}
    mask_highs::Observable{Vector{Float64}}
    source_format::Symbol
    source_path::String
end

function _component_color(name)
    s = lowercase(String(name))
    return startswith(s, "e") ? TK_BLUE : TK_BLACK
end

function _ensure_datetime(times)
    first(times) isa DateTime && return collect(DateTime, times)
    first(times) isa Date && return [DateTime(d) for d in times]
    error("Unsupported time axis element type $(eltype(times)). Need DateTime or Date.")
end

function _seconds_since(origin::DateTime, times::Vector{DateTime})
    return Float64[Dates.value(t - origin) / 1000.0 for t in times]
end

function _nice_step(target::Float64)
    nice = (1.0, 2.0, 5.0, 10.0, 30.0, 60.0, 120.0, 300.0, 600.0, 1800.0,
            3600.0, 7200.0, 10800.0, 21600.0, 43200.0, 86400.0, 172800.0,
            432000.0, 604800.0, 1209600.0, 2592000.0)
    target <= 0 && return nice[1]
    return nice[argmin(abs.(log.(nice) .- log(target)))]
end

function _datetime_xticks_for_range(origin::DateTime, x_lo::Float64, x_hi::Float64; target = 6)
    span = x_hi - x_lo
    span <= 0 && return (Float64[x_lo], [Dates.format(origin + Millisecond(round(Int, x_lo * 1000)), "yyyy-mm-dd HH:MM:SS")])
    step = _nice_step(span / (target - 1))
    first_tick = ceil(x_lo / step) * step
    positions = Float64[]
    v = first_tick
    while v <= x_hi + 1e-9
        push!(positions, v)
        v += step
    end
    isempty(positions) && (positions = [x_lo, x_hi])
    fmt = span < 2 * 86400 ? "yyyy-mm-dd\nHH:MM:SS" : "yyyy-mm-dd"
    labels = String[Dates.format(origin + Millisecond(round(Int, v * 1000)), fmt) for v in positions]
    return (positions, labels)
end

function _interactive_timearray()
    times = [DateTime(1970, 1, 1)]
    names = collect(LEMI424_MT_COMPONENTS)
    vals = fill(NaN, 1, length(names))
    metadata = Dict{Symbol, Any}(
        :site => "interactive",
        :instrument => "LEMI-424",
        :source_format => :lemi424,
        :sample_rate => 1.0,
        :start_time => first(times),
        :end_time => first(times),
        :n_samples => 1,
        :units => Dict(name => component_units(name) for name in names),
    )
    return TimeArray(times, vals, names, metadata)
end

function _load_data_file(path::AbstractString)
    ext = lowercase(splitext(path)[2])
    if ext == ".xyz"
        return _load_lemi_xyz(path), :lemi_xyz
    end
    return load_lemi424(path), :lemi424
end

function _load_lemi_xyz(path::AbstractString)
    raw = readlines(path)
    rows = [strip(l) for l in raw if !isempty(strip(l))]
    n = length(rows)
    n > 0 || error("XYZ file is empty: $path")
    times = Vector{DateTime}(undef, n)
    vals = Matrix{Float64}(undef, n, 5)
    for (i, line) in enumerate(rows)
        parts = split(line)
        length(parts) >= 7 || error("XYZ line $i has $(length(parts)) columns, expected 7")
        times[i] = DateTime(parts[1] * "T" * parts[2])
        for j in 1:5
            vals[i, j] = parse(Float64, parts[2 + j])
        end
    end
    names = [:Bx, :By, :Bz, :Ex, :Ey]
    fs = n > 1 ? 1.0 / (Dates.value(times[2] - times[1]) / 1000.0) : 1.0
    metadata = Dict{Symbol, Any}(
        :site => _site_from_path(path),
        :instrument => "LEMI (xyz)",
        :source_format => :lemi_xyz,
        :sample_rate => fs,
        :start_time => first(times),
        :end_time => last(times),
        :n_samples => n,
        :units => Dict(:Bx => "nT", :By => "nT", :Bz => "nT", :Ex => "mV/km", :Ey => "mV/km"),
    )
    return TimeArray(times, vals, names, metadata)
end

function _write_lemi_xyz(path::AbstractString, ta::TimeArray)
    times = _ensure_datetime(_ta_timestamps(ta))
    vals = _ta_values(ta)
    n = length(times)
    n_cols = min(5, size(vals, 2))
    open(path, "w") do io
        for i in 1:n
            t = times[i]
            row = ntuple(j -> j <= n_cols ? vals[i, j] : NaN, 5)
            @printf(io,
                "%04d-%02d-%02d %02d:%02d:%02d  %8.2f  %8.2f  %8.2f  %8.3f  %8.3f\n",
                year(t), month(t), day(t), hour(t), minute(t), second(t),
                row[1], row[2], row[3], row[4], row[5])
        end
    end
    return path
end

function _write_data_file(path::AbstractString, ta::TimeArray, source_format::Symbol)
    if source_format === :lemi_xyz
        return _write_lemi_xyz(path, ta)
    end
    return write_lemi424(path, ta)
end

function _summary_text(ta::TimeArray)
    metadata = _ta_meta(ta)
    site = metadata isa AbstractDict ? get(metadata, :site, "—") : "—"
    instrument = metadata isa AbstractDict ? get(metadata, :instrument, "—") : "—"
    n = length(_ta_timestamps(ta))
    fs = _sample_rate_from_timearray(ta)
    return "$(site)  ·  $(instrument)  ·  $(n) samples  ·  $(round(fs; digits = 4)) Hz"
end

function _refresh_line_obs!(app::TKApp)
    n_channels = length(app.line_clean)
    n = size(app.raw_values, 1)
    masked = app.mask.masked
    @assert length(masked) == n "Mask length $(length(masked)) != sample count $n"
    for j in 1:n_channels
        col = @view app.raw_values[:, j]
        clean = Vector{Float32}(undef, n)
        masked_y = Vector{Float32}(undef, n)
        @inbounds for i in 1:n
            v = Float32(col[i])
            if masked[i]
                clean[i] = NaN32
                masked_y[i] = v
            else
                clean[i] = v
                masked_y[i] = NaN32
            end
        end
        app.line_clean[j][] = clean
        app.line_masked[j][] = masked_y
    end
    return app
end

function _refresh_mask_overlay!(app::TKApp)
    lows = Float64[]
    highs = Float64[]
    for (a, b) in app.mask.intervals
        push!(lows, Dates.value(a - app.origin) / 1000.0)
        push!(highs, Dates.value(b - app.origin) / 1000.0)
    end
    app.mask_lows[] = lows
    app.mask_highs[] = highs
    _refresh_line_obs!(app)
    _refresh_status!(app)
    return app
end

function _format_dt(app::TKApp, secs::Float64)
    return Dates.format(app.origin + Millisecond(round(Int, secs * 1000)), "yyyy-mm-dd HH:MM:SS")
end

function _refresh_status!(app::TKApp)
    n_masked = masked_samples(app.mask)
    n_intervals = length(app.mask.intervals)
    if app.selection_visible[]
        lo, hi = app.selection[]
        sel_text = "Selection " * _format_dt(app, lo) * "  →  " * _format_dt(app, hi)
    else
        sel_text = "Left-drag on any panel to select a time range  ·  Right-drag = pan  ·  Scroll = zoom y"
    end
    app.status_label.text[] = "$(n_masked) masked samples in $(n_intervals) intervals    ·    $(sel_text)"
    return app
end

function _selection_to_datetimes(app::TKApp)
    lo, hi = app.selection[]
    return (
        app.origin + Millisecond(round(Int, lo * 1000)),
        app.origin + Millisecond(round(Int, hi * 1000)),
    )
end

function _apply_selection_mask!(app::TKApp, value::Bool)
    app.selection_visible[] || return
    t0, t1 = _selection_to_datetimes(app)
    if value
        mask_interval!(app.mask, t0, t1)
    else
        unmask_interval!(app.mask, t0, t1)
    end
    app.selection_visible[] = false
    _refresh_mask_overlay!(app)
end

function _clear_all_masks!(app::TKApp)
    clear_mask!(app.mask)
    _refresh_mask_overlay!(app)
end

function _autoscale_y!(app::TKApp, x_lo::Float64, x_hi::Float64)
    secs = app.time_seconds
    isempty(secs) && return
    idx_lo = searchsortedfirst(secs, x_lo)
    idx_hi = searchsortedlast(secs, x_hi)
    idx_lo = clamp(idx_lo, 1, length(secs))
    idx_hi = clamp(idx_hi, 1, length(secs))
    idx_lo > idx_hi && return
    for (j, ax) in enumerate(app.axes)
        col = @view app.raw_values[idx_lo:idx_hi, j]
        ymin = Inf
        ymax = -Inf
        @inbounds for v in col
            if isfinite(v)
                v < ymin && (ymin = v)
                v > ymax && (ymax = v)
            end
        end
        if isfinite(ymin) && isfinite(ymax)
            if ymin == ymax
                pad = max(abs(ymax) * 0.02, 1.0)
            else
                pad = (ymax - ymin) * 0.08
            end
            ylims!(ax, ymin - pad, ymax + pad)
        end
    end
end

function _update_x_window!(app::TKApp; force = false)
    isempty(app.axes) && return
    ws = app.window_seconds[]
    span = app.span_seconds
    visible = isfinite(ws) ? min(ws, max(span, 1.0)) : max(span, 1.0)
    x_lo = app.window_start[]
    if span > visible
        x_lo = clamp(x_lo, 0.0, span - visible)
    else
        x_lo = 0.0
    end
    x_hi = x_lo + visible
    ticks = _datetime_xticks_for_range(app.origin, x_lo, x_hi)
    for ax in app.axes
        ax.xticks[] = ticks
    end
    last_ax = last(app.axes)
    xlims!(last_ax, x_lo, x_hi)
    _autoscale_y!(app, x_lo, x_hi)
    return
end

mutable struct DragSelect
    app::TKApp
    dragging::Bool
    anchor::Float64
end

function Makie.process_interaction(s::DragSelect, event::Makie.MouseEvent, ax::Axis)
    et = event.type
    if et === Makie.MouseEventTypes.leftdragstart
        s.dragging = true
        s.anchor = event.data[1]
        s.app.selection[] = (event.data[1], event.data[1])
        s.app.selection_visible[] = true
        _refresh_status!(s.app)
        return Consume(true)
    elseif et === Makie.MouseEventTypes.leftdrag && s.dragging
        x = event.data[1]
        lo = min(s.anchor, x)
        hi = max(s.anchor, x)
        s.app.selection[] = (lo, hi)
        _refresh_status!(s.app)
        return Consume(true)
    elseif (et === Makie.MouseEventTypes.leftdragstop || et === Makie.MouseEventTypes.leftup) && s.dragging
        s.dragging = false
        _refresh_status!(s.app)
        return Consume(true)
    end
    return Consume(false)
end

function _build_axes!(app::TKApp, ta::TimeArray, existing::Vector{Axis})
    for ax in existing
        delete!(ax)
    end
    empty!(existing)
    empty!(app.line_clean)
    empty!(app.line_masked)

    names = _ta_colnames(ta)
    vals = _ta_values(ta)
    times = _ensure_datetime(_ta_timestamps(ta))
    metadata = _ta_meta(ta)
    units_map = metadata isa AbstractDict ? get(metadata, :units, Dict{Symbol, String}()) : Dict{Symbol, String}()

    origin = first(times)
    secs = _seconds_since(origin, times)
    span = isempty(secs) ? 0.0 : last(secs) - first(secs)

    app.origin = origin
    app.time_seconds = secs
    app.span_seconds = span
    app.raw_values = Matrix{Float64}(vals)
    app.line_x[] = secs

    n = length(names)
    axes = Axis[]
    for (i, name) in enumerate(names)
        is_last = i == n
        unit_str = get(units_map, _symbolize(name), component_units(_symbolize(name)))
        ax = Axis(app.plot_layout[i, 1];
            ylabel = "$(String(name))\n[$(unit_str)]",
            ylabelrotation = 0.0,
            ylabelpadding = 8.0,
            ylabelsize = 13,
            ylabelfont = :bold,
            yticklabelsize = 10,
            xticklabelsize = 10,
            xticklabelsvisible = is_last,
            xticksvisible = is_last,
            xgridvisible = false,
            ygridvisible = false,
            topspinevisible = false,
            rightspinevisible = false,
            yticks = LinearTicks(3),
            spinewidth = 0.9,
            bottomspinecolor = TK_FRAME,
            leftspinecolor = TK_FRAME,
            xtickcolor = TK_FRAME,
            ytickcolor = TK_FRAME,
            xticklabelcolor = TK_BLACK,
            yticklabelcolor = TK_BLACK,
            backgroundcolor = TK_PANEL_BG,
            xzoomlock = true,
            xpanlock = true,
            tellheight = false,
            tellwidth = false,
        )
        if !is_last
            hidexdecorations!(ax; ticks = true, ticklabels = true, grid = false)
        end

        clean_obs = Observable{Vector{Float32}}(Float32.(@view vals[:, i]))
        masked_obs = Observable{Vector{Float32}}(fill(NaN32, length(secs)))
        push!(app.line_clean, clean_obs)
        push!(app.line_masked, masked_obs)

        anchor = isempty(secs) ? 0.0 : first(secs)

        mask_lows_padded = lift(ls -> isempty(ls) ? Float64[anchor] : ls, app.mask_lows)
        mask_highs_padded = lift(hs -> isempty(hs) ? Float64[anchor] : hs, app.mask_highs)
        vspan!(ax, mask_lows_padded, mask_highs_padded; color = TK_MASK_FILL)

        sel_lows = lift((vis, sel) -> vis ? Float64[sel[1]] : Float64[anchor], app.selection_visible, app.selection)
        sel_highs = lift((vis, sel) -> vis ? Float64[sel[2]] : Float64[anchor], app.selection_visible, app.selection)
        vspan!(ax, sel_lows, sel_highs; color = TK_SEL_FILL, strokecolor = TK_SEL_EDGE, strokewidth = 0.8)

        col = _component_color(name)
        if length(secs) == 1
            scatter!(ax, secs, lift(y -> y, clean_obs); color = col, markersize = 4)
            scatter!(ax, secs, lift(y -> y, masked_obs); color = TK_MUTED, markersize = 4)
        else
            lines!(ax, app.line_x, clean_obs; color = col, linewidth = 1.4, joinstyle = :round)
            lines!(ax, app.line_x, masked_obs; color = TK_MUTED, linewidth = 1.4, joinstyle = :round)
        end

        deregister_interaction!(ax, :rectanglezoom)
        register_interaction!(ax, :tk_select, DragSelect(app, false, 0.0))

        push!(axes, ax)
    end
    if length(axes) > 1
        linkxaxes!(axes...)
        rowgap!(app.plot_layout, 8)
    end
    append!(existing, axes)
    return axes
end

function _refresh_slider_range!(app::TKApp)
    ws = app.window_seconds[]
    span = app.span_seconds
    visible = isfinite(ws) ? min(ws, max(span, 1.0)) : max(span, 1.0)
    max_start = max(0.0, span - visible)
    step = max(visible / 200.0, 1.0)
    if max_start <= 0
        rng = 0.0:1.0:0.0
    else
        rng = 0.0:step:max_start
    end
    app.slider.range[] = rng
    new_start = clamp(app.window_start[], 0.0, max_start)
    set_close_to!(app.slider, new_start)
    return
end

function _load_into_app!(app::TKApp, path::AbstractString)
    ta, fmt = _load_data_file(path)
    app.data = ta
    app.mask = TimekeeperMask(ta)
    app.source_format = fmt
    app.source_path = path
    app.selection_visible[] = false
    app.window_start[] = 0.0
    _build_axes!(app, ta, app.axes)
    _refresh_slider_range!(app)
    app.summary_label.text[] = _summary_text(ta)
    _refresh_mask_overlay!(app)
    _update_x_window!(app)
    return app
end

function _ext_for_format(fmt::Symbol)
    fmt === :lemi_xyz && return ".xyz"
    return ".txt"
end

function TKApp(
    ta::TimeArray;
    size = (1400, 1100),
    source_format::Symbol = :lemi424,
    source_path::AbstractString = "",
)
    GLMakie.activate!()
    fig = Figure(; size = size, backgroundcolor = :white, fontsize = 12,
        figure_padding = (14, 14, 8, 8))

    toolbar = GridLayout(fig[1, 1]; tellheight = true)
    summary_label = Label(toolbar[1, 1], _summary_text(ta);
        fontsize = 12, color = TK_GREY, halign = :left, tellwidth = false)

    actions = GridLayout(toolbar[1, 2]; tellheight = false, halign = :right)
    load_btn = Button(actions[1, 1]; label = "Load…", fontsize = 11)
    mask_btn = Button(actions[1, 2]; label = "Mask", fontsize = 11)
    unmask_btn = Button(actions[1, 3]; label = "Unmask", fontsize = 11)
    clear_btn = Button(actions[1, 4]; label = "Clear", fontsize = 11)
    write_btn = Button(actions[1, 5]; label = "Write", fontsize = 11)
    window_label = Label(actions[1, 6], "Window:"; fontsize = 11, color = TK_GREY)
    window_menu = Menu(actions[1, 7];
        options = WINDOW_OPTIONS, default = "1 hour",
        fontsize = 11, width = 110)
    colgap!(actions, 6)

    colsize!(toolbar, 1, Auto(true, 1.0))
    colsize!(toolbar, 2, Auto(false))

    plot_layout = GridLayout(fig[2, 1]; tellheight = false)

    slider_grid = GridLayout(fig[3, 1]; tellheight = true)
    Label(slider_grid[1, 1], "Scroll"; fontsize = 11, color = TK_GREY, halign = :right)
    slider = Slider(slider_grid[1, 2]; range = 0.0:1.0:0.0, startvalue = 0.0,
        color_active = TK_BLUE, color_active_dimmed = TK_SEL_FILL)
    colsize!(slider_grid, 1, Fixed(50))
    colsize!(slider_grid, 2, Auto(true, 1.0))
    colgap!(slider_grid, 8)

    status_label = Label(fig[4, 1], "";
        fontsize = 10.5, color = TK_GREY, halign = :left, tellwidth = false)

    rowsize!(fig.layout, 1, Fixed(38))
    rowsize!(fig.layout, 2, Auto(true, 1.0))
    rowsize!(fig.layout, 3, Fixed(28))
    rowsize!(fig.layout, 4, Fixed(20))
    rowgap!(fig.layout, 6)

    selection = Observable((0.0, 0.0))
    selection_visible = Observable(false)
    mask_lows = Observable(Float64[])
    mask_highs = Observable(Float64[])
    window_seconds_obs = Observable(3600.0)
    window_start_obs = Observable(0.0)
    line_x_obs = Observable(Float64[])

    app = TKApp(
        ta,
        TimekeeperMask(ta),
        fig,
        plot_layout,
        summary_label,
        status_label,
        Axis[],
        DateTime(1970),
        Float64[],
        0.0,
        Matrix{Float64}(undef, 0, 0),
        Observable{Vector{Float32}}[],
        Observable{Vector{Float32}}[],
        line_x_obs,
        window_seconds_obs,
        window_start_obs,
        slider,
        window_menu,
        selection,
        selection_visible,
        mask_lows,
        mask_highs,
        source_format,
        String(source_path),
    )

    _build_axes!(app, ta, app.axes)
    _refresh_mask_overlay!(app)
    _refresh_slider_range!(app)
    _update_x_window!(app)

    on(slider.value) do v
        app.window_start[] = Float64(v)
        _update_x_window!(app)
    end
    on(window_menu.selection) do secs
        secs === nothing && return
        app.window_seconds[] = Float64(secs)
        _refresh_slider_range!(app)
        _update_x_window!(app)
    end

    on(load_btn.clicks) do _
        path = ""
        try
            path = pick_file(; filterlist = "txt,dat,lem,xyz")
        catch err
            @warn "Could not open file picker" exception = err
            return
        end
        isempty(path) && return
        try
            _load_into_app!(app, path)
        catch err
            @warn "Could not load $path" exception = err
        end
    end
    on(mask_btn.clicks) do _
        _apply_selection_mask!(app, true)
    end
    on(unmask_btn.clicks) do _
        _apply_selection_mask!(app, false)
    end
    on(clear_btn.clicks) do _
        _clear_all_masks!(app)
    end
    on(write_btn.clicks) do _
        if isempty(app.source_path)
            @warn "Load a file first; nothing to write"
            app.status_label.text[] = "Load a file before writing"
            return
        end
        dir = dirname(app.source_path)
        stem, ext = splitext(basename(app.source_path))
        clean_path = joinpath(dir, stem * "_clean" * ext)
        mask_path = joinpath(dir, stem * "_mask.csv")
        try
            cleaned = cleaned_timearray(app)
            _write_data_file(clean_path, cleaned, app.source_format)
            write_mask(mask_path, app)
            @info "Wrote cleaned data and mask" clean_path mask_path
            app.status_label.text[] = "Wrote $(basename(clean_path))  and  $(basename(mask_path))"
        catch err
            @warn "Could not write outputs" exception = err
            app.status_label.text[] = "Write failed: $(sprint(showerror, err))"
        end
    end

    return app
end

function TKApp(path::AbstractString; kwargs...)
    ta, fmt = _load_data_file(path)
    return TKApp(ta; source_format = fmt, source_path = path, kwargs...)
end

function TKApp(; kwargs...)
    return TKApp(_interactive_timearray(); kwargs...)
end

function Base.display(app::TKApp)
    return display(app.figure)
end

function run_tkapp(app::TKApp)
    screen = display(app.figure)
    try
        wait(screen)
    catch
    end
    return app
end

run_tkapp(; kwargs...) = run_tkapp(TKApp(; kwargs...))
run_tkapp(path::AbstractString; kwargs...) = run_tkapp(TKApp(path; kwargs...))
run_tkapp(ta::TimeArray; kwargs...) = run_tkapp(TKApp(ta; kwargs...))

function cleaned_timearray(tk::TKApp; mode = :nan)
    return cleaned_timearray(tk.data, tk.mask; mode = mode)
end

function good_segments(tk::TKApp; min_samples = 1)
    return good_segments(tk.data, tk.mask; min_samples = min_samples)
end

function sample_weights(tk::TKApp; good = 1.0, bad = 0.0)
    return sample_weights(tk.mask; good = good, bad = bad)
end

function write_cleaned(path::AbstractString, tk::TKApp; mode = :nan, delimiter = ',')
    return write_cleaned(path, tk.data, tk.mask; mode = mode, delimiter = delimiter)
end

function write_mask(path::AbstractString, tk::TKApp; delimiter = ',')
    return write_mask(path, tk.mask; delimiter = delimiter)
end
