# Dashboard.jl - the survey dashboard, TKDash.
# Author: @pankajkmishra
#
# A GLMakie window for a scanned Survey (Survey.jl). It has a header to go
# through the sites, then a body with three columns:
# - the site map
# - the button that collapses the map
# - the Gantt charts of the runs
#
# At the start, the window shows the full survey: one chart of all the runs.
# When you select a site, all the views use that site:
# - The map shows the site as a magenta star. It shows its base sites (they
#   recorded with it, within the base distance) in blue, and its remote sites
#   (they recorded with it, farther) in amber. A marker is darker when its
#   site recorded with the site for a longer time.
# - The charts become two: the site above its base sites, and the site above
#   its remote sites. The time that each site recorded with the site has the
#   color of that site.
# The window shows only the best sites of each list. A right click drops a
# base or remote site from the lists of the site, or brings it back. Export
# writes the lists of each site as a text table that TKApp can read
#
# The survey is small (tens of sites). Thus, each change draws the map and the
# charts again from the start. The code does not update the plots

# Magenta for the site, steel blue for base sites and amber for remote sites.
# Persons with each type of color vision can tell them apart. Each bar and
# marker also has a black outline and a name
const DASH_SITE_COLOR = parse(RGBf, "#c2185b")
const DASH_BASE_RAMP = parse.(RGBf, ["#9dc3e6", "#7aa9d6", "#5a8fc4", "#3e74ae", "#2c5f94", "#1f4e79"])
const DASH_REMOTE_RAMP = parse.(RGBf, ["#f2d49b", "#ebc47f", "#e2b263", "#d49e4a", "#c68b34", "#b7791f"])
const DASH_DOT = RGBf(0.29, 0.33, 0.39)
const DASH_IDLE = RGBf(0.80, 0.82, 0.85)
const DASH_RUN = RGBf(0.84, 0.85, 0.87)
const DASH_RUN_ALL = parse(RGBf, "#3c8d93")              # the overview: not a color of a role
const DASH_RUN_FAINT = RGBf(0.93, 0.94, 0.95)
const DASH_PANEL = parse(RGBf, "#f7f8fa")
const DASH_GUIDE = RGBAf(0.0, 0.0, 0.0, 0.07)
const DASH_OUTLINE = 0.75
const DASH_BAR_HALF = 0.175                             # bars fill 35 % of a row
const DASH_EDGE_PX = 22
const DASH_NAVY = parse(RGBf, "#1b2f5b")                  # site names in the menu and the charts
const DASH_COMMON = RGBAf(0.235, 0.553, 0.576, 0.14)  # light teal: the common window
const DASH_COMMON_TEXT = parse(RGBf, "#2b6f74")
const DASH_STATUS_OK = parse(RGBf, "#1098ad")             # cyan: the action was successful
const DASH_HINT = "Click a site on the map or a chart row to see its base and remote sites · " *
                  "right-click one to drop it · Export saves the table."
const DASH_STATUS_ERROR = parse(RGBf, "#e03131")
const DASH_SHOWN = 5

# The shade of a ramp for an overlap `h`, if the overlaps go up to `hmax`. The
# longest overlap is the darkest
function _ramp(ramp, h, hmax)
    f = hmax > 0 ? clamp(h / hmax, 0, 1) : 1.0
    return ramp[1 + round(Int, f * (length(ramp) - 1))]
end

"""
    TKDash

The survey dashboard. It holds:
- a [`Survey`](@ref);
- the site in focus (`0` for the overview);
- the rate and the limits on the screen;
- the base and remote sites that you dropped by hand;
- the GLMakie figure.

To make one, use `TKDash(root)` and `display` it. Or use [`run_tkdash`](@ref),
which opens it and waits.

After you close the window, your choices stay. `reference_plan(dash)` returns
the plan as the window shows it, and `write_reference_plan(path, dash)`
writes it.
"""
mutable struct TKDash
    figure::Figure
    header::GridLayout
    survey::Survey
    focus::Int
    rate::Union{Nothing, Symbol, Float64}
    base_km::Float64
    remote_km::Float64
    min_overlap_hours::Float64
    shown::Int
    exclude::Dict{String, Set{String}}
    refs::Any
    map_open::Bool
    map_axis::Union{Nothing, Axis}
    map_legend::Any
    map_grid::GridLayout
    zoom_button::Any
    body::GridLayout
    chart_grid::GridLayout
    charts::Vector{Axis}
    all_charts::Vector{Axis}
    chart_rows::Vector{Vector{Int}}
    edge::Button
    site_menu::Menu
    rate_menu::Menu
    status::Label
    common::Vector{Tuple{DateTime, DateTime}}
    t0::DateTime
    tips::Vector{Any}
    updating::Bool
end

function reference_plan(d::TKDash)
    return reference_plan(d.survey; rate = d.rate, base_km = d.base_km, remote_km = d.remote_km,
                          min_overlap_hours = d.min_overlap_hours, exclude = d.exclude)
end

write_reference_plan(path::AbstractString, d::TKDash) = write_reference_plan(path, reference_plan(d))

#---------- formatting -----

# 42732 s is "11.9h"
_hours(seconds::Real) = @sprintf("%.1fh", seconds / 3600)

# The label at the end of the bar of a site: "[own]" for the site in focus and
# for the overview, "[overlap/own]" for a base or remote site
_bar_label(own_seconds) = "[" * _hours(own_seconds) * "]"
_bar_label(overlap_seconds, own_seconds) = "[" * _hours(overlap_seconds) * "/" * _hours(own_seconds) * "]"

_number_text(v::Real) = isinteger(v) ? string(Int(v)) : string(v)

# Each sampling rate in the survey, then all the rates together
_rate_menu_options(s::Survey) =
    vcat([(_fs_label(fs), fs) for fs in survey_rates(s)], [("All rates", :all)])

# The rate of a survey at the start:
# 1. 128 Hz, the band that most surveys share, if a site recorded it.
# 2. If not, all rates for a survey with different instruments, because their
#    rates are usually different.
# 3. If not, the rate with the most hours recorded together, for all pairs of
#    sites.
# 4. If there is no overlap, the rate with the longest recording
function _default_rate(s::Survey)
    rates = survey_rates(s)
    any(fs -> _same_rate(fs, 128), rates) && return 128.0
    length(unique(x.format for x in s.sites)) > 1 && return :all
    isempty(rates) && return :all
    shared = [(m = overlap_matrix(s; rate = fs); sum(m) - sum(m[i, i] for i in axes(m, 1))) for fs in rates]
    maximum(shared) > 0 && return rates[argmax(shared)]
    return rates[argmax([sum(recording_seconds(x; rate = fs) for x in s.sites) for fs in rates])]
end

# Ticks at clock times, for an axis in hours from the midnight `t0`. The steps
# go from minutes to years
function _time_ticks(t0::DateTime)
    day = 24.0
    steps = [1 / 60, 5 / 60, 15 / 60, 0.5, 1, 2, 3, 6, 12, day, 2day, 4day, 7day, 14day,
             30day, 91day, 182day, 365day, 730day, 1826day]
    return function (vmin, vmax)
        span = vmax - vmin
        step = steps[something(findfirst(s -> span / s <= 7, steps), length(steps))]
        if step >= 30day
            # whole months or years, from the month of the first tick
            months = step >= 365day ? 12 * round(Int, step / 365day) : round(Int, step / 30day)
            t = floor(DateTime(Date(t0 + Millisecond(round(Int, vmin * 3.6e6)))), Month)
            ts = DateTime[]
            while Dates.value(t - t0) / 3.6e6 <= vmax
                Dates.value(t - t0) / 3.6e6 >= vmin && push!(ts, t)
                t += Month(months)
            end
            return ([Dates.value(t - t0) / 3.6e6 for t in ts],
                    [Dates.format(t, months >= 12 ? "yyyy" : "u yyyy") for t in ts])
        end
        vals = collect((ceil(vmin / step) * step):step:vmax)
        labels = String[]
        last_day = nothing
        for v in vals
            t = t0 + Millisecond(round(Int, v * 3.6e6))
            if step >= day
                push!(labels, Dates.format(t, "d u"))
            elseif Date(t) != last_day
                push!(labels, Dates.format(t, "HH:MM") * "\n" * Dates.format(t, "d u yyyy"))
            else
                push!(labels, Dates.format(t, "HH:MM"))
            end
            last_day = Date(t)
        end
        return (vals, labels)
    end
end

#---------- building -----

"""
    TKDash(root; size = (1600, 900), kwargs...) -> TKDash
    TKDash(survey::Survey; size = (1600, 900), base_km = 5.0, remote_km = 20.0,
           min_overlap_hours = 1.0, rate = :auto, shown = 5, show_map = true) -> TKDash

Make the dashboard for a survey directory (the function scans it with
[`scan_survey`](@ref)) or for a [`Survey`](@ref) that you scanned. The
function does not open a window.

- Sites that recorded with the site in focus for `min_overlap_hours` or more
  are its base sites within `base_km`, and its remote sites at `remote_km` or
  more.
- The window shows the best `shown` sites of each list.
- `rate` has the same meaning as for [`overlap_intervals`](@ref). With
  `:auto`, the window opens at 128 Hz if a site recorded that rate. If not, a
  survey with different instruments opens at `:all`, and other surveys open
  at the rate with the most overlap.
- With `show_map = false`, the map starts in the collapsed state.
"""
TKDash(root::AbstractString; kwargs...) = TKDash(scan_survey(root); kwargs...)

function TKDash(survey::Survey; size = (1600, 900), base_km::Real = 5.0, remote_km::Real = 20.0,
                min_overlap_hours::Real = 1.0, rate = :auto, shown::Integer = DASH_SHOWN,
                show_map::Bool = true)
    GLMakie.activate!(title = "TKDash")
    fig = Figure(; size = size, figure_padding = (16, 16, 10, 10))

    # header: go through the sites, set the rules, and do the file actions
    header = GridLayout(fig[1, 1]; tellwidth = false)
    b_first = Button(header[1, 1]; label = "|<")
    b_prev = Button(header[1, 2]; label = "< Prev")
    site_menu = _dash_menu(header[1, 3]; options = [("Overview", 0)], width = 140, textcolor = DASH_NAVY)
    b_next = Button(header[1, 4]; label = "Next >")
    b_last = Button(header[1, 5]; label = ">|")
    b_overview = Button(header[1, 6]; label = "Overview")
    b_restore = Button(header[1, 7]; label = "Restore")
    Label(header[1, 8], "Rate"; padding = (14, 0, 0, 0))
    rate_menu = _dash_menu(header[1, 9]; options = [("All rates", :all)], width = 110)
    Label(header[1, 10], "Base ≤"; padding = (10, 0, 0, 0))
    base_box = _dash_number_box(header[1, 11], base_km)
    Label(header[1, 12], "km   Remote ≥")
    remote_box = _dash_number_box(header[1, 13], remote_km)
    Label(header[1, 14], "km   Overlap ≥")
    overlap_box = _dash_number_box(header[1, 15], min_overlap_hours)
    Label(header[1, 16], "h")
    Box(header[1, 17]; visible = false)
    b_open = Button(header[1, 18]; label = "Open…")
    b_export = Button(header[1, 19]; label = "Export…")
    colsize!(header, 17, Auto(true, 1.0))
    colgap!(header, 6)

    # body: map | edge button | charts. The code makes the map when it opens.
    # It makes the charts one time and shows or hides them for each view
    body = GridLayout(fig[2, 1]; tellwidth = false)
    edge = Button(body[1, 2]; label = "‹", width = DASH_EDGE_PX, height = Relative(1),
                  tellheight = false, cornerradius = 2, buttoncolor = RGBf(0.93, 0.93, 0.93),
                  fontsize = 26, padding = (0, 0, 0, 0))
    map_grid = GridLayout(body[1, 1]; tellwidth = false)
    chart_grid = GridLayout(body[1, 3])
    colgap!(body, 8)

    status = Label(fig[3, 1], ""; fontsize = 12, halign = :left, tellwidth = false)
    colsize!(fig.layout, 1, Relative(1))
    rowsize!(fig.layout, 2, Auto(true, 1.0))
    rowgap!(fig.layout, 8)

    d = TKDash(fig, header, survey, 0, nothing, Float64(base_km), Float64(remote_km), Float64(min_overlap_hours), Int(shown),
               Dict{String, Set{String}}(), nothing, show_map, nothing, nothing, map_grid, nothing,
               body, chart_grid,
               Axis[], _make_charts!(chart_grid), Vector{Int}[], edge, site_menu, rate_menu, status,
               Tuple{DateTime, DateTime}[], DateTime(2000), Any[], false)
    _set_survey!(d, survey; rate)

    n() = length(d.survey.sites)
    on(_ -> n() > 0 && _focus!(d, 1), b_first.clicks)
    on(_ -> n() > 0 && _focus!(d, d.focus <= 1 ? n() : d.focus - 1), b_prev.clicks)
    on(_ -> n() > 0 && _focus!(d, d.focus % n() + 1), b_next.clicks)
    on(_ -> n() > 0 && _focus!(d, n()), b_last.clicks)
    on(_ -> _focus!(d, 0), b_overview.clicks)
    on(_ -> _restore!(d), b_restore.clicks)
    on(site_menu.selection) do i
        d.updating || i === nothing || i == d.focus || _focus!(d, i)
    end
    on(rate_menu.selection) do fs
        d.updating && return
        d.rate = fs
        _refresh!(d)
    end
    on(base_box.stored_string) do s
        v = tryparse(Float64, s)
        v === nothing || (d.base_km = v; _refresh!(d))
    end
    on(remote_box.stored_string) do s
        v = tryparse(Float64, s)
        v === nothing || (d.remote_km = v; _refresh!(d))
    end
    on(overlap_box.stored_string) do s
        v = tryparse(Float64, s)
        v === nothing || (d.min_overlap_hours = v; _refresh!(d))
    end
    on(edge.clicks) do _
        d.map_open = !d.map_open
        _layout_map!(d)
        _refresh!(d)
    end
    on(b_open.clicks) do _
        dir = try
            pick_folder()
        catch err
            _dash_status!(d, "Could not open a folder dialog: $(sprint(showerror, err))"; error = true)
            ""
        end
        isempty(dir) || _rescan!(d, dir)
    end
    on(_ -> _export_plan!(d), b_export.clicks)
    _install_clicks!(d)
    _install_hover!(d)

    _layout_map!(d)
    _refresh!(d)
    _dash_status!(d, "")                                # the instruction line
    return d
end

# A plain white menu: the cells have no tint, and the list is opaque. Thus, an
# open menu is easy to read above the charts and axes. A long list scrolls,
# and text that you type filters the list
function _dash_menu(pos; textcolor = :black, kwargs...)
    return Menu(pos;
        cell_color_inactive_even = :white, cell_color_inactive_odd = :white,
        cell_color_hover = RGBf(0.93, 0.95, 0.98), cell_color_active = RGBf(0.87, 0.91, 0.97),
        selection_cell_color_inactive = :white,
        textcolor = textcolor, textcolor_active = textcolor, textcolor_hover = textcolor,
        dropdown_arrow_color = textcolor, kwargs...)
end

function _dash_number_box(pos, value::Real)
    return Textbox(pos; stored_string = _number_text(value), width = 52,
        validator = s -> (v = tryparse(Float64, s); v !== nothing && v >= 0), halign = :left)
end

function Base.display(d::TKDash)
    display(d.figure)
    return d
end

# Open or close the map column. The edge button stays. The map keeps the shape
# of the survey at true scale, and it becomes larger and smaller with the
# window:
# - Scroll to zoom.
# - Draw a box with the left button to zoom to the box.
# - Drag with the right button to pan.
# - Click a site to select it.
# The zoom stays when you change the site. The Reset Zoom button below the map
# shows the full survey again
function _layout_map!(d::TKDash)
    d.map_legend === nothing || (delete!(d.map_legend); d.map_legend = nothing)
    d.map_axis === nothing || (delete!(d.map_axis); d.map_axis = nothing)
    d.zoom_button === nothing || (delete!(d.zoom_button); d.zoom_button = nothing)
    if d.map_open
        # the map is directly above the button. Thus, the button is part of the map
        d.map_axis = Axis(d.map_grid[1, 1]; tellheight = false, tellwidth = false, valign = :bottom,
            backgroundcolor = DASH_PANEL, xgridvisible = false, ygridvisible = false,
            xlabelfont = :regular, ylabelfont = :regular,
            xlabel = "Longitude (°)", ylabel = "Latitude (°)",
            xticks = WilkinsonTicks(3), yticks = WilkinsonTicks(4))
        d.zoom_button = Button(d.map_grid[2, 1]; label = "Reset Zoom", halign = :right, tellwidth = false)
        on(_ -> _reset_zoom!(d), d.zoom_button.clicks)
        rowgap!(d.map_grid, 6)
        # the legend never changes. The code makes it with the map and shows it
        # when a site is in focus
        d.map_legend = axislegend(d.map_axis,
            [MarkerElement(; marker = :star5, color = DASH_SITE_COLOR, strokecolor = :black, strokewidth = 0.8, markersize = 16),
             MarkerElement(; marker = :circle, color = DASH_BASE_RAMP[4], strokecolor = :black, strokewidth = 0.8, markersize = 11),
             MarkerElement(; marker = :circle, color = DASH_REMOTE_RAMP[4], strokecolor = :black, strokewidth = 0.8, markersize = 11),
             MarkerElement(; marker = :circle, color = DASH_IDLE, strokecolor = RGBf(0.6, 0.6, 0.62), strokewidth = 0.6, markersize = 7)],
            ["site", "base", "remote", "other"];
            position = :lt, labelsize = 11, rowgap = 0, padding = (6, 6, 4, 4), patchsize = (14, 14),
            framecolor = RGBAf(0, 0, 0, 0.15), backgroundcolor = RGBAf(1, 1, 1, 0.9))
        _map_home!(d)
    end
    # the map and the charts share the width 1 : 2, for all contents
    colsize!(d.body, 1, d.map_open ? Auto(false, 1.0) : Fixed(0))
    colsize!(d.body, 3, Auto(false, 2.0))
    d.edge.label[] = d.map_open ? "‹" : "›"
    return d
end

# The full survey at true scale. Reset Zoom goes back to this frame
function _map_home!(d::TKDash)
    d.map_axis === nothing && return d
    lon, lat, aspect = _map_frame(d.survey.sites)
    d.map_axis.aspect = AxisAspect(aspect)
    d.map_axis.limits = (lon, lat)
    reset_limits!(d.map_axis)
    return d
end

function _reset_zoom!(d::TKDash)
    _map_home!(d)
    _refresh!(d)                                        # the redraw resets the charts
    return d
end

# The three charts, made one time: the overview, and the base chart above the
# remote chart of a site. A change of view shows some charts and hides the
# others. It does not make new axes, because new axes made the change to a
# site slow
function _make_charts!(grid::GridLayout)
    charts = map(1:3) do k
        ax = Axis(grid[k, 1]; yreversed = true, ygridvisible = false, yticksvisible = false,
                  backgroundcolor = DASH_PANEL, xgridcolor = RGBAf(0, 0, 0, 0.06),
                  titlealign = :left, titlefont = :regular,
                  xrectzoom = false, yrectzoom = false, yzoomlock = true, ypanlock = true)
        deactivate_interaction!(ax, :dragpan)          # the right button drops sites
        ax
    end
    linkxaxes!(charts[2], charts[3])
    return charts
end

# One chart for the overview (n = 1). For a site, its base chart above its
# remote chart (n = 2). No chart for an empty survey. Hidden charts use no
# space
function _layout_charts!(d::TKDash, n::Integer)
    on = n == 1 ? (true, false, false) : n == 2 ? (false, true, true) : (false, false, false)
    for (k, (ax, show)) in enumerate(zip(d.all_charts, on))
        ax.blockscene.visible[] = show
        ax.scene.visible[] = show
        rowsize!(d.chart_grid, k, show ? Auto(1.0) : Fixed(0))
    end
    rowgap!(d.chart_grid, 1, 0)
    rowgap!(d.chart_grid, 2, n == 2 ? 14 : 0)
    d.charts = [ax for (ax, show) in zip(d.all_charts, on) if show]
    return d
end

#---------- state changes -----

function _set_survey!(d::TKDash, survey::Survey; rate = :auto)
    d.survey = survey
    d.exclude = Dict{String, Set{String}}()
    d.focus = 0
    d.rate = rate === :auto ? _default_rate(survey) : rate
    d.updating = true
    try
        d.site_menu.options[] = vcat([("Overview", 0)], [(s.name, i) for (i, s) in enumerate(survey.sites)])
        d.site_menu.i_selected[] = 1
        opts = _rate_menu_options(survey)
        d.rate_menu.options[] = opts
        d.rate_menu.i_selected[] = something(findfirst(o -> o[2] == d.rate, opts), 1)
    finally
        d.updating = false
    end
    return d
end

function _focus!(d::TKDash, i::Integer)
    (0 <= i <= length(d.survey.sites)) || return d
    d.focus = i
    _dash_status!(d, "")                                # the feedback was for the previous site
    d.updating = true
    try
        d.site_menu.i_selected[] = i + 1
    finally
        d.updating = false
    end
    _refresh!(d)
    return d
end

# Bring back the sites dropped for the site in focus. In the overview, bring
# back the dropped sites of all the sites
function _restore!(d::TKDash)
    if d.focus == 0
        n = sum(length, values(d.exclude); init = 0)
        empty!(d.exclude)
        what = "every site"
    else
        name = d.survey.sites[d.focus].name
        n = length(pop!(d.exclude, name, Set{String}()))
        what = name
    end
    _refresh!(d)
    _dash_status!(d, n == 0 ? "Nothing dropped for $what." :
                     "Restored $n dropped site$(n == 1 ? "" : "s") for $what.")
    return d
end

# Drop a base or remote site from the lists of the site, or bring it back
function _toggle!(d::TKDash, i::Integer)
    (d.focus == 0 || i == d.focus || d.refs === nothing) && return d
    name = d.survey.sites[i].name
    any(c -> c.site == name, vcat(d.refs.base, d.refs.remote)) || return d
    ex = get!(d.exclude, d.survey.sites[d.focus].name, Set{String}())
    name in ex ? delete!(ex, name) : push!(ex, name)
    _refresh!(d)
    _dash_status!(d, (name in ex ? "Dropped $name from" : "Restored $name to") *
                     " the lists of $(d.survey.sites[d.focus].name).")
    return d
end

function _rescan!(d::TKDash, dir::AbstractString)
    _dash_status!(d, "Scanning $dir …")
    try
        t = @elapsed survey = scan_survey(dir)
        _set_survey!(d, survey)
        _map_home!(d)
        _refresh!(d)
        _dash_status!(d, @sprintf("Scanned %d sites in %.1f s.", length(survey.sites), t))
    catch err
        _dash_status!(d, "Could not scan $dir: $(sprint(showerror, err))"; error = true)
    end
    return d
end

function _export_plan!(d::TKDash)
    isempty(d.survey.sites) && return _dash_status!(d, "No sites to export."; error = true)
    default = joinpath(d.survey.root, "reference_plan.txt")
    path = try
        save_file(default; filterlist = "txt")
    catch
        default
    end
    isempty(path) && return d
    endswith(lowercase(path), ".txt") || (path *= ".txt")
    try
        write_reference_plan(path, d)
        _dash_status!(d, "Wrote the base and remote sites of $(length(d.survey.sites)) sites to $path")
    catch err
        _dash_status!(d, "Could not write $path: $(sprint(showerror, err))"; error = true)
    end
    return d
end

# The line below the charts. It tells what to do, in grey, until there is a
# new message. A message is cyan when the action was successful and red when
# it failed. An empty message shows the instruction again
function _dash_status!(d::TKDash, text::AbstractString; error::Bool = false)
    d.status.color[] = isempty(text) ? TK_GREY : error ? DASH_STATUS_ERROR : DASH_STATUS_OK
    d.status.text[] = isempty(text) ? DASH_HINT : text
    return d
end


#---------- drawing -----

# The base and remote sites that the charts compare with the site: the best
# `shown` sites of each list that are not dropped. A dropped site leaves the
# comparison, and the next best site takes its position. The dropped site
# stays on the map as a hollow marker. A right click there brings it back
_shown(d::TKDash, list) = first([c for c in list if !c.excluded], d.shown)

function _refresh!(d::TKDash)
    sites = d.survey.sites
    d.map_legend === nothing || (d.map_legend.blockscene.visible[] = d.focus != 0)
    # a redraw keeps the zoom of the map. Only Reset Zoom changes it
    zoom = d.map_axis === nothing ? nothing : d.map_axis.targetlimits[]
    d.map_axis === nothing || empty!(d.map_axis)
    if isempty(sites)
        _layout_charts!(d, 0)
        d.refs = nothing
        d.chart_rows = Vector{Int}[]
        _dash_status!(d, "No sites found under $(d.survey.root); open a survey directory holding Metronix, LEMI-424 or GEOMAG sites.")
        return d
    end
    d.refs = d.focus == 0 ? nothing :
             site_references(d.survey, d.focus; rate = d.rate, base_km = d.base_km, remote_km = d.remote_km,
                             min_overlap_hours = d.min_overlap_hours,
                             exclude = get(d.exclude, sites[d.focus].name, ()))
    # the part that the site shares with all the base and remote sites that it
    # keeps, also the sites that are not on the screen
    d.common = d.focus == 0 ? Tuple{DateTime, DateTime}[] :
               common_window(d.survey, d.focus,
                             [c.site for c in vcat(d.refs.base, d.refs.remote) if !c.excluded]; rate = d.rate)
    _draw_charts!(d)
    if d.map_axis !== nothing
        _draw_map!(d)
        d.map_axis.targetlimits[] = zoom
    end
    return d
end

# The role of each site on the screen for the site in focus, by site index:
# (:site | :base | :remote, color, its site_references entry or nothing)
function _site_roles(d::TKDash)
    roles = Dict{Int, Tuple{Symbol, RGBf, Any}}()
    d.focus == 0 && return roles
    roles[d.focus] = (:site, DASH_SITE_COLOR, nothing)
    index = Dict(s.name => i for (i, s) in enumerate(d.survey.sites))
    for (role, list, ramp) in ((:base, d.refs.base, DASH_BASE_RAMP), (:remote, d.refs.remote, DASH_REMOTE_RAMP))
        shown = _shown(d, list)
        hmax = maximum((c.overlap_hours for c in shown); init = 0.0)
        foreach(c -> roles[index[c.site]] = (role, _ramp(ramp, c.overlap_hours, hmax), c), shown)
        # dropped sites keep a role for the hollow marker of the map and the
        # tooltip
        foreach(c -> c.excluded && (roles[index[c.site]] = (role, DASH_IDLE, c)), list)
    end
    return roles
end

function _draw_charts!(d::TKDash)
    sites = d.survey.sites
    n = length(sites)
    focus = d.focus
    roles = _site_roles(d)
    index = Dict(s.name => i for (i, s) in enumerate(sites))
    at_rate(r) = _at_rate(r, d.rate)

    # rows of each chart: all the sites for the overview. For a site: the site
    # above its base sites on the screen, and the site above its remote sites
    # on the screen
    groups = focus == 0 ? [(nothing, collect(1:n))] :
             [(role, vcat(focus, [index[c.site] for c in _shown(d, list)]))
              for (role, list) in ((:base, d.refs.base), (:remote, d.refs.remote))]
    length(d.charts) == length(groups) || _layout_charts!(d, length(groups))
    d.chart_rows = [rows for (_, rows) in groups]

    # one time frame for the two charts: the rows on the screen at the selected
    # rate
    t0 = d.t0 = DateTime(Date(minimum(r.start for s in sites for r in s.runs)))
    hrs(t) = Dates.value(t - t0) / 3.6e6
    near = [r for (_, rows) in groups for i in rows for r in sites[i].runs if at_rate(r)]
    isempty(near) && (near = [r for s in sites for r in s.runs])
    x0, x1 = hrs(minimum(r.start for r in near)), hrs(maximum(r.stop for r in near))
    span = max(x1 - x0, 1.0e-3)
    minw = 0.002 * span                                 # thus, a short run is visible
    bar(a, b, y) = Rect2f(hrs(a), y - DASH_BAR_HALF, max(hrs(b) - hrs(a), minw), 2DASH_BAR_HALF)

    for (ax, (role, rows)) in zip(d.charts, groups)
        empty!(ax)
        m = length(rows)
        ax.title = role === nothing ? "" : role === :base ? "base" : "remote"
        hlines!(ax, 1:m; color = DASH_GUIDE, linewidth = 1)    # row guides
        rows_on_screen = max(m, (focus == 0 ? 2 : 1) * (d.shown + 1))
        if !isempty(d.common)
            for (a, b) in d.common
                poly!(ax, Rect2f(hrs(a), 0.4, hrs(b) - hrs(a), rows_on_screen + 0.2);
                      color = DASH_COMMON, strokewidth = 0)
            end
            text!(ax, hrs(first(d.common)[1]), 0.42;
                  text = "common " * _bar_label(3600 * _window_hours(d.common)),
                  align = (:left, :top), fontsize = 11, color = DASH_COMMON_TEXT, offset = (4, 0))
        end
        faint, full, colors = Rect2f[], Rect2f[], RGBf[]
        for (y, i) in enumerate(rows), r in sites[i].runs
            if at_rate(r)
                push!(full, bar(r.start, r.stop, y))
                push!(colors, focus == 0 ? DASH_RUN_ALL : i == focus ? DASH_SITE_COLOR : DASH_RUN)
            else
                push!(faint, bar(r.start, r.stop, y))
            end
        end
        isempty(faint) || poly!(ax, faint; color = DASH_RUN_FAINT, strokewidth = DASH_OUTLINE,
                                strokecolor = RGBAf(0, 0, 0, 0.25))
        isempty(full) || poly!(ax, full; color = colors, strokewidth = DASH_OUTLINE, strokecolor = :black)

        # the time that each base or remote site recorded with the site, in its
        # shade
        for (y, i) in enumerate(rows)
            r = get(roles, i, nothing)
            (r === nothing || r[1] === :site) && continue
            poly!(ax, [bar(a, b, y) for (a, b) in r[3].windows]; color = r[2],
                  strokewidth = DASH_OUTLINE, strokecolor = :black)
        end

        # durations immediately after the last bar of each row: the recording of
        # the site, after the overlap for a base or remote site
        for (y, i) in enumerate(rows)
            isempty(sites[i].runs) && continue
            r = get(roles, i, nothing)
            own = recording_seconds(sites[i]; rate = d.rate)
            txt = r === nothing || r[1] === :site ? _bar_label(own) :
                  _bar_label(3600r[3].overlap_hours, own)
            drawn = [r for r in sites[i].runs if at_rate(r)]
            isempty(drawn) && (drawn = sites[i].runs)
            xend = maximum(hrs(r.stop) for r in drawn)
            text!(ax, xend, y; text = txt, align = (:left, :center), fontsize = 12, offset = (6, 0))
        end
        m == 1 && focus != 0 &&
            text!(ax, x0 + 0.5span, 1.75; text = "no $(role) site recorded with $(sites[focus].name)",
                  align = (:center, :center), fontsize = 12)

        # the site in focus is navy, as in the site menu
        ax.yticks = (1:m, [rich(_row_name(sites[i]); color = i == focus ? DASH_NAVY : :black)
                           for i in rows])
        ax.xticks = _time_ticks(t0)
        xlims!(ax, x0 - 0.02span, x1 + 0.12span)            # room for the durations
        # space for the site and its best `shown` sites in each chart. Thus, a
        # bar has the same thickness with one site or six sites on the screen.
        # The one chart of the overview is two times as tall. Thus, it keeps
        # space for two times as many rows
        ylims!(ax, rows_on_screen + 0.6, 0.4)
    end
    return d
end

# The longitude and latitude limits around each site with a position, with 8 %
# on each side, and the aspect of that frame at true scale
function _map_frame(sites)
    located = [s for s in sites if isfinite(s.latitude)]
    isempty(located) && return ((-1.0, 1.0), (-1.0, 1.0), 1.0)
    pad(v) = (lo = minimum(v); hi = maximum(v); m = max(0.08 * (hi - lo), 0.01); (lo - m, hi + m))
    lon = pad([s.longitude for s in located])
    lat = pad([s.latitude for s in located])
    lat0 = sum(s.latitude for s in located) / length(located)
    return (lon, lat, (lon[2] - lon[1]) * cosd(lat0) / (lat[2] - lat[1]))
end

function _draw_map!(d::TKDash)
    ax = d.map_axis
    sites = d.survey.sites
    located = [i for (i, s) in enumerate(sites) if isfinite(s.latitude)]
    pt(i) = Point2f(sites[i].longitude, sites[i].latitude)
    roles = _site_roles(d)

    if d.focus == 0
        isempty(located) || scatter!(ax, pt.(located); color = DASH_DOT, strokecolor = :black,
                                     strokewidth = 0.8, markersize = 9)
    else
        idle = [i for i in located if !haskey(roles, i)]
        isempty(idle) || scatter!(ax, pt.(idle); color = DASH_IDLE, strokecolor = RGBf(0.6, 0.6, 0.62),
                                  strokewidth = 0.6, markersize = 7)
        # remote below base, the shortest overlap below the longest, the site on
        # top
        linked = sort([(i, r) for (i, r) in roles if r[1] !== :site && i in located];
                      by = ((i, r),) -> (r[1] === :base, r[3].excluded ? -1.0 : r[3].overlap_hours))
        for (i, (_, color, c)) in linked
            if c.excluded
                scatter!(ax, [pt(i)]; color = :white, strokecolor = TK_MUTED, strokewidth = 1.2,
                         markersize = 11)
            else
                scatter!(ax, [pt(i)]; color = color, strokecolor = :black, strokewidth = 0.8,
                         markersize = 13)
            end
        end
        d.focus in located && scatter!(ax, [pt(d.focus)]; color = DASH_SITE_COLOR, marker = :star5,
                                       strokecolor = :black, strokewidth = 0.8, markersize = 22)
    end

    return ax
end

#---------- clicks -----

function _install_clicks!(d::TKDash)
    on(events(d.figure).mousebutton; priority = 1) do ev
        ev.action == Mouse.press || return Consume(false)
        ev.button in (Mouse.left, Mouse.right) || return Consume(false)
        isempty(d.survey.sites) && return Consume(false)
        i = nothing
        if d.map_axis !== nothing && is_mouseinside(d.map_axis.scene)
            i = _map_hit(d, mouseposition(d.map_axis.scene))
        else
            for (ax, rows) in zip(d.charts, d.chart_rows)
                is_mouseinside(ax.scene) || continue
                r = round(Int, mouseposition(ax.scene)[2])
                1 <= r <= length(rows) && (i = rows[r])
            end
        end
        i === nothing && return Consume(false)
        ev.button == Mouse.left ? (i == d.focus || _focus!(d, i)) : _toggle!(d, i)
        return Consume(true)
    end
    return d
end

#---------- hover -----

# One tooltip for the full window. The code draws it above all the panels, in
# pixels. Thus, no axis clips it. It follows the cursor and moves away from
# the edges of the window
function _install_hover!(d::TKDash)
    pos, txt, vis = Observable(Point2f(0, 0)), Observable(" "), Observable(false)
    placement = Observable(:above)
    tip = tooltip!(d.figure.scene, pos, txt; visible = vis, placement = placement, fontsize = 12,
                   textpadding = (8, 8, 5, 5), outline_linewidth = 0.75, justification = :left,
                   backgroundcolor = RGBAf(1, 1, 1, 0.96), overdraw = true, depth_shift = -1.0f0)
    translate!(tip, 0, 0, 1000)
    d.tips = Any[(pos, txt, vis, placement)]
    on(events(d.figure).mouseposition) do mp
        hit = nothing
        for ax in vcat(d.charts, d.map_axis === nothing ? Axis[] : [d.map_axis])
            is_mouseinside(ax.scene) || continue
            hit = _hover(d, ax, mouseposition(ax.scene))
            break
        end
        if hit === nothing
            vis[] && (vis[] = false)
        else
            w, h = widths(d.figure.scene.viewport[])
            placement[] = mp[1] > 0.7w ? :left : mp[2] > 0.75h ? :below : :above
            pos[] = Point2f(mp)
            txt[] = hit[2]
            vis[] || (vis[] = true)
        end
        return Consume(false)
    end
    return d
end

# (anchor, text) for the cursor at `p` on `ax`, or nothing if the cursor is not
# on a site
function _hover(d::TKDash, ax::Axis, p)
    sites = d.survey.sites
    if ax === d.map_axis
        i = _map_hit(d, p)
        i === nothing && return nothing
        return (Point2f(sites[i].longitude, sites[i].latitude), _site_tip(d, i))
    end
    k = findfirst(c -> c === ax, d.charts)
    (k === nothing || k > length(d.chart_rows)) && return nothing
    rows = d.chart_rows[k]
    y = round(Int, p[2])
    (1 <= y <= length(rows) && abs(p[2] - y) <= DASH_BAR_HALF + 0.1) || return nothing
    i = rows[y]
    tip = _site_tip(d, i)
    t = d.t0 + Millisecond(round(Int, p[1] * 3.6e6))
    run = findfirst(r -> r.start <= t <= r.stop, sites[i].runs)
    if run !== nothing
        r = sites[i].runs[run]
        tip *= "\nrun $(basename(r.path)), $(_fs_label(r.sample_rate))\n" *
               Dates.format(r.start, "d u HH:MM") * " – " * Dates.format(r.stop, "d u HH:MM")
    end
    return (Point2f(p[1], y - DASH_BAR_HALF), tip)
end

# The name of a chart row. The row of a telluric site tells that it has no
# magnetic field
_row_name(s::SurveySite) = has_magnetic(s) ? s.name : s.name * " · E only"

function _site_tip(d::TKDash, i::Integer)
    s = d.survey.sites[i]
    own = recording_seconds(s; rate = d.rate)
    pos = isfinite(s.latitude) ? @sprintf("%.4f°, %.4f°", s.latitude, s.longitude) : "no position"
    channels = "channels: " * _channel_names(site_components(s)) *
               (has_magnetic(s) ? "" : " (no Hx, Hy: never a base or remote site)")
    (d.focus == 0 || i == d.focus) && return "$(s.name)\nown $(_hours(own)) · $pos\n$channels"
    focus = d.survey.sites[d.focus]
    km = site_distance(s, focus)
    dist = isfinite(km) ? @sprintf("%.2f km from %s", km, focus.name) : "no position"
    for (role, list) in (("base", d.refs.base), ("remote", d.refs.remote))
        c = findfirst(c -> c.site == s.name, list)
        c === nothing && continue
        c = list[c]
        return "$(s.name) · $role$(c.excluded ? " (dropped)" : "")\n$dist\n" *
               "overlap/own $(_bar_label(3600c.overlap_hours, own))\n$channels"
    end
    return "$(s.name)\n$dist\noverlap/own $(_bar_label(overlap_seconds(s, focus; rate = d.rate), own))\n$channels"
end

# The site below the cursor on the map, within a few percent of the frame
function _map_hit(d::TKDash, p)
    lim = d.map_axis.finallimits[]
    lat0 = lim.origin[2] + lim.widths[2] / 2
    tol = 0.04 * max(lim.widths[1] * cosd(lat0), lim.widths[2])
    best, best_d = nothing, tol
    for (i, s) in enumerate(d.survey.sites)
        isfinite(s.latitude) || continue
        dist = hypot((s.longitude - p[1]) * cosd(lat0), s.latitude - p[2])
        dist < best_d && ((best, best_d) = (i, dist))
    end
    return best
end

#---------- running -----

"""
    run_tkdash(root; kwargs...) -> TKDash
    run_tkdash(survey::Survey; kwargs...) -> TKDash
    run_tkdash(; kwargs...) -> TKDash

Open the survey dashboard and wait until you close its window. Then return the
[`TKDash`](@ref) with the choices that you made in it. The function scans
`root` with [`scan_survey`](@ref). It gives the keywords to [`TKDash`](@ref).
If you do not give `root`, a folder dialog asks for it.

```julia
dash = run_tkdash("data/survey")       # pick sites, drop references, close the window
write_reference_plan("plan.txt", dash)
```
"""
function run_tkdash(d::TKDash)
    @info "Opening TKDash" sites = length(d.survey.sites)
    screen = display(GLMakie.Screen(; title = "TKDash", focus_on_show = true), d.figure)
    _place_dash_window!(screen, _header_width(d))
    try
        wait(screen)
    catch
    end
    @info "TKDash closed"
    return d
end

# The width that the header needs to show all its buttons and boxes: the
# widths of its controls, the gaps between them and the padding of the figure
_header_width(d::TKDash) = _header_width(d.header)

function _header_width(h::GridLayout)
    w = 0.0
    for c in h.content
        inner = c.content.layoutobservables.reporteddimensions[].inner[1]
        inner === nothing || (w += inner)
    end
    return w + 6 * (h.size[2] - 1) + 32
end

# A window, not the full screen: 80 % of the monitor in each direction, but
# never narrower than the header needs (up to the width of the monitor), in the
# center. The layout changes when you resize the window
function _place_dash_window!(screen, min_width::Real = 0)
    try
        glwin = screen.glscreen
        vmode = GLMakie.GLFW.GetVideoMode(GLMakie.GLFW.GetPrimaryMonitor())
        mw, mh = Int(vmode.width), Int(vmode.height)
        w = round(Int, clamp(max(0.8mw, min_width), 0, 0.98mw))
        h = round(Int, 0.8mh)
        GLMakie.GLFW.SetWindowSize(glwin, w, h)
        GLMakie.GLFW.SetWindowPos(glwin, max((mw - w) ÷ 2, 0), max((mh - h) ÷ 2, 0))
    catch err
        @debug "Could not place the TKDash window" exception = err
    end
    return screen
end

run_tkdash(survey::Survey; kwargs...) = run_tkdash(TKDash(survey; kwargs...))

function run_tkdash(root::AbstractString; kwargs...)
    @info "Scanning survey" root
    return run_tkdash(TKDash(root; kwargs...))
end

function run_tkdash(; kwargs...)
    dir = pick_folder()
    isempty(dir) && error("No survey directory chosen")
    return run_tkdash(dir; kwargs...)
end
