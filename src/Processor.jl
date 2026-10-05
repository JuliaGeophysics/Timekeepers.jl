# Processor.jl - the transfer function window, TKProc.
# Author: @pankajkmishra
#
# A GLMakie window that estimates the transfer function of one site at a time
# (Processing.jl) and writes it as EDI, ModEM and a plot (EDI.jl, ModEM.jl).
# It has:
# - a header: Load Site, the site, its base site, its remote site and the rate
# - a row of processing options
# - a row with the Full tensor switch and the legend
# - the apparent resistivity above the phase, and Tzx above Tzy, against the
#   period, with the error bars
# - the polarity check of the estimate (Polarity.jl)
# - a status line that follows the processing
#
# Load Site opens a site directory. The window scans the survey around it (the
# directory above the site) and fills the base and remote menus with the sites
# that recorded with it, nearest and longest first. It selects the set-up of
# default_references: the best remote site, and for a site without Hx and Hy,
# the best base site.
#
# If you selected the base and remote sites of each site in TKDash and
# exported them, the window uses that plan: the reference_plan.txt in the
# survey directory or up to two directories above it. With the TKDash plan
# switch on, the menus hold only the sites of the plan for the site, and the
# first base site and the first remote site are selected. Each Process uses
# one combination: one base site (or none) and one remote site (or none). The
# switch off, or no plan for the site, gives the lists that the window
# calculates. Process estimates the site in a task of its own. Thus,
# the window stays live and the status line shows each step. Below the plots,
# the check of the channels (Polarity.jl) tells if a channel looks flipped.
# After Process, it reports only what the quadrants of Zxy and Zyx show.
# FlipCheck then compares H of the inputs with the remote sites and up to two
# base sites (witnesses, which take no part in the estimate) to tell which
# channel is flipped.
# Export writes the EDI, the ModEM file, a plot and a record of the processing
# (Report.jl) of the site
#
# The plots: Zxy red, Zyx blue, Zxx green and Zyy lilac, circles with a black
# edge over black error bars. Each phase is as recorded, on -200° to 200°. The
# real tipper is red and the imaginary tipper is blue

const PROC_ZCOLOURS = (RGBf(0.62, 0.84, 0.60), RGBf(0.84, 0.30, 0.10),
                       RGBf(0.12, 0.38, 0.72), RGBf(0.78, 0.66, 0.88))     # xx, xy, yx, yy
const PROC_RE = RGBf(0.84, 0.30, 0.10)
const PROC_IM = RGBf(0.12, 0.38, 0.72)
const PROC_MARKER = (marker = :circle, markersize = 10, strokecolor = :black, strokewidth = 1.6)
const PROC_PHASE_RANGE = (-200.0, 200.0)
const PROC_HINT = "Load Site opens a site with its base and remote sites · Process estimates it · " *
                  "FlipCheck finds a flipped channel · Export writes its EDI, ModEM file, plot and record · " *
                  "Clear empties the screen."

# The text of the ⓘ next to each control
const PROC_HELP = Dict(
    :base => "The site that gives the magnetic inputs Hx, Hy. Use a base site for a site without " *
             "magnetic channels, or with poor ones. With none, the site gives its own Hx, Hy.",
    :remote => "The far site whose Hx, Hy are the reference. Noise in the local magnetic field is not " *
               "at the remote site. Thus, the estimate is free of the downward bias of a single site. " *
               "With none, the estimate is single site.",
    :rate => "The sampling rate to process. All rates estimates each rate that all the sites recorded " *
             "and keeps the estimate with the smallest error at each period. Bursts at a high rate " *
             "then widen the spectrum at the short periods.",
    :plan => "Use the base and remote sites that you selected in TKDash (reference_plan.txt from its " *
             "Export). The menus then hold only these sites, one combination at a time. Off: the " *
             "lists that TKProc calculates.",
    :window => "Samples in each FFT window, at every decimation level (a power of 2). Larger: finer " *
               "frequency resolution and more harmonics in each band, but fewer windows, so the " *
               "robust weights and the errors have less data. Smaller: more windows, coarser " *
               "resolution. 256 suits most records.",
    :overlap => "The part of each window that the next window shares. Higher: more windows, but they are " *
                "not independent, so the gain is small above 0.5. Lower: fewer windows. 0.5 suits the " *
                "Hann window.",
    :bands => "Periods in each decade of the output. Higher: more periods and a finer curve, with fewer " *
              "harmonics in each, so larger errors. Lower: fewer periods with smaller errors.",
    :prewhiten => "The order of the AR filter that flattens each spectrum before the FFT. The " *
                  "correction after the FFT is exact. Higher: less leakage from the strong long periods " *
                  "into the short ones. 0: no filter. 2 to 5 suits most records.",
    :top => "The highest frequency, as a fraction of the Nyquist frequency of the record. auto: " *
           "where the anti-alias filter of the recorder ends the signal, found from the coherence " *
           "of E with H (it falls below 0.5 there), less 5 %, from 0.5 to 0.9. A number sets it: " *
           "higher gives shorter periods, but closer to the filter, where the data can be biased.",
    :harmonic => "The lowest FFT harmonic that the last decimation level uses. Lower (2 or 3): longer " *
                 "periods, but closer to the trend of each window, with more leakage. Higher: safer " *
                 "long periods, a shorter spectrum.",
    :windows => "The windows that a decimation level needs. Lower: more levels and longer periods, " *
                "with few windows, so larger errors and less robust weights. Higher: only well " *
                "determined periods.",
    :periods => "Estimate only the periods in this range, in seconds. Inf: no upper limit.",
    :huber => "Where the Huber weights start: residuals larger than this many robust standard " *
              "deviations get less weight. Lower: more outliers removed, but also more good data, so " *
              "larger errors. Higher: closer to least squares. 1.5 is the usual value.",
    :jackknife => "The blocks of windows that the jackknife deletes one at a time to estimate the " *
                  "errors. More: smoother errors when there are many windows. Fewer: each block is " *
                  "longer and more independent. The number of windows limits it.",
    :method => "CT2004 (after Chave & Thomson 2004): two stages with bounded influence. Stage 1 " *
               "fits the local H to the reference and removes noise bursts in H; stage 2 uses Huber, " *
               "then redescending weights, and caps windows of high leverage. Best with noise in the " *
               "local H or extreme windows. EB1986 (after Egbert & Booker 1986): one robust " *
               "M-estimate with Huber weights only, the reference as it is. Keeps more data with a far " *
               "or weakly coherent remote site, and keeps the strong events. On clean data both agree.",
    :leverage => "Give less weight to windows with extreme magnetic inputs (storms, noise in H) that " *
                 "would control the fit alone. On: bounded influence, safer. Off: all the windows count.",
    :full => "Show Zxx and Zyy too, not only Zxy and Zyx. The export plot follows this switch.",
    :bars => "Show the error bars (one standard deviation). The export plot follows this switch.",
)

# Break a help text into lines of up to `n` characters
function _wrap(text::AbstractString, n::Integer = 64)
    lines, line = String[], ""
    for w in split(text)
        if !isempty(line) && length(line) + 1 + length(w) > n
            push!(lines, line)
            line = String(w)
        else
            line = isempty(line) ? String(w) : line * " " * w
        end
    end
    isempty(line) || push!(lines, line)
    return join(lines, "\n")
end

"""
    TKProc

The transfer function window. It holds:
- the [`Survey`](@ref) around the loaded site;
- the loaded site;
- the [`TransferFunction`](@ref) of each site that you processed, in
  `results`;
- the GLMakie figure.

To make one, use `TKProc(site_dir)` and `display` it, or use
[`run_tkproc`](@ref). After you close the window, `proc.results` holds the
estimates.
"""
mutable struct TKProc
    figure::Figure
    header::GridLayout
    survey::Survey
    focus::Int
    results::Dict{String, TransferFunction}
    site_label::Label
    base_menu::Menu
    remote_menu::Menu
    rate_menu::Menu
    boxes::Dict{Symbol, Textbox}
    method_menu::Menu
    leverage_toggle::Toggle
    full_toggle::Toggle
    bars_toggle::Toggle
    plan_toggle::Toggle
    plan_source::Any
    plan::Vector{Any}
    plan_path::String
    process_button::Button
    flip_button::Button
    axes::Vector{Axis}
    status::Label
    polarity::Label
    updating::Bool
    busy::Bool
end

"""
    TKProc(site_dir; plan = nothing, size = (1600, 950)) -> TKProc
    TKProc(survey::Survey, site; plan = nothing, size = (1600, 950)) -> TKProc
    TKProc(; plan = nothing, size = (1600, 950)) -> TKProc

Make the transfer function window for one site. With `site_dir`, the function
scans the directory above it with [`scan_survey`](@ref) to find the base and
remote sites. With a [`Survey`](@ref), give the site by its name or index.
Without a site, the window opens empty: use its Load Site button. The function
does not open a window.

`plan` selects the base and remote sites of a TKDash plan
([`write_reference_plan`](@ref)):
- `nothing` -- the `reference_plan.txt` in the survey directory or up to two
  directories above it, if there is one;
- a path -- that plan file;
- `false` -- no plan: the window calculates the lists.
"""
function TKProc(site_dir::AbstractString; kwargs...)
    p = TKProc(; kwargs...)
    _load_proc_site!(p, site_dir) || error(p.status.text[])
    return p
end

function TKProc(survey::Survey, site; kwargs...)
    p = TKProc(; kwargs...)
    p.survey = survey
    _load_plan!(p, survey.root)
    t = _target_site(survey, site)
    _proc_focus!(p, findfirst(s -> s === t, survey.sites))
    return p
end

function TKProc(; plan = nothing, size = (1600, 950))
    GLMakie.activate!(title = "TKProc")
    fig = Figure(; size = size, figure_padding = (16, 16, 10, 10))

    header = GridLayout(fig[1, 1]; tellwidth = false)
    b_load = Button(header[1, 1]; label = "Load Site…")
    site_label = Label(header[1, 2], "no site"; color = DASH_NAVY, font = :bold, padding = (8, 8, 0, 0))
    help = Tuple{Any, String}[]
    hint(label, key) = (push!(help, (label, _wrap(PROC_HELP[key]))); label)
    hint(Label(header[1, 3], "Base ⓘ"; padding = (10, 0, 0, 0)), :base)
    base_menu = _dash_menu(header[1, 4]; options = [("none", nothing)], width = 150)
    hint(Label(header[1, 5], "Remote ⓘ"; padding = (10, 0, 0, 0)), :remote)
    remote_menu = _dash_menu(header[1, 6]; options = [("none", nothing)], width = 170)
    hint(Label(header[1, 7], "Rate ⓘ"; padding = (10, 0, 0, 0)), :rate)
    rate_menu = _dash_menu(header[1, 8]; options = [("All rates", :all)], width = 100)
    plan_toggle = Toggle(header[1, 9]; active = false)
    hint(Label(header[1, 10], "TKDash plan ⓘ"), :plan)
    Box(header[1, 11]; visible = false)
    b_run = Button(header[1, 12]; label = "Process", width = 110)
    b_flip = Button(header[1, 13]; label = "FlipCheck", width = 100)
    b_export = Button(header[1, 14]; label = "Export…")
    b_clear = Button(header[1, 15]; label = "Clear")
    colsize!(header, 11, Auto(true, 1.0))
    colgap!(header, 6)

    # the processing options, with the keywords of estimate_tf: the spectrum in
    # the first row, the robust fit in the second
    opts = GridLayout(fig[2, 1]; tellwidth = false, halign = :left)
    rows = [GridLayout(opts[1, 1]; halign = :left), GridLayout(opts[2, 1]; halign = :left)]
    boxes = Dict{Symbol, Textbox}()
    col = [0, 0]
    next!(row) = (col[row] += 1)
    function field!(row, key, label, value, check, help_key; width = 52)
        c = next!(row)
        hint(Label(rows[row][1, c], label * " ⓘ"; padding = (c == 1 ? 0 : 12, 0, 0, 0), halign = :right), help_key)
        boxes[key] = Textbox(rows[row][1, next!(row)]; stored_string = value, width = width, validator = check,
                             halign = :left)
    end
    function switch!(row, label, help_key; active)
        t = Toggle(rows[row][1, next!(row)]; active = active)
        hint(Label(rows[row][1, next!(row)], label * " ⓘ"; halign = :left), help_key)
        return t
    end
    isint(lo, hi) = s -> (v = tryparse(Int, s); v !== nothing && lo <= v <= hi)
    isreal(lo, hi) = s -> (v = tryparse(Float64, s); v !== nothing && lo <= v <= hi)
    field!(1, :window, "Window", "256", s -> (v = tryparse(Int, s); v !== nothing && 64 <= v <= 65536 && ispow2(v)), :window)
    field!(1, :overlap, "Overlap", "0.5", isreal(0.0, 0.9), :overlap)
    field!(1, :bands, "Bands/decade", "8", isint(2, 20), :bands)
    field!(1, :prewhiten, "AR order", "3", isint(0, 20), :prewhiten)
    field!(1, :top, "Top f ×Nyquist", "auto", s -> lowercase(strip(s)) == "auto" || isreal(0.1, 0.9)(s), :top)
    field!(1, :harmonic, "Lowest harmonic", "4", isint(1, 64), :harmonic)
    field!(1, :windows, "Min windows", "8", isint(3, 1000), :windows)
    field!(1, :min_period, "Periods (s)", "0", isreal(0.0, Inf), :periods; width = 60)
    Label(rows[1][1, next!(1)], "to")
    boxes[:max_period] = Textbox(rows[1][1, next!(1)]; stored_string = "Inf", width = 60,
                                 validator = isreal(0.0, Inf), halign = :left)
    hint(Label(rows[2][1, next!(2)], "Method ⓘ"; halign = :right), :method)
    method_menu = _dash_menu(rows[2][1, next!(2)];
                             options = [("CT2004", :ct2004), ("EB1986", :eb1986)], width = 100)
    field!(2, :huber, "Huber", "1.5", isreal(0.5, 10.0), :huber)
    field!(2, :jackknife, "Jackknife blocks", "50", isint(3, 1000), :jackknife)
    leverage_toggle = switch!(2, "Leverage weights", :leverage; active = true)
    foreach(r -> colgap!(r, 6), rows)
    rowgap!(opts, 6)

    view = GridLayout(fig[3, 1]; tellwidth = false, halign = :left)
    full_toggle = Toggle(view[1, 1]; active = false)
    hint(Label(view[1, 2], "Full tensor ⓘ"; padding = (0, 16, 0, 0)), :full)
    bars_toggle = Toggle(view[1, 3]; active = true)
    hint(Label(view[1, 4], "Error bars ⓘ"; padding = (0, 16, 0, 0)), :bars)
    _proc_legend!(view[1, 5])
    colgap!(view, 6)

    body = GridLayout(fig[4, 1])
    axes = _tf_axes!(body)
    polarity = Label(fig[5, 1], ""; fontsize = 12, halign = :left, tellwidth = false, justification = :left)
    status = Label(fig[6, 1], ""; fontsize = 12, halign = :left, tellwidth = false)
    colsize!(fig.layout, 1, Relative(1))
    rowsize!(fig.layout, 4, Auto(true, 1.0))
    rowgap!(fig.layout, 8)

    p = TKProc(fig, header, Survey("", SurveySite[]), 0, Dict{String, TransferFunction}(), site_label,
               base_menu, remote_menu, rate_menu, boxes, method_menu, leverage_toggle,
               full_toggle, bars_toggle, plan_toggle, plan, Any[], "", b_run, b_flip, axes, status, polarity, false, false)

    on(b_load.clicks) do _
        p.busy && return
        dir = try
            pick_folder()
        catch err
            _proc_status!(p, "Could not open a folder dialog: $(sprint(showerror, err))"; error = true)
            ""
        end
        isempty(dir) || _load_proc_site!(p, dir)
    end
    on(_ -> _process_focus!(p), b_run.clicks)
    on(_ -> _export_site!(p), b_export.clicks)
    on(_ -> _flip_check!(p), b_flip.clicks)
    on(_ -> _clear_proc!(p), b_clear.clicks)
    on(_ -> _draw_proc!(p), full_toggle.active)
    on(_ -> _draw_proc!(p), bars_toggle.active)
    _install_proc_help!(fig, help)
    on(plan_toggle.active) do active
        p.updating && return
        if active && isempty(p.plan)
            p.updating = true
            plan_toggle.active[] = false
            p.updating = false
            return _proc_status!(p, "No TKDash plan: export one from TKDash (reference_plan.txt in the survey directory)."; error = true)
        end
        p.focus == 0 || _proc_focus!(p, p.focus)
    end
    _draw_proc!(p)
    _proc_status!(p, "")
    return p
end

function Base.display(p::TKProc)
    display(p.figure)
    return p
end

# FlipCheck: compare H of the inputs with the remote sites and up to two base
# sites, in a task of its own, then show the full check
function _flip_check!(p::TKProc)
    (p.busy || p.focus == 0) && return nothing
    name = p.survey.sites[p.focus].name
    tf = get(p.results, name, nothing)
    tf === nothing && return _proc_status!(p, "Process $name before FlipCheck."; error = true)
    p.busy = true
    p.flip_button.label[] = "Checking…"
    return @async begin
        messages = Channel{String}(Inf)
        t0 = time()
        try
            task = Threads.@spawn flip_check!(p.survey, tf; progress = msg -> (put!(messages, msg); yield()))
            step = "comparing H with the other sites"
            while !istaskdone(task)
                while isready(messages)
                    step = take!(messages)
                end
                _proc_status!(p, @sprintf("FlipCheck %s · %s · %.0f s", name, step, time() - t0); busy = true)
                sleep(0.25)
            end
            try
                fetch(task)
            catch err
                throw(err isa TaskFailedException ? err.task.exception : err)
            end
            _show_polarity!(p, tf)
            others = vcat(tf.remote, get(tf.metadata, :witnesses, String[]))
            _proc_status!(p, isempty(others) ? "FlipCheck: no other site recorded H with $name." :
                             "FlipCheck: compared H of $(isempty(tf.base) ? name : tf.base) with " *
                             join(others, ", ") * @sprintf(" in %.1f s.", time() - t0))
        catch err
            _proc_status!(p, "FlipCheck failed: $(sprint(showerror, err))"; error = true)
        finally
            p.busy = false
            p.flip_button.label[] = "FlipCheck"
        end
    end
end

# The polarity check of the estimate on the screen, below the plots: grey when
# it finds nothing, amber when a channel looks reversed
function _show_polarity!(p::TKProc, tf)
    if tf === nothing
        p.polarity.text[] = ""
        return p
    end
    c = check_polarity(tf; magnetic = get(tf.metadata, :flipcheck, false))
    p.polarity.color[] = c.ok ? TK_GREY : DASH_REMOTE_RAMP[end]
    p.polarity.text[] = _wrap((get(tf.metadata, :flipcheck, false) ? "FlipCheck: " : "Channels: ") * c.message, 190)
    return p
end

# One help box for the window: a white box with the help of the label below
# the cursor, drawn in pixels above all the panels. It opens below the cursor
# and moves left near the right edge and up near the bottom of the window
function _install_proc_help!(fig, help)
    pos, txt = Observable(Point2f(0, 0)), Observable(" ")
    box = Observable(Rect2f(0, 0, 1, 1))
    bg = poly!(fig.scene, box; color = RGBAf(1, 1, 1, 0.97), strokecolor = RGBAf(0, 0, 0, 0.45),
               strokewidth = 0.75, space = :pixel, visible = false)
    tx = text!(fig.scene, pos; text = txt, space = :pixel, fontsize = 12, color = TK_BLACK,
               align = (:left, :top), visible = false)
    translate!(bg, 0, 0, 1000)
    translate!(tx, 0, 0, 1001)
    shown = Ref(false)
    on(events(fig).mouseposition) do mp
        hit = nothing
        for (label, text) in help
            bb = label.layoutobservables.computedbbox[]
            (bb.origin[1] <= mp[1] <= bb.origin[1] + bb.widths[1] &&
             bb.origin[2] <= mp[2] <= bb.origin[2] + bb.widths[2]) || continue
            hit = text
            break
        end
        if hit === nothing
            shown[] && (bg.visible[] = tx.visible[] = shown[] = false)
            return Consume(false)
        end
        txt[] = hit
        pos[] = Point2f(mp[1] + 12, mp[2] - 18)
        tb = boundingbox(tx, :pixel)
        w, h = widths(fig.scene.viewport[])
        x = mp[1] + 12 + tb.widths[1] + 10 > w ? mp[1] - tb.widths[1] - 22 : mp[1] + 12
        y = mp[2] - 18 - tb.widths[2] - 8 < 0 ? mp[2] + tb.widths[2] + 18 : mp[2] - 18
        pos[] = Point2f(x, y)
        box[] = Rect2f(x - 8, y - tb.widths[2] - 6, tb.widths[1] + 16, tb.widths[2] + 12)
        shown[] || (bg.visible[] = tx.visible[] = shown[] = true)
        return Consume(false)
    end
    return fig
end

#---------- the plots -----

# The four period axes: ρa above φ on the left, Tzx above Tzy on the right.
# They keep the light grey panel and have no grid
function _tf_axes!(grid)
    axis(pos; kw...) = Axis(pos; xscale = log10, backgroundcolor = DASH_PANEL, xgridvisible = false,
                            ygridvisible = false, xlabelfont = :regular, ylabelfont = :regular, kw...)
    ax_rho = axis(grid[1, 1]; yscale = log10, ylabel = "Apparent resistivity (Ω·m)",
                  xticklabelsvisible = false)
    ax_phi = axis(grid[2, 1]; ylabel = "Phase (°)", xlabel = "Period (s)", yticks = -180:90:180)
    ax_tzx = axis(grid[1, 2]; ylabel = "Tzx", xticklabelsvisible = false)
    ax_tzy = axis(grid[2, 2]; ylabel = "Tzy", xlabel = "Period (s)")
    axes = [ax_rho, ax_phi, ax_tzx, ax_tzy]
    linkxaxes!(axes...)
    colgap!(grid, 24)
    rowgap!(grid, 8)
    return axes
end

function _proc_legend!(pos)
    items = Any[MarkerElement(; PROC_MARKER..., color = PROC_ZCOLOURS[i]) for i in (2, 3, 1, 4)]
    append!(items, [MarkerElement(; PROC_MARKER..., color = PROC_RE), MarkerElement(; PROC_MARKER..., color = PROC_IM)])
    return Legend(pos, items, ["Zxy", "Zyx", "Zxx", "Zyy", "Re T", "Im T"]; orientation = :horizontal,
                  framevisible = false, labelsize = 12, patchsize = (14, 12))
end

# Decades as 10ⁿ when two or more fall in the range, otherwise 1-2-5 steps in
# plain numbers
function _decade_ticks(lo, hi)
    decades = ceil(Int, log10(lo)):floor(Int, log10(hi))
    length(decades) >= 2 && return (10.0 .^ decades, [rich("10", superscript(string(n))) for n in decades])
    v = [m * 10.0^n for n in floor(Int, log10(lo)):ceil(Int, log10(hi)) for m in (1, 2, 5) if lo <= m * 10.0^n <= hi]
    return (v, [x >= 10 ? @sprintf("%.0f", x) : @sprintf("%.2g", x) for x in v])
end

_proc_note!(ax, s) = text!(ax, 0.02, 0.96; text = s, space = :relative, align = (:left, :top), fontsize = 12,
                           color = :grey30)

# One series: black error bars under circles with a black edge
function _tf_series!(ax, T, y, lo, hi, colour; bars::Bool = true)
    ok = findall(isfinite, y)
    isempty(ok) && return Float64[]
    e = filter(i -> isfinite(lo[i]) && isfinite(hi[i]), ok)
    (isempty(e) || !bars) || errorbars!(ax, T[e], y[e], lo[e], hi[e]; color = :black, linewidth = 1, whiskerwidth = 6)
    scatter!(ax, T[ok], y[ok]; PROC_MARKER..., color = colour)
    return y[ok]
end

# The impedance and the tipper of `tf` on the four axes. `full_tensor` adds
# Zxx and Zyy to Zxy and Zyx
function _draw_tf!(axes, tf::TransferFunction; full_tensor::Bool = false, bars::Bool = true)
    ax_rho, ax_phi, ax_tzx, ax_tzy = axes
    foreach(empty!, axes)
    T = tf.periods
    isempty(T) && return nothing
    σ = sqrt.(tf.Z_var ./ 2)
    seen = Float64[]
    for (i, j, c) in (full_tensor ? ((1, 2, 2), (2, 1, 3), (1, 1, 1), (2, 2, 4)) : ((1, 2, 2), (2, 1, 3)))
        z = tf.Z[i, j, :]
        rho = [isfinite(v) && abs(v) > 0 ? 0.2 * t * abs2(v) : NaN for (v, t) in zip(z, T)]
        phi = [isfinite(r) ? rad2deg(angle(v)) : NaN for (v, r) in zip(z, rho)]
        rel = σ[i, j, :] ./ abs.(z)
        δρ, δφ = 2 .* rho .* rel, min.(rad2deg.(rel), 90)
        append!(seen, _tf_series!(ax_rho, T, rho, min.(δρ, 0.95 .* rho), δρ, PROC_ZCOLOURS[c]; bars))
        _tf_series!(ax_phi, T, phi, δφ, δφ, PROC_ZCOLOURS[c]; bars)
    end
    xl = (first(T) / 1.5, last(T) * 1.5)
    rl = isempty(seen) ? (1.0, 1.0e4) : (10^(log10(minimum(seen)) - 0.5), 10^(log10(maximum(seen)) + 0.5))
    isempty(seen) && _proc_note!(ax_rho, "no impedance at this site")
    ax_rho.yticks = _decade_ticks(rl...)
    ax_rho.limits[] = (xl, rl)
    ax_phi.limits[] = (xl, PROC_PHASE_RANGE)
    tv = sqrt.(tf.T_var ./ 2)
    for (j, ax) in ((1, ax_tzx), (2, ax_tzy))
        t = tf.T[j, :]
        if !any(isfinite, t)
            _proc_note!(ax, "no tipper at this site")
            ax.limits[] = (xl, (-0.5, 0.5))
            continue
        end
        hlines!(ax, [0.0]; color = :grey70, linewidth = 1)
        top = 0.0
        for (part, colour) in ((real, PROC_RE), (imag, PROC_IM))
            y = [isfinite(v) ? part(v) : NaN for v in t]
            shown = _tf_series!(ax, T, y, tv[j, :], tv[j, :], colour; bars)
            top = max(top, maximum(abs, shown; init = 0.0))
        end
        top = max(top, 0.2)
        ax.limits[] = (xl, (-1.1top, 1.1top))
    end
    for ax in axes
        ax.xticks = _decade_ticks(xl...)
        reset_limits!(ax)
    end
    return nothing
end

"""
    plot_tf(tf::TransferFunction; full_tensor = false, errors = true, size = (1300, 800)) -> Figure

A figure of `tf`: the apparent resistivity above the phase on the left, the
tipper Tzx above Tzy on the right, against the period, with the error bars.
Zxy is red and Zyx blue. With `full_tensor`, Zxx (green) and Zyy (lilac) are
there too. With `errors = false`, no error bars. The real tipper is red and
the imaginary tipper is blue. Save it
with `save("site.png", plot_tf(tf))`.
"""
function plot_tf(tf::TransferFunction; full_tensor::Bool = false, errors::Bool = true, size = (1300, 800))
    fig = Figure(; size = size, figure_padding = (16, 20, 10, 10))
    top = GridLayout(fig[1, 1]; tellwidth = false, halign = :left)
    Label(top[1, 1], isempty(tf.periods) ? tf.site : "$(tf.site) · $(_setup_text(tf))"; font = :bold,
          padding = (0, 20, 0, 0))
    _proc_legend!(top[1, 2])
    axes = _tf_axes!(GridLayout(fig[2, 1]))
    _draw_tf!(axes, tf; full_tensor, bars = errors)
    return fig
end

function _draw_proc!(p::TKProc)
    name = p.focus == 0 ? "" : p.survey.sites[p.focus].name
    tf = get(p.results, name, nothing)
    _show_polarity!(p, tf)
    if tf === nothing
        foreach(empty!, p.axes)
        _proc_note!(p.axes[1], isempty(name) ? "no site loaded" : "$name: not processed")
        for ax in p.axes
            ax.limits[] = ((1.0e-3, 1.0e3), ax === p.axes[1] ? (0.1, 1000.0) : ax === p.axes[2] ? PROC_PHASE_RANGE : (-0.5, 0.5))
            ax.xticks = _decade_ticks(1.0e-3, 1.0e3)
            reset_limits!(ax)
        end
        p.axes[1].yticks = _decade_ticks(0.1, 1000.0)
        return p
    end
    _draw_tf!(p.axes, tf; full_tensor = p.full_toggle.active[], bars = p.bars_toggle.active[])
    return p
end

#---------- state changes -----

# Load the site in `dir` and the survey around it: the directory above the
# site, or the survey that is loaded if it holds the site already
function _load_proc_site!(p::TKProc, dir::AbstractString)
    path = _norm_path(dir)
    site = isdir(path) ? _scan_site(path) : nothing
    if site === nothing
        _proc_status!(p, "$path is not a site: give a directory with Metronix, LEMI-424 or GEOMAG recordings."; error = true)
        return false
    end
    if !any(c -> c in site_components(site), (:e1, :e2))
        _proc_status!(p, "$(site.name) has no electric channels: it can be a base or remote site, not a target."; error = true)
        return false
    end
    i = findfirst(s -> s.path == site.path, p.survey.sites)
    if i === nothing
        root = dirname(site.path)
        _proc_status!(p, "Scanning $root for the base and remote sites of $(site.name) …")
        survey = try
            scan_survey(root)
        catch err
            _proc_status!(p, "Could not scan $root: $(sprint(showerror, err))"; error = true)
            return false
        end
        p.survey = survey
        _load_plan!(p, root)
        i = findfirst(s -> s.path == site.path, survey.sites)
        if i === nothing
            p.survey = Survey(root, vcat(survey.sites, site))
            i = length(p.survey.sites)
        end
    end
    _proc_focus!(p, i)
    return true
end

# The TKDash plan: the file that `plan_source` names, or the first
# reference_plan.txt in `root` and up to two directories above it. The switch
# is on when the window finds a plan
function _load_plan!(p::TKProc, root::AbstractString)
    p.plan, p.plan_path = Any[], ""
    if p.plan_source !== false && !isempty(root)
        candidates = p.plan_source isa AbstractString ? [String(p.plan_source)] :
                     [joinpath(d, "reference_plan.txt") for d in unique([root, dirname(root), dirname(dirname(root))])]
        for f in candidates
            isfile(f) || continue
            try
                p.plan = collect(Any, read_reference_plan(f))
                p.plan_path = f
                break
            catch err
                @warn "Could not read the TKDash plan $f" exception = err
            end
        end
    end
    p.updating = true
    try
        # a Toggle animates at each change of `active`, also to the same value
        want = !isempty(p.plan)
        p.plan_toggle.active[] == want || (p.plan_toggle.active[] = want)
    finally
        p.updating = false
    end
    return p
end

# The plan entry of a site, or nothing if the switch is off or the plan does
# not list the site
function _plan_entry(p::TKProc, name)
    (p.plan_toggle.active[] && !isempty(p.plan)) || return nothing
    k = findfirst(r -> r.site == name, p.plan)
    return k === nothing ? nothing : p.plan[k]
end

function _proc_focus!(p::TKProc, i::Integer)
    p.busy && return p
    p.focus = i
    site = p.survey.sites[i]
    entry = _plan_entry(p, site.name)
    if entry === nothing
        refs = site_references(p.survey, site)
        defaults = default_references(p.survey, site)
        km(c) = @sprintf("%s · %.1f km", c.site, c.distance_km)
        base_opts = vcat(Any[("none", nothing)], Any[(km(c), c.site) for c in refs.base])
        remote_opts = vcat(Any[("none", nothing)], Any[(km(c), c.site) for c in refs.remote])
        length(refs.remote) > 1 && push!(remote_opts, ("all remotes", :all))
        nb, nr = length(refs.base), length(refs.remote)
        source = "recorded with it"
    else
        # one combination at a time: the menus hold the sites of the plan, no
        # "all remotes"
        function label(name, hours)
            j = findfirst(s -> s.name == name, p.survey.sites)
            d = j === nothing ? NaN : site_distance(site, p.survey.sites[j])
            return isfinite(d) ? @sprintf("%s · %.1f km · %.1f h", name, d, hours) : @sprintf("%s · %.1f h", name, hours)
        end
        base_opts = vcat(Any[("none", nothing)], Any[(label(n, h), n) for (n, h) in zip(entry.base, entry.base_hours)])
        remote_opts = vcat(Any[("none", nothing)], Any[(label(n, h), n) for (n, h) in zip(entry.remote, entry.remote_hours)])
        defaults = (base = isempty(entry.base) ? nothing : first(entry.base),
                    remote = isempty(entry.remote) ? nothing : first(entry.remote))
        nb, nr = length(entry.base), length(entry.remote)
        source = "in the TKDash plan $(basename(p.plan_path))"
    end
    rate_opts = vcat(Any[("All rates", :all)], Any[(_fs_label(fs), fs) for fs in site_rates(site)])
    pick(opts, v) = something(findfirst(o -> o[2] == v, opts), 1)
    p.updating = true
    try
        p.site_label.text[] = has_magnetic(site) ? site.name : site.name * " · E only"
        p.base_menu.options[] = base_opts
        p.base_menu.i_selected[] = pick(base_opts, defaults.base)
        p.remote_menu.options[] = remote_opts
        p.remote_menu.i_selected[] = pick(remote_opts, defaults.remote)
        p.rate_menu.options[] = rate_opts
        p.rate_menu.i_selected[] = 1
    finally
        p.updating = false
    end
    _draw_proc!(p)
    nothing_in_plan = entry === nothing && p.plan_toggle.active[] && !isempty(p.plan)
    _proc_status!(p, haskey(p.results, site.name) ? "Showing the estimate of $(site.name)." :
                     "Loaded $(site.name): $nb base site$(nb == 1 ? "" : "s") and $nr remote site$(nr == 1 ? "" : "s") " *
                     "$source$(nothing_in_plan ? " (the TKDash plan does not list it)" : ""). " *
                     "Select one base and one remote site, then Process.")
    return p
end

# The keywords of estimate_tf from the option row
function _proc_options(p::TKProc)
    get(key, T, default) = something(tryparse(T, p.boxes[key].stored_string[]), default)
    return (; window = get(:window, Int, 256), overlap = get(:overlap, Float64, 0.5),
            bands_per_decade = get(:bands, Int, 8), prewhiten = get(:prewhiten, Int, 3),
            huber = get(:huber, Float64, 1.5), jackknife_groups = get(:jackknife, Int, 50),
            nyquist_fraction = something(tryparse(Float64, p.boxes[:top].stored_string[]), :auto), min_harmonic = get(:harmonic, Int, 4),
            min_windows = get(:windows, Int, 8),
            min_period = get(:min_period, Float64, 0.0), max_period = get(:max_period, Float64, Inf),
            method = something(p.method_menu.selection[], :ct2004), leverage = p.leverage_toggle.active[])
end

"""
    _process_focus!(p) -> Union{Task, Nothing}

Estimate the loaded site with the menus and the options. The estimate runs in
a task on another thread. A task on the main thread shows its steps and the
time in the status line, then draws the result. The function returns that
task (`wait` on it to block).
"""
function _process_focus!(p::TKProc)
    (p.busy || p.focus == 0) && return nothing
    site = p.survey.sites[p.focus]
    base = p.base_menu.selection[]
    remote = p.remote_menu.selection[]
    remote === :all && (remote = [c.site for c in site_references(p.survey, site).remote])
    rate = something(p.rate_menu.selection[], :all)
    p.busy = true
    p.process_button.label[] = "Processing…"
    # a new estimate starts from an empty screen: the old one of the site
    # does not stay in view if this one fails
    delete!(p.results, site.name)
    _draw_proc!(p)
    return @async _run_estimate!(p, site.name; base, remote, rate)
end

# Clear: forget every estimate of the window and empty the plots, the check of
# the channels and the status line. The loaded site, its menus and the options
# stay
function _clear_proc!(p::TKProc)
    p.busy && return p
    n = length(p.results)
    empty!(p.results)
    _draw_proc!(p)
    _proc_status!(p, n == 0 ? "Nothing to clear." :
                     "Cleared $n estimate$(n == 1 ? "" : "s"). Press Process to estimate the site again.")
    return p
end

function _run_estimate!(p::TKProc, name; base, remote, rate = :all)
    messages = Channel{String}(Inf)
    t0 = time()
    setup = base === nothing && remote === nothing ? "single site" :
            join(filter(!isempty, [base === nothing ? "" : "base $base",
                                   remote === nothing ? "" : "remote $(join(vcat(remote), " + "))"]), ", ")
    try
        opts = _proc_options(p)
        task = Threads.@spawn estimate_tf(p.survey, name; base, remote, rate, opts...,
                                          progress = msg -> (put!(messages, msg); yield()))
        step = "starting"
        while !istaskdone(task)
            while isready(messages)
                step = take!(messages)
            end
            _proc_status!(p, @sprintf("Processing %s (%s) · %s · %.0f s", name, setup, step, time() - t0); busy = true)
            sleep(0.25)
        end
        tf = try
            fetch(task)
        catch err
            throw(err isa TaskFailedException ? err.task.exception : err)
        end
        p.results[name] = tf
        _draw_proc!(p)
        _proc_status!(p, @sprintf("%s: %d periods, %.3g–%.3g s, in %.1f s (%s).", name, length(tf.periods),
                                  first(tf.periods), last(tf.periods), time() - t0, _setup_text(tf)))
    catch err
        _proc_status!(p, "Could not process $name: $(sprint(showerror, err))"; error = true)
    finally
        p.busy = false
        p.process_button.label[] = "Process"
    end
    return p
end

_setup_text(tf::TransferFunction) =
    tf.mode === :single ? "single site" :
    tf.mode === :base ? "base $(tf.base)" :
    tf.mode === :remote ? "remote $(join(tf.remote, " + "))" :
    "base $(tf.base), remote $(join(tf.remote, " + "))"

# Export the site in focus: its EDI, its ModEM file and a PNG of the view on
# the screen, into a directory that you select
function _export_site!(p::TKProc)
    p.busy && return p
    p.focus == 0 && return _proc_status!(p, "Load and process a site first."; error = true)
    name = p.survey.sites[p.focus].name
    tf = get(p.results, name, nothing)
    tf === nothing && return _proc_status!(p, "Process $name before the export."; error = true)
    dir = try
        pick_folder(dirname(p.survey.sites[p.focus].path))
    catch err
        return _proc_status!(p, "Could not open a folder dialog: $(sprint(showerror, err))"; error = true)
    end
    (isempty(dir) || !isdir(dir)) && return p
    try
        written = export_tf(dir, tf; full_tensor = p.full_toggle.active[], errors = p.bars_toggle.active[])
        _proc_status!(p, "Wrote " * join(basename.(written), ", ") * " to $dir")
    catch err
        _proc_status!(p, "Could not export $name: $(sprint(showerror, err))"; error = true)
    end
    return p
end

"""
    export_tf(dir, tf::TransferFunction; full_tensor = false, errors = true) -> Vector{String}

Write the files of one site into `dir`: `<site>.edi` ([`write_edi`](@ref)),
`<site>.dat` ([`write_modem`](@ref)), `<site>.png` ([`plot_tf`](@ref),
with or without Zxx and Zyy and the error bars) and `<site>.md`, the record of the processing
([`write_tf_report`](@ref)). The errors are the estimated errors, with no
floor. The function returns the paths.
"""
function export_tf(dir::AbstractString, tf::TransferFunction; full_tensor::Bool = false, errors::Bool = true)
    mkpath(dir)
    edi = write_edi(joinpath(dir, tf.site * ".edi"), tf)
    dat = write_modem(joinpath(dir, tf.site * ".dat"), tf)
    png = joinpath(dir, tf.site * ".png")
    save(png, plot_tf(tf; full_tensor, errors); px_per_unit = 2)
    report = write_tf_report(joinpath(dir, tf.site * ".md"), tf; files = [edi, dat, png])
    return [edi, dat, png, report]
end

# The line below the plots. It tells what to do, in grey, until there is a new
# message. A message is cyan when the action was successful, red when it
# failed and navy while the processing runs. An empty message shows the
# instruction again
function _proc_status!(p::TKProc, text::AbstractString; error::Bool = false, busy::Bool = false)
    p.status.color[] = isempty(text) ? TK_GREY : error ? DASH_STATUS_ERROR : busy ? DASH_NAVY : DASH_STATUS_OK
    p.status.text[] = isempty(text) ? PROC_HINT : text
    return p
end

#---------- running -----

"""
    run_tkproc(site_dir; kwargs...) -> TKProc
    run_tkproc(survey::Survey, site; kwargs...) -> TKProc
    run_tkproc(; kwargs...) -> TKProc

Open the transfer function window and wait until you close it. Then return the
[`TKProc`](@ref) with the estimates that you made in it (`proc.results`).
With `site_dir`, the window opens on that site with its base and remote sites
(refer to [`TKProc`](@ref)). Without it, the window opens empty: use Load
Site. Start Julia with more than one thread (`julia -t auto`). Thus, the
window stays live while a site is processed.

```julia
proc = run_tkproc("data/survey/site004")   # process, export, close the window
export_tf("out", proc.results["site004"])
```
"""
function run_tkproc(p::TKProc)
    @info "Opening TKProc" site = p.focus == 0 ? "none" : p.survey.sites[p.focus].name threads = Threads.nthreads()
    screen = display(GLMakie.Screen(; title = "TKProc", focus_on_show = true), p.figure)
    _place_dash_window!(screen, _header_width(p.header))
    try
        wait(screen)
    catch
    end
    @info "TKProc closed" processed = length(p.results)
    return p
end

run_tkproc(survey::Survey, site; kwargs...) = run_tkproc(TKProc(survey, site; kwargs...))
run_tkproc(site_dir::AbstractString; kwargs...) = run_tkproc(TKProc(site_dir; kwargs...))
run_tkproc(; kwargs...) = run_tkproc(TKProc(; kwargs...))
