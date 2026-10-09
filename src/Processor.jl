# Processor.jl - the transfer function window, TKProc.
# Author: @pankajkmishra
#
# A GLMakie window that estimates the transfer function of one site at a time
# (Processing.jl) and writes it as EDI and a plot (EDI.jl).
# It has:
# - a header: Survey, the site menu, its base site, its remote site, the rate,
#   the FFT window, its overlap and the method
# - two rows of processing options, with the Full tensor, Error bars and
#   Phase 0-90° switches and the legend at their right
# - the apparent resistivity above the phase, and Tzx above Tzy, against the
#   period, with the error bars
# - the polarity check of the estimate (Polarity.jl)
# - a progress bar and a status line that follow the processing
#
# Survey opens a survey directory. The window scans it once and lists its
# sites with electric channels in the site menu. When you select a site, the
# window fills the base and remote menus with the sites that recorded with it,
# nearest and longest first. It selects the set-up of
# default_references: the best remote site, and for a site without Hx and Hy,
# the best base site.
#
# If you selected the base and remote sites of each site in TKDash and
# exported them, the window uses that plan: the reference_plan.txt in the
# survey directory or up to two directories above it. The menus then hold
# only the sites of the plan for the site, and the first base site and the
# first remote site are selected. Each Process uses one combination: one base
# site (or none) and one remote site (or none). Without a plan, or for a site
# that the plan does not list, the window makes the lists itself and the
# status line notes that a plan from TKDash can choose them. Process estimates the site in a task of its own. Thus,
# the window stays live and the status line shows each step. Below the plots,
# the check of the channels (Polarity.jl) tells if a channel looks flipped.
# After Process, it reports only what the quadrants of Zxy and Zyx show.
# FlipCheck then compares H of the inputs with the remote sites and up to two
# base sites (witnesses, which take no part in the estimate) to tell which
# channel is flipped.
# Export writes the EDI and a plot of the site. Their name holds the site, its
# base and remote sites and the values of the options (tf_filename)
#
# The plots: Zxy red, Zyx blue, Zxx green and Zyy lilac, circles with a black
# edge over grey error bars. Each phase is as recorded, on -200° to 200°, or
# wrapped to 0°-90°. The real tipper is red and the imaginary tipper is blue

const PROC_ZCOLOURS = (RGBf(0.62, 0.84, 0.60), RGBf(0.84, 0.30, 0.10),
                       RGBf(0.12, 0.38, 0.72), RGBf(0.78, 0.66, 0.88))     # xx, xy, yx, yy
const PROC_RE = RGBf(0.84, 0.30, 0.10)
const PROC_IM = RGBf(0.12, 0.38, 0.72)
const PROC_MARKER = (marker = :circle, markersize = 10, strokecolor = :black, strokewidth = 1.6)
const PROC_PHASE_RANGE = (-200.0, 200.0)
const PROC_LABEL = parse(RGBf, "#0f766e")        # teal: the selection in the site, base, remote and method menus
const PROC_TEXT = parse(RGBf, "#0f766e")         # teal: the lines below the plots
const PROC_SURVEY = parse.(RGBf, ("#7cc4bd", "#6ab8b0", "#58aca4"))   # light teal: the Survey… button, hover, pressed
const PROC_OK = parse(RGBf, "#0f766e")           # teal: the processing runs or ran as expected
const PROC_BAR = RGBf(0.74, 0.76, 0.79)        # light grey: the progress bar
const PROC_ERROR = parse(RGBf, "#d55e00")        # vermillion (Okabe-Ito): an action failed, apart from teal for colour-blind eyes
const PROC_HINT = "Survey… opens a survey · Site selects a site with its base and remote sites · Process estimates it · " *
                  "FlipCheck finds a flipped channel · Export writes its EDI and plot · " *
                  "Clear empties the screen."

# The text of the ⓘ next to each control
const PROC_HELP = Dict(
    :base => "The site that gives the magnetic inputs Hx, Hy. Use a base site for a site without " *
             "magnetic channels, or with poor ones. The first entry is the site itself: it gives its own Hx, Hy.",
    :remote => "The far site whose Hx, Hy are the reference. Noise in the local magnetic field is not " *
               "at the remote site. Thus, the estimate is free of the downward bias of a single site. " *
               "With none, the estimate is single site.",
    :rate => "The sampling rate to process. All rates estimates each rate that all the sites recorded " *
             "and keeps the estimate with the smallest error at each period. Bursts at a high rate " *
             "then widen the spectrum at the short periods.",
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
    :levels => "The most decimation levels. Each level is a factor of 4 below the one above it and adds " *
               "longer periods, if the record gives Min windows windows there. Higher: longer periods " *
               "from long records. Lower: the estimate stops at shorter periods.",
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
    :wrap => "Wrap the phases to 0°-90° as in conventional MT plots. The export plot follows this switch.",
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
- the loaded [`Survey`](@ref);
- the selected site;
- the [`TransferFunction`](@ref) of each site that you processed, in
  `results`;
- the GLMakie figure.

To make one, use `TKProc(survey_dir)` and `display` it, or use
[`run_tkproc`](@ref). After you close the window, `proc.results` holds the
estimates.
"""
mutable struct TKProc
    figure::Figure
    header::GridLayout
    survey::Survey
    focus::Int
    results::Dict{String, TransferFunction}
    site_menu::Menu
    base_menu::Menu
    remote_menu::Menu
    rate_menu::Menu
    boxes::Dict{Symbol, Textbox}
    method_menu::Menu
    leverage_toggle::Toggle
    full_toggle::Toggle
    bars_toggle::Toggle
    wrap_toggle::Toggle
    plan_source::Any
    plan::Vector{Any}
    plan_path::String
    process_button::Button
    flip_button::Button
    axes::Vector{Axis}
    status::Label
    polarity::Label
    bar::Makie.Poly
    bar_span::Observable{Tuple{Float64, Float64}}
    updating::Bool
    busy::Bool
end

"""
    TKProc(survey_dir, site = nothing; plan = nothing, size = (1600, 950)) -> TKProc
    TKProc(survey::Survey, site; plan = nothing, size = (1600, 950)) -> TKProc
    TKProc(; plan = nothing, size = (1600, 950)) -> TKProc

Make the transfer function window for a survey. With `survey_dir`, the
function scans the survey with [`scan_survey`](@ref) and lists its sites with
electric channels in the Site menu. It selects `site` (a name or an index), or
the first of these sites. With a [`Survey`](@ref), give the site by its name or
index. Each site that you select comes with its base and remote sites. Without
a survey, the window opens empty: use its Survey button. The function does not
open a window.

`plan` selects the base and remote sites of a TKDash plan
([`write_reference_plan`](@ref)):
- `nothing` -- the `reference_plan.txt` in the survey directory or up to two
  directories above it, if there is one;
- a path -- that plan file;
- `false` -- no plan: the window calculates the lists.
"""
function TKProc(survey_dir::AbstractString, site = nothing; kwargs...)
    p = TKProc(; kwargs...)
    _load_proc_survey!(p, survey_dir; site) || error(p.status.text[])
    return p
end

function TKProc(survey::Survey, site; kwargs...)
    p = TKProc(; kwargs...)
    _set_proc_survey!(p, survey)
    t = _target_site(survey, site)
    _proc_focus!(p, findfirst(s -> s === t, survey.sites))
    return p
end

function TKProc(; plan = nothing, size = (1600, 950))
    GLMakie.activate!(title = "TKProc")
    fig = Figure(; size = size, figure_padding = (16, 16, 10, 10))

    # 4 px between a label and its control, 14 px between the pairs
    pair_gaps!(grid, gaps) = foreach(((j, g),) -> colgap!(grid, j, g), enumerate(gaps))
    help = Tuple{Any, String}[]
    hint(label, key) = (push!(help, (label, _wrap(PROC_HELP[key]))); label)

    # the header: the site and its set-up on the left, the method and the
    # actions on the right
    header = GridLayout(fig[1, 1]; tellwidth = false)
    b_load = Button(header[1, 1]; label = "Survey…", font = :bold, labelcolor = :black, buttoncolor = PROC_SURVEY[1],
                    buttoncolor_hover = PROC_SURVEY[2], buttoncolor_active = PROC_SURVEY[3])
    site_menu = _bold_selection!(_dash_menu(header[1, 2]; options = [("no site", 0)], width = 130))
    hint(Label(header[1, 3], "Base ⓘ"), :base)
    base_menu = _bold_selection!(_dash_menu(header[1, 4]; options = [("none", nothing)], width = 110))
    hint(Label(header[1, 5], "Remote ⓘ"), :remote)
    remote_menu = _bold_selection!(_dash_menu(header[1, 6]; options = [("none", nothing)], width = 110))
    hint(Label(header[1, 7], "Rate ⓘ"), :rate)
    rate_menu = _dash_menu(header[1, 8]; options = [("All rates", :all)], width = 90)
    # the FFT window and its overlap, next to the rate (the other options are
    # in the rows below)
    isint(lo, hi) = s -> (v = tryparse(Int, s); v !== nothing && lo <= v <= hi)
    isreal(lo, hi) = s -> (v = tryparse(Float64, s); v !== nothing && lo <= v <= hi)
    hint(Label(header[1, 9], "Window ⓘ"), :window)
    window_box = Textbox(header[1, 10]; stored_string = "256", width = 52,
                         validator = s -> (v = tryparse(Int, s); v !== nothing && 64 <= v <= 65536 && ispow2(v)))
    hint(Label(header[1, 11], "Overlap ⓘ"), :overlap)
    overlap_box = Textbox(header[1, 12]; stored_string = "0.5", width = 52, validator = isreal(0.0, 0.9))
    Box(header[1, 13]; visible = false)
    hint(Label(header[1, 14], "Method ⓘ"), :method)
    method_menu = _bold_selection!(_dash_menu(header[1, 15]; options = [("CT2004", :ct2004), ("EB1986", :eb1986)], width = 90))
    b_run = Button(header[1, 16]; label = "Process")
    b_flip = Button(header[1, 17]; label = "FlipCheck")
    b_export = Button(header[1, 18]; label = "Export…")
    b_clear = Button(header[1, 19]; label = "Clear")
    colsize!(header, 13, Auto(true, 1.0))
    pair_gaps!(header, (8, 10, 4, 10, 4, 10, 4, 10, 4, 10, 4, 0, 0, 4, 10, 4, 4, 4))

    # the processing options, with the keywords of estimate_tf: the spectrum in
    # the first row, the robust fit in the second. The view switches are at the
    # right of the first row and the legend at the right of the second
    opts = GridLayout(fig[2, 1]; tellwidth = false)
    rows = [GridLayout(opts[1, 1]; halign = :left), GridLayout(opts[2, 1]; halign = :left)]
    boxes = Dict{Symbol, Textbox}(:window => window_box, :overlap => overlap_box)
    function field!(row, j, key, label, value, check)
        hint(Label(rows[row][1, 2j - 1], label * " ⓘ"), key)
        boxes[key] = Textbox(rows[row][1, 2j]; stored_string = value, width = 52, validator = check)
    end
    function switch!(g, j, label, key; active)
        hint(Label(g[1, 2j - 1], label * " ⓘ"), key)
        return Toggle(g[1, 2j]; active = active)
    end
    field!(1, 1, :bands, "Bands/decade", "8", isint(2, 20))
    field!(1, 2, :prewhiten, "AR order", "3", isint(0, 20))
    field!(1, 3, :top, "Top f ×Nyquist", "auto", s -> lowercase(strip(s)) == "auto" || isreal(0.1, 0.9)(s))
    field!(1, 4, :harmonic, "Lowest harmonic", "4", isint(1, 64))
    field!(1, 5, :windows, "Min windows", "8", isint(3, 1000))
    field!(1, 6, :levels, "Max levels", "12", isint(1, 20))
    field!(2, 1, :huber, "Huber", "1.5", isreal(0.5, 10.0))
    field!(2, 2, :jackknife, "Jackknife blocks", "50", isint(3, 1000))
    leverage_toggle = switch!(rows[2], 3, "Leverage weights", :leverage; active = true)
    pair_gaps!(rows[1], ntuple(j -> isodd(j) ? 4 : 14, 11))
    pair_gaps!(rows[2], (4, 14, 4, 14, 4))

    Box(opts[1:2, 2]; visible = false)
    view = GridLayout(opts[1, 3]; halign = :right)
    full_toggle = switch!(view, 1, "Full tensor", :full; active = false)
    bars_toggle = switch!(view, 2, "Error bars", :bars; active = true)
    wrap_toggle = switch!(view, 3, "Phase 0-90°", :wrap; active = false)
    pair_gaps!(view, (4, 14, 4, 14, 4))
    _proc_legend!(opts[2, 3]; halign = :right)
    colsize!(opts, 2, Auto(true, 1.0))
    rowgap!(opts, 6)

    # Outside: the axis labels stay in the column, so the controls above use
    # the full width of the window
    body = GridLayout(fig[3, 1]; alignmode = Outside())
    axes = _tf_axes!(body)
    polarity = Label(fig[4, 1], ""; fontsize = 12, halign = :left, tellwidth = false, justification = :left)
    # a translucent progress bar above the status line, as in the loading window
    bar_place = Box(fig[5, 1]; visible = false, height = 6)
    track = bar_place.layoutobservables.computedbbox
    poly!(fig.scene, track; color = PROGRESS_TRACK, strokewidth = 0, space = :pixel)
    bar_span = Observable((0.0, 0.0))
    bar = poly!(fig.scene, lift((bb, (x, w)) -> Rect2f(bb.origin[1] + x * bb.widths[1], bb.origin[2],
                                                        w * bb.widths[1], bb.widths[2]), track, bar_span);
                color = PROC_BAR, strokewidth = 0, space = :pixel, visible = false)
    status = Label(fig[6, 1], ""; fontsize = 13, halign = :left, tellwidth = false, color = PROC_TEXT)
    colsize!(fig.layout, 1, Relative(1))
    rowsize!(fig.layout, 3, Auto(true, 1.0))
    rowgap!(fig.layout, 8)

    p = TKProc(fig, header, Survey("", SurveySite[]), 0, Dict{String, TransferFunction}(), site_menu,
               base_menu, remote_menu, rate_menu, boxes, method_menu, leverage_toggle,
               full_toggle, bars_toggle, wrap_toggle, plan, Any[], "", b_run, b_flip, axes, status, polarity, bar, bar_span, false, false)

    on(b_load.clicks) do _
        p.busy && return
        dir = try
            pick_folder()
        catch err
            _proc_status!(p, "Could not open a folder dialog: $(sprint(showerror, err))"; error = true)
            ""
        end
        isempty(dir) || _load_proc_survey!(p, dir)
    end
    # a site from the menu comes with its base and remote sites. While a site
    # is processed, the menu goes back to the site in focus
    on(site_menu.selection) do i
        (p.updating || i === nothing || i == 0 || i == p.focus) && return
        p.busy ? _select_proc_site!(p) : _proc_focus!(p, i)
    end
    on(_ -> _process_focus!(p), b_run.clicks)
    on(_ -> _export_site!(p), b_export.clicks)
    on(_ -> _flip_check!(p), b_flip.clicks)
    on(_ -> _clear_proc!(p), b_clear.clicks)
    on(_ -> _draw_proc!(p), full_toggle.active)
    on(_ -> _draw_proc!(p), bars_toggle.active)
    on(_ -> _draw_proc!(p), wrap_toggle.active)
    _install_proc_help!(fig, help)
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
                _proc_status!(p, @sprintf("FlipCheck » %s · %.0f s", step, time() - t0); busy = true)
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
                             join(others, ", ") * @sprintf(" in %.1f s.", time() - t0); ok = true)
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
    p.polarity.color[] = c.ok ? PROC_TEXT : DASH_REMOTE_RAMP[end]
    p.polarity.text[] = (get(tf.metadata, :flipcheck, false) ? "FlipCheck: " : "Channels: ") * _polarity_summary(c)
    return p
end

# The selection of a site or method menu in bold teal; the list stays regular. Menu has
# no font, so this sets the text plot of its selection
function _bold_selection!(menu)
    for pl in menu.blockscene.plots
        Makie.plotfunc(pl) === Makie.editabletext || continue
        pl.font = :bold
        pl.color = PROC_LABEL
    end
    return menu
end

# The check of the channels in a few plain words; check_polarity has the
# full reasoning
function _polarity_summary(c)
    c.ok && return "all channels look fine."
    c.swapped && return "the x and y channels look swapped, or the layout is turned about 90°."
    c.parallel_e && return "Ex and Ey record the same signal: check the electric dipoles."
    c.parallel_h && return "Hx and Hy record the same signal: check the magnetic sensors."
    isempty(c.flipped) || return "reversed (wrong sign): $(join(c.flipped, ", ")). Swap the wires of " *
                                 (length(c.flipped) == 1 ? "this channel." : "these channels.")
    pairs = [z == :xy ? "Ex or Hy" : "Ey or Hx" for z in (:xy, :yx) if getfield(c, Symbol(:z, z)) === :flipped]
    isempty(pairs) || return "a channel is reversed ($(join(pairs, ", and "))). " *
                             "Run FlipCheck, or process with a remote site, to tell which."
    return "no clear answer: the data are too noisy, or the ground is 3D near the surface."
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
# White with no grid, as in MTGeophysics.jl
function _tf_axes!(grid)
    axis(pos; kw...) = Axis(pos; xscale = log10, backgroundcolor = :white, xgridvisible = false,
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

function _proc_legend!(pos; kwargs...)
    items = Any[MarkerElement(; PROC_MARKER..., color = PROC_ZCOLOURS[i]) for i in (2, 3, 1, 4)]
    append!(items, [MarkerElement(; PROC_MARKER..., color = PROC_RE), MarkerElement(; PROC_MARKER..., color = PROC_IM)])
    return Legend(pos, items, ["Zxy", "Zyx", "Zxx", "Zyy", "Re T", "Im T"]; orientation = :horizontal,
                  framevisible = false, labelsize = 12, patchsize = (14, 12), kwargs...)
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

# One series: grey error bars under circles with a black edge
function _tf_series!(ax, T, y, lo, hi, colour; bars::Bool = true)
    ok = findall(isfinite, y)
    isempty(ok) && return Float64[]
    e = filter(i -> isfinite(lo[i]) && isfinite(hi[i]), ok)
    (isempty(e) || !bars) || errorbars!(ax, T[e], y[e], lo[e], hi[e]; color = :grey40, linewidth = 1, whiskerwidth = 6)
    scatter!(ax, T[ok], y[ok]; PROC_MARKER..., color = colour)
    return y[ok]
end

# The phase axis: -200° to 200°, or 0° to 90° (180° if a wrapped phase is above 90°)
function _phase_axis!(ax, wrap::Bool, seen = Float64[])
    hi = maximum(seen; init = 0.0) > 90 ? 180.0 : 90.0
    ax.yticks = wrap ? (0:(hi > 90 ? 45 : 15):hi) : (-180:90:180)
    return wrap ? (0.0, hi) : PROC_PHASE_RANGE
end

# The impedance and the tipper of `tf` on the four axes. `full_tensor` adds
# Zxx and Zyy to Zxy and Zyx, `wrap_phase` wraps the phases to 0°-90°
function _draw_tf!(axes, tf::TransferFunction; full_tensor::Bool = false, bars::Bool = true,
                   wrap_phase::Bool = false)
    ax_rho, ax_phi, ax_tzx, ax_tzy = axes
    foreach(empty!, axes)
    T = tf.periods
    isempty(T) && return nothing
    σ = sqrt.(tf.Z_var ./ 2)
    seen, φseen = Float64[], Float64[]
    for (i, j, c) in (full_tensor ? ((1, 2, 2), (2, 1, 3), (1, 1, 1), (2, 2, 4)) : ((1, 2, 2), (2, 1, 3)))
        z = tf.Z[i, j, :]
        rho = [isfinite(v) && abs(v) > 0 ? 0.2 * t * abs2(v) : NaN for (v, t) in zip(z, T)]
        phi = [isfinite(r) ? rad2deg(angle(v)) : NaN for (v, r) in zip(z, rho)]
        wrap_phase && (phi = mod.(phi, 180.0))           # Zyx from the third quadrant to the first
        rel = σ[i, j, :] ./ abs.(z)
        δρ, δφ = 2 .* rho .* rel, min.(rad2deg.(rel), 90)
        append!(seen, _tf_series!(ax_rho, T, rho, min.(δρ, 0.95 .* rho), δρ, PROC_ZCOLOURS[c]; bars))
        append!(φseen, _tf_series!(ax_phi, T, phi, δφ, δφ, PROC_ZCOLOURS[c]; bars))
    end
    xl = (first(T) / 1.5, last(T) * 1.5)
    rl = isempty(seen) ? (1.0, 1.0e4) : (10^(log10(minimum(seen)) - 0.5), 10^(log10(maximum(seen)) + 0.5))
    isempty(seen) && _proc_note!(ax_rho, "no impedance at this site")
    ax_rho.yticks = _decade_ticks(rl...)
    ax_rho.limits[] = (xl, rl)
    ax_phi.limits[] = (xl, _phase_axis!(ax_phi, wrap_phase, φseen))
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
    plot_tf(tf::TransferFunction; full_tensor = false, errors = true, wrap_phase = false,
            size = (1300, 800)) -> Figure

A figure of `tf`: the apparent resistivity above the phase on the left, the
tipper Tzx above Tzy on the right, against the period, with the error bars.
Zxy is red and Zyx blue. With `full_tensor`, Zxx (green) and Zyy (lilac) are
there too. With `errors = false`, no error bars. With `wrap_phase`, the
phases are wrapped to 0°-90° (modulo 180°). The real tipper is red and
the imaginary tipper is blue. Save it
with `save("site.png", plot_tf(tf))`.
"""
function plot_tf(tf::TransferFunction; full_tensor::Bool = false, errors::Bool = true, wrap_phase::Bool = false,
                 size = (1300, 800))
    fig = Figure(; size = size, figure_padding = (16, 20, 10, 10))
    top = GridLayout(fig[1, 1]; tellwidth = false, halign = :left)
    Label(top[1, 1], isempty(tf.periods) ? tf.site : "$(tf.site) · $(_setup_text(tf))"; font = :bold,
          padding = (0, 20, 0, 0))
    _proc_legend!(top[1, 2])
    axes = _tf_axes!(GridLayout(fig[2, 1]))
    _draw_tf!(axes, tf; full_tensor, bars = errors, wrap_phase)
    return fig
end

function _draw_proc!(p::TKProc)
    name = p.focus == 0 ? "" : p.survey.sites[p.focus].name
    tf = get(p.results, name, nothing)
    _show_polarity!(p, tf)
    if tf === nothing
        foreach(empty!, p.axes)
        _proc_note!(p.axes[1], isempty(name) ? "no site loaded" : "$name: not processed")
        φl = _phase_axis!(p.axes[2], p.wrap_toggle.active[])
        for ax in p.axes
            ax.limits[] = ((1.0e-3, 1.0e3), ax === p.axes[1] ? (0.1, 1000.0) : ax === p.axes[2] ? φl : (-0.5, 0.5))
            ax.xticks = _decade_ticks(1.0e-3, 1.0e3)
            reset_limits!(ax)
        end
        p.axes[1].yticks = _decade_ticks(0.1, 1000.0)
        return p
    end
    _draw_tf!(p.axes, tf; full_tensor = p.full_toggle.active[], bars = p.bars_toggle.active[],
              wrap_phase = p.wrap_toggle.active[])
    return p
end

#---------- state changes -----

# Load the survey in `dir`: scan it, read its TKDash plan, list its sites with
# electric channels in the site menu and select `site` (a name or an index),
# or the first of them
function _load_proc_survey!(p::TKProc, dir::AbstractString; site = nothing)
    root = _norm_path(dir)
    if !isdir(root)
        _proc_status!(p, "$root is not a directory: give the directory of a survey."; error = true)
        return false
    end
    _proc_status!(p, "Scanning $root for its sites …"; busy = true)
    survey = try
        scan_survey(root)
    catch err
        _proc_status!(p, "Could not scan $root: $(sprint(showerror, err))"; error = true)
        return false
    end
    if length(survey.sites) == 1 && survey.sites[1].path == root
        _proc_status!(p, "$root is a site, not a survey: give the survey directory $(dirname(root)) " *
                         "and select $(survey.sites[1].name) in the Site menu."; error = true)
        return false
    end
    targets = _proc_targets(survey)
    if isempty(targets)
        _proc_status!(p, "$root has no site with electric channels: give a directory with Metronix, " *
                         "LEMI-424 or GEOMAG recordings."; error = true)
        return false
    end
    i = site === nothing ? first(targets) : _proc_site_index(survey, site)
    if i === nothing || !(i in targets)
        _proc_status!(p, "$site is not a site with electric channels in $root."; error = true)
        return false
    end
    _set_proc_survey!(p, survey)
    _proc_focus!(p, i)
    return true
end

# The index of a site (a name or an index) in the survey, or nothing
function _proc_site_index(survey::Survey, site)
    t = try
        _target_site(survey, site)
    catch
        return nothing
    end
    return findfirst(s -> s === t, survey.sites)
end

# The sites that can be processed: the ones with an electric channel. The
# others are base or remote sites only
_proc_targets(survey::Survey) =
    [i for (i, s) in enumerate(survey.sites) if any(c -> c in site_components(s), (:e1, :e2))]

# Make `survey` the survey of the window: its plan and its site menu
function _set_proc_survey!(p::TKProc, survey::Survey)
    p.survey, p.focus = survey, 0
    _load_plan!(p, survey.root)
    opts = Any[(has_magnetic(survey.sites[i]) ? survey.sites[i].name : survey.sites[i].name * " · E only", i)
               for i in _proc_targets(survey)]
    p.updating = true
    try
        p.site_menu.options[] = isempty(opts) ? Any[("no site", 0)] : opts
        p.site_menu.i_selected[] = 1
    finally
        p.updating = false
    end
    return p
end

# Show the site in focus in the site menu
function _select_proc_site!(p::TKProc)
    k = findfirst(o -> o[2] == p.focus, p.site_menu.options[])
    k === nothing && return p
    p.updating = true
    try
        p.site_menu.i_selected[] = k
    finally
        p.updating = false
    end
    return p
end

# The TKDash plan: the file that `plan_source` names, or the first
# reference_plan.txt in `root` and up to two directories above it
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
    return p
end

# The plan entry of a site, or nothing if there is no plan or the plan does
# not list the site
function _plan_entry(p::TKProc, name)
    isempty(p.plan) && return nothing
    k = findfirst(r -> r.site == name, p.plan)
    return k === nothing ? nothing : p.plan[k]
end

function _proc_focus!(p::TKProc, i::Integer)
    p.busy && return p
    p.focus = i
    site = p.survey.sites[i]
    entry = _plan_entry(p, site.name)
    # a site with Hx, Hy is its own base site: the first entry of the base
    # menu is the site, and it takes the inputs from its own Hx, Hy
    local_base = (has_magnetic(site) ? site.name : "none", nothing)
    if entry === nothing
        refs = site_references(p.survey, site)
        defaults = default_references(p.survey, site)
        base_opts = vcat(Any[local_base], Any[(c.site, c.site) for c in refs.base])
        remote_opts = vcat(Any[("none", nothing)], Any[(c.site, c.site) for c in refs.remote])
        length(refs.remote) > 1 && push!(remote_opts, ("all remotes", :all))
        nb, nr = length(refs.base), length(refs.remote)
        source = "recorded with it"
    else
        # one combination at a time: the menus hold the sites of the plan, no
        # "all remotes"
        label(name, hours) = @sprintf("%s · %.1f h", name, hours)
        base_opts = vcat(Any[local_base], Any[(label(n, h), n) for (n, h) in zip(entry.base, entry.base_hours)])
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
        p.base_menu.options[] = base_opts
        p.base_menu.i_selected[] = pick(base_opts, defaults.base)
        p.remote_menu.options[] = remote_opts
        p.remote_menu.i_selected[] = pick(remote_opts, defaults.remote)
        p.rate_menu.options[] = rate_opts
        p.rate_menu.i_selected[] = 1
    finally
        p.updating = false
    end
    _select_proc_site!(p)
    _draw_proc!(p)
    loaded = "Loaded $(site.name): $nb base site$(nb == 1 ? "" : "s") and $nr remote site$(nr == 1 ? "" : "s") $source."
    # a note, not an error: without a plan for the site, the lists are the
    # ones that TKProc makes
    note = p.plan_source === false ? "" :
           isempty(p.plan) ? " No TKDash plan found, so TKProc made the lists of base and remote sites." :
           entry === nothing ? " The TKDash plan does not list $(site.name), so TKProc made its lists." : ""
    isempty(note) || (note *= " To choose them yourself, make a plan in TKDash (Export…).")
    _proc_status!(p, haskey(p.results, site.name) ? "Showing the estimate of $(site.name)." :
                     loaded * note * " Select one base and one remote site, then Process.")
    return p
end

# The keywords of estimate_tf from the option row
function _proc_options(p::TKProc)
    get(key, T, default) = something(tryparse(T, p.boxes[key].stored_string[]), default)
    return (; window = get(:window, Int, 256), overlap = get(:overlap, Float64, 0.5),
            bands_per_decade = get(:bands, Int, 8), prewhiten = get(:prewhiten, Int, 3),
            huber = get(:huber, Float64, 1.5), jackknife_groups = get(:jackknife, Int, 50),
            nyquist_fraction = something(tryparse(Float64, p.boxes[:top].stored_string[]), :auto), min_harmonic = get(:harmonic, Int, 4),
            min_windows = get(:windows, Int, 8), max_levels = get(:levels, Int, 12),
            min_period = 0.0, max_period = Inf,
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
    try
        opts = _proc_options(p)
        task = Threads.@spawn estimate_tf(p.survey, name; base, remote, rate, opts...,
                                          progress = msg -> (put!(messages, msg); yield()))
        step = "starting"
        while !istaskdone(task)
            while isready(messages)
                step = take!(messages)
            end
            _proc_status!(p, @sprintf("Processing » %s · %.0f s", step, time() - t0); busy = true)
            sleep(0.25)
        end
        tf = try
            fetch(task)
        catch err
            throw(err isa TaskFailedException ? err.task.exception : err)
        end
        p.results[name] = tf
        _draw_proc!(p)
        _proc_status!(p, @sprintf("Done » %d periods, %.3g–%.3g s, in %.1f s.", length(tf.periods),
                                  first(tf.periods), last(tf.periods), time() - t0); ok = true)
    catch err
        _proc_status!(p, "Processing failed » $(sprint(showerror, err))"; error = true)
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

# Export the site in focus: its EDI and a PNG of the view on the screen, into
# a directory that you select
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
        written = export_tf(dir, tf; full_tensor = p.full_toggle.active[], errors = p.bars_toggle.active[],
                            wrap_phase = p.wrap_toggle.active[])
        _proc_status!(p, "Wrote " * join(basename.(written), ", ") * " to $dir")
    catch err
        _proc_status!(p, "Could not export $name: $(sprint(showerror, err))"; error = true)
    end
    return p
end

"""
    tf_filename(tf::TransferFunction) -> String

The name (without extension) of the export of `tf`: the site, the base site,
the remote sites, the rate and the values of the options, joined by `-`:

    site-base-remote-rate-window-overlap-bands-prewhiten-top-harmonic-windows-levels-method-huber-leverage-jackknife

e.g. `site002-site002-site099-all-256-0.5-8-3-auto-4-8-12-ct2004-1.5-1-50`.
The base is the site itself when the site gives its own Hx, Hy; the remote
is `none` without a remote site, and remote sites join with `+`. `rate` is
`all` or the rate that was processed (`128Hz`, `16s`); `top` is `auto` or the
fraction of Nyquist; `leverage` is 1 or 0. The name holds all you need to
give [`estimate_tf`](@ref) to make the estimate again. A period range other
than all periods adds `minperiod-maxperiod` at the end.
"""
function tf_filename(tf::TransferFunction)
    md = tf.metadata
    clean(x) = replace(string(x), r"[/\\:*?\"<>| ]" => "_")
    num(x) = isinteger(x) ? string(Int(x)) : string(x)
    rate = get(md, :rate, :all)
    parts = String[tf.site, isempty(tf.base) ? tf.site : tf.base,
                   isempty(tf.remote) ? "none" : join(tf.remote, "+"),
                   rate === :all ? "all" : replace(_fs_label(rate), " " => "")]
    if haskey(md, :options)
        o = md[:options]
        append!(parts, [string(o.window), num(o.overlap), string(o.bands_per_decade), string(o.prewhiten),
                        o.nyquist_fraction === :auto ? "auto" : num(o.nyquist_fraction), string(o.min_harmonic),
                        string(o.min_windows), string(o.max_levels), string(o.method), num(o.huber),
                        o.leverage ? "1" : "0", string(o.jackknife_groups)])
        (o.min_period > 0 || isfinite(o.max_period)) && append!(parts, [num(o.min_period), num(o.max_period)])
    end
    return join(clean.(parts), "-")
end

"""
    export_tf(dir, tf::TransferFunction; full_tensor = false, errors = true, wrap_phase = false) -> Vector{String}

Write the estimate of one site into `dir` as `<name>.edi` ([`write_edi`](@ref))
and `<name>.png` ([`plot_tf`](@ref), with or without Zxx and Zyy and the
error bars). `<name>` is [`tf_filename`](@ref): the site, its base and remote
sites and the values of the options. The INFO block of the EDI file holds
each sensor with its calibration file or its dipole length and where the
length came from. The errors are the estimated errors, with no floor. The
function returns the paths.
"""
function export_tf(dir::AbstractString, tf::TransferFunction; full_tensor::Bool = false, errors::Bool = true,
                   wrap_phase::Bool = false)
    mkpath(dir)
    name = tf_filename(tf)
    edi = write_edi(joinpath(dir, name * ".edi"), tf)
    png = joinpath(dir, name * ".png")
    save(png, plot_tf(tf; full_tensor, errors, wrap_phase); px_per_unit = 2)
    return [edi, png]
end

# The line below the plots and the bar above it. The line tells what to do
# until there is a new message. A message is teal while a job runs (`busy`,
# the light grey bar sweeps) and when it ran as expected (`ok`, the bar
# fills), vermillion when an action failed (the bar fills vermillion). An empty
# message shows the instruction again
function _proc_status!(p::TKProc, text::AbstractString; error::Bool = false, busy::Bool = false, ok::Bool = false)
    p.status.color[] = error ? PROC_ERROR : (busy || ok) ? PROC_OK : PROC_TEXT
    p.status.text[] = isempty(text) ? PROC_HINT : text
    x, w = busy ? _progress_bar_span(:busy, 0, 0, time()) : (0.0, 1.0)
    p.bar.visible[] = busy || ok || error
    p.bar.color[] = error ? RGBAf(PROC_ERROR, 0.55) : PROC_BAR
    p.bar_span[] = (x, w)
    return p
end

#---------- running -----

"""
    run_tkproc(survey_dir, site = nothing; kwargs...) -> TKProc
    run_tkproc(survey::Survey, site; kwargs...) -> TKProc
    run_tkproc(; kwargs...) -> TKProc

Open the transfer function window and wait until you close it. Then return the
[`TKProc`](@ref) with the estimates that you made in it (`proc.results`).
With `survey_dir`, the window opens on that survey with `site` (or its first
site) and its base and remote sites selected (refer to [`TKProc`](@ref)).
Select the other sites in the Site menu. Without it, the window opens empty:
use Survey…. Start Julia with more than one thread (`julia -t auto`). Thus, the
window stays live while a site is processed.

```julia
proc = run_tkproc("data/survey", "site004")   # process, export, close the window
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
run_tkproc(survey_dir::AbstractString, site = nothing; kwargs...) = run_tkproc(TKProc(survey_dir, site; kwargs...))
run_tkproc(; kwargs...) = run_tkproc(TKProc(; kwargs...))
