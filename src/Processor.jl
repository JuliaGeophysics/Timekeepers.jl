# Processor.jl - the transfer function window, TKProc.
# Author: @pankajkmishra
#
# A GLMakie window that estimates the transfer function of one site at a time
# (Processing.jl) and writes it as EDI and a plot (EDI.jl).
# It has:
# - a header: Survey, the site menu, its base site, its remote site, the rate,
#   the FFT window, its overlap and the method
# - two rows of processing options, with the Full tensor, Error bars and
#   Phase 0-90° switches and the legend at their right, and Save
# - the set-up of each site in procset.dat in its directory: the window
#   fills the menus and the options from it when you select the site, and
#   Save writes it (read_procset, write_procset). The minimum distance of
#   remote sites, the period range, the time span and the rotation of the
#   plots are in procset.dat only
# - the apparent resistivity above the phase, and Tzx above Tzy, against the
#   period, with the error bars
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
# the window stays live and the status line shows each step. To reverse a
# channel, set -1 for it in d/orientations.dat of the survey.
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
const PROC_HINT = "Survey… opens a survey · Site selects a site with its base and remote sites · Process estimates it and writes its set-up to procset.dat of the site · " *
                  "Save set-up writes the set-up without processing · Export writes its EDI, figures and notes to proc/ of the site · " *
                  "Clear empties the screen."

# The text of the ⓘ next to each control
const PROC_HELP = Dict(
    :base => "The site that gives the magnetic inputs Hx, Hy. Use a base site for a site without " *
             "magnetic channels, or with poor ones. The first entry is the site itself: it gives its own Hx, Hy.",
    :remote => "The far site whose Hx, Hy are the reference. Noise in the local magnetic field is not " *
               "at the remote site. Thus, the estimate is free of the downward bias of a single site. " *
               "With none, the estimate is single site. The menu lists every site with Hx, Hy that recorded " *
               "at the same time, nearest first, with its distance. \"all ≥ …\" uses all those at the minimum " *
               "distance (min_distance in procset.dat, 100 m if none) or farther that recorded 75 % or more of the time of the site (a remote site that " *
               "recorded less would cut the data of the others; select it by itself).",
    :rate => "The sampling rate to process. All rates estimates each rate that all the sites recorded " *
             "and keeps the estimate with the smallest error at each period. Bursts at a high rate " *
             "then widen the spectrum at the short periods.",
    :section_length => "Section length: the samples in each data section that is tapered and Fourier transformed, as " *
                       "in the Welch overlapped section averaging (WOSA) of Chave & Thomson (2004, section 2). A " *
                       "power of 2. Chave & Thomson shorten the section length as the frequency rises; Timekeepers " *
                       "keeps it fixed and decimates the record by 4 at each level, which does the same. Longer: " *
                       "finer frequency resolution, but fewer sections. 256 suits most records.",
    :section_overlap => "The fraction of each data section that the next one shares. Chave & Thomson (2004, section " *
                        "2): data sections may be overlapped to improve statistical efficiency. Timekeepers tapers " *
                        "with a Hann window, not their Slepian sequence; 0.5 suits it.",
    :bands_per_decade => "Frequency bands per decade of the output. Each band averages the Fourier transforms of the " *
                         "frequencies in it. Chave & Thomson (2004, section 3) process single frequencies and avoid " *
                         "band averaging; Timekeepers averages bands. Higher: more periods with fewer data each, so " *
                         "larger errors.",
    :ar_terms => "The terms of the AR (autoregressive) pre-whitening filter of each time series, fitted by " *
                 "least squares (the Yule-Walker equations). Chave & Thomson (2004, section 2) use a robust AR " *
                 "filter instead; their examples have five terms. The pre-whitening is corrected for after the " *
                 "Fourier transform. 0: no pre-whitening.",
    :top_nyquist => "The highest frequency, as a fraction of the Nyquist frequency of the record. auto: where the " *
                    "anti-alias filter of the recorder ends the signal, found from the coherence of E with H (it " *
                    "falls below 0.5 there), less 5 %, from 0.5 to 0.9. A number sets it: higher gives shorter " *
                    "periods, but closer to the filter, where the data can be biased.",
    :min_frequency_index => "The lowest frequency index of a data section that the last decimation level uses. Chave & " *
                            "Thomson (2004, section 3) keep the frequency of interest a few over the section length " *
                            "(their second example uses the fourth and sixth frequencies of each section). Lower: longer " *
                            "periods, with more leakage from unresolved low-frequency components.",
    :min_sections => "The data sections that a decimation level needs to give an estimate. Lower: longer periods " *
                     "from fewer sections, so larger errors and less reliable robust weights.",
    :decimation_levels => "The most decimation levels. Each level is a factor of 4 below the one above it: it gives the " *
                          "longer periods that Chave & Thomson (2004) reach with longer sections.",
    :huber_a => "The parameter a of the Huber weights (Chave & Thomson 2004, eq. 7): downweighting begins " *
                "when a residual over the scale (its median magnitude over the median of a Rayleigh model) exceeds a. " *
                "a = 1.5 gives better than 95 per cent efficiency with outlier-free Gaussian data. The final " *
                "estimate then uses the more severe weight of their eq. (10).",
    :jackknife_blocks => "The jackknife confidence limits of Chave & Thomson (2004, section 3): the estimate is made " *
                         "again with one block of data sections deleted at a time. Timekeepers deletes blocks of " *
                         "sections, not single sections; this is their number.",
    :method => "CT2004: the two-stage bounded influence estimator of Chave & Thomson (2004, section 7). With " *
               "a remote reference, the first stage is the bounded influence regression of the local " *
               "magnetic field on the reference magnetic field; the second stage is the bounded influence " *
               "regression of the electric field on the predicted magnetic field, with the first-stage " *
               "weights kept. Without a remote reference it is the bounded influence estimator of their " *
               "section 6. EB1986: the robust remote reference M-estimate of Egbert & Booker (1986), with " *
               "Huber weights only.",
    :leverage_weights => "On: the bounded influence (BI) estimator of Chave & Thomson (2004, section 6), whose weights " *
                         "are the M-estimator weights times the leverage weights of their eq. (26), from the hat " *
                         "matrix diagonal of eq. (27). It limits the influence of leverage points: unusual magnetic " *
                         "field data. Off: the M-estimator only, the ordinary robust (OR) estimate. CT2004 method " *
                         "only.",
    :full => "Show Zxx and Zyy too, not only Zxy and Zyx. The export plot follows this switch.",
    :bars => "Show the error bars (one standard deviation). The export plot follows this switch.",
    :wrap => "Wrap the phases to 0°-90° as in conventional MT plots. The export plot follows this switch.",
    :save => "Write the set-up on the screen (base, remote, rate, method and all the options) to procset.dat " *
             "in the directory of the site, with its min_distance, min_period, max_period, span and rotate. " *
             "Process writes it too. The window reads it each time you select the site and when you press Process.",
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
    axes::Vector{Axis}
    status::Label
    bar::Makie.Poly
    bar_span::Observable{Tuple{Float64, Float64}}
    updating::Bool
    busy::Bool
    orient_label::Label
    defaults::Dict{Symbol, String}
    file_only::Dict{Symbol, String}
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
    remote_menu = _bold_selection!(_dash_menu(header[1, 6]; options = [("none", nothing)], width = 130))
    isint(lo, hi) = s -> (v = tryparse(Int, s); v !== nothing && lo <= v <= hi)
    isreal(lo, hi) = s -> (v = tryparse(Float64, s); v !== nothing && lo <= v <= hi)
    hint(Label(header[1, 7], "Rate ⓘ"), :rate)
    rate_menu = _dash_menu(header[1, 8]; options = [("All rates", :all)], width = 90)
    hint(Label(header[1, 9], "Section length ⓘ"), :section_length)
    window_box = Textbox(header[1, 10]; stored_string = "256", width = 52,
                         validator = s -> (v = tryparse(Int, s); v !== nothing && 64 <= v <= 65536 && ispow2(v)))
    hint(Label(header[1, 11], "Section overlap ⓘ"), :section_overlap)
    overlap_box = Textbox(header[1, 12]; stored_string = "0.5", width = 52, validator = isreal(0.0, 0.9))
    Box(header[1, 13]; visible = false)
    hint(Label(header[1, 14], "Method ⓘ"), :method)
    method_menu = _bold_selection!(_dash_menu(header[1, 15]; options = [("CT2004", :ct2004), ("EB1986", :eb1986)], width = 90))
    b_run = Button(header[1, 16]; label = "Process")
    b_export = Button(header[1, 17]; label = "Export…")
    b_clear = Button(header[1, 18]; label = "Clear")
    colsize!(header, 13, Auto(true, 1.0))
    pair_gaps!(header, (8, 10, 4, 10, 4, 10, 4, 10, 4, 10, 4, 0, 0, 4, 10, 4, 4))

    # the processing options, with the keywords of estimate_tf: the spectrum in
    # the first row, the robust fit in the second. The view switches are at the
    # right of the first row and the legend at the right of the second
    opts = GridLayout(fig[2, 1]; tellwidth = false)
    rows = [GridLayout(opts[1, 1]; halign = :left), GridLayout(opts[2, 1]; halign = :left)]
    boxes = Dict{Symbol, Textbox}(:section_length => window_box, :section_overlap => overlap_box)
    function field!(row, j, key, label, value, check)
        hint(Label(rows[row][1, 2j - 1], label * " ⓘ"), key)
        boxes[key] = Textbox(rows[row][1, 2j]; stored_string = value, width = 52, validator = check)
    end
    function switch!(g, j, label, key; active)
        hint(Label(g[1, 2j - 1], label * " ⓘ"), key)
        return Toggle(g[1, 2j]; active = active)
    end
    field!(1, 1, :bands_per_decade, "Bands/decade", "8", isint(2, 20))
    field!(1, 2, :ar_terms, "AR filter terms", "3", isint(0, 20))
    field!(1, 3, :top_nyquist, "Top f ×Nyquist", "auto", s -> lowercase(strip(s)) == "auto" || isreal(0.1, 0.9)(s))
    field!(1, 4, :min_frequency_index, "Lowest frequency index", "4", isint(1, 64))
    field!(2, 1, :min_sections, "Min sections", "8", isint(3, 1000))
    field!(2, 2, :decimation_levels, "Decimation levels", "12", isint(1, 20))
    field!(2, 3, :huber_a, "Huber a", "1.5", isreal(0.5, 10.0))
    field!(2, 4, :jackknife_blocks, "Jackknife blocks", "50", isint(3, 1000))
    leverage_toggle = switch!(rows[2], 5, "Leverage weights", :leverage_weights; active = true)
    b_save = hint(Button(rows[2][1, 11]; label = "Save set-up"), :save)
    pair_gaps!(rows[1], ntuple(j -> isodd(j) ? 4 : 14, 7))
    pair_gaps!(rows[2], (4, 14, 4, 14, 4, 14, 4, 14, 4, 14))

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
    # a translucent progress bar above the status line, as in the loading window
    bar_place = Box(fig[4, 1]; visible = false, height = 6)
    track = bar_place.layoutobservables.computedbbox
    poly!(fig.scene, track; color = PROGRESS_TRACK, strokewidth = 0, space = :pixel)
    bar_span = Observable((0.0, 0.0))
    bar = poly!(fig.scene, lift((bb, (x, w)) -> Rect2f(bb.origin[1] + x * bb.widths[1], bb.origin[2],
                                                        w * bb.widths[1], bb.widths[2]), track, bar_span);
                color = PROC_BAR, strokewidth = 0, space = :pixel, visible = false)
    bottom = GridLayout(fig[5, 1]; tellwidth = false)
    status = Label(bottom[1, 1], ""; fontsize = 13, halign = :left, tellwidth = false, color = PROC_TEXT)
    orient_label = Label(bottom[1, 2], ""; fontsize = 13, halign = :right, color = PROC_TEXT)
    colsize!(bottom, 1, Auto(true, 1.0))
    colgap!(bottom, 20)
    colsize!(fig.layout, 1, Relative(1))
    rowsize!(fig.layout, 3, Auto(true, 1.0))
    rowgap!(fig.layout, 8)

    p = TKProc(fig, header, Survey("", SurveySite[]), 0, Dict{String, TransferFunction}(), site_menu,
               base_menu, remote_menu, rate_menu, boxes, method_menu, leverage_toggle,
               full_toggle, bars_toggle, wrap_toggle, plan, Any[], "", b_run, axes, status, bar, bar_span,
               false, false, orient_label,
               merge(Dict(k => String(b.stored_string[]) for (k, b) in boxes), Dict(:method => "ct2004", :leverage_weights => "on")),
               copy(PROC_FILE_ONLY))
    on(_ -> _save_procset!(p), b_save.clicks)

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
# White with no grid
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
    _draw_tf!(p.axes, _proc_rotated(p, tf); full_tensor = p.full_toggle.active[], bars = p.bars_toggle.active[],
              wrap_phase = p.wrap_toggle.active[])
    return p
end

_proc_rotation(p::TKProc) = something(tryparse(Float64, strip(get(p.file_only, :rotate, "0"))), 0.0)

function _proc_rotated(p::TKProc, tf::TransferFunction)
    θ = _proc_rotation(p)
    return θ == 0 ? tf : rotate_tf(tf, θ)
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
    _read_file_only!(p, site)
    entry = _plan_entry(p, site.name)
    # a site with Hx, Hy is its own base site: the first entry of the base
    # menu is the site, and it takes the inputs from its own Hx, Hy
    local_base = (has_magnetic(site) ? site.name : "none", nothing)
    # every site with Hx, Hy that recorded with this one is a possible remote
    # site, nearest first; "all ≥ …" takes those at the minimum distance or
    # farther, and is the default when there is one
    concurrent = _concurrent_remotes(p.survey, site)
    far = _far_remotes(p, site)
    if entry === nothing
        refs = site_references(p.survey, site)
        defaults = (base = default_references(p.survey, site).base, remote = isempty(far) ? nothing : :all)
        base_opts = vcat(Any[local_base], Any[(c.site, c.site) for c in refs.base])
        remote_opts = vcat(Any[("none", nothing), (_all_remotes_label(p), :all)],
                           Any[(_distance_label(c.site, c.distance_km), c.site) for c in concurrent])
        nb, nr = length(refs.base), length(concurrent)
        source = "recorded with it"
    else
        # the sites of the plan, with "all ≥ …" of the sites that recorded with it
        label(name, hours) = @sprintf("%s · %.1f h", name, hours)
        base_opts = vcat(Any[local_base], Any[(label(n, h), n) for (n, h) in zip(entry.base, entry.base_hours)])
        remote_opts = vcat(Any[("none", nothing), (_all_remotes_label(p), :all)],
                           Any[(label(n, h), n) for (n, h) in zip(entry.remote, entry.remote_hours)])
        plan_base = isempty(entry.base) ? nothing : first(entry.base)
        plan_remotes = filter(!=(plan_base), entry.remote)
        defaults = (base = plan_base, remote = isempty(plan_remotes) ? nothing : first(plan_remotes))
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
    saved = _load_procset!(p, site)
    _select_proc_site!(p)
    _draw_proc!(p)
    _show_site_tables!(p, saved)
    loaded = "Loaded $(site.name): $nb base site$(nb == 1 ? "" : "s") and $nr remote site$(nr == 1 ? "" : "s") $source."
    # a note, not an error: without a plan for the site, the lists are the
    # ones that TKProc makes
    note = p.plan_source === false ? "" :
           isempty(p.plan) ? " No TKDash plan found, so TKProc made the lists of base and remote sites." :
           entry === nothing ? " The TKDash plan does not list $(site.name), so TKProc made its lists." : ""
    isempty(note) || (note *= " To choose them yourself, make a plan in TKDash (Export…).")
    if saved.error !== nothing
        _proc_status!(p, "Could not read $(saved.file): $(saved.error)"; error = true)
    elseif haskey(p.results, site.name)
        _proc_status!(p, "Showing the estimate of $(site.name).")
    elseif !isempty(saved.file)
        _proc_status!(p, "Loaded $(site.name) with its set-up from $(saved.file)" *
                         (isempty(saved.missing) ? "" : " ($(join(saved.missing, ", ")) not in the menus now)") *
                         ". Press Process.")
    else
        _proc_status!(p, loaded * note * " Select the base and remote sites, then Process; Save set-up keeps them.")
    end
    return p
end

_proc_min_km(p::TKProc) = something(tryparse(Float64, strip(get(p.file_only, :min_distance, "100"))), 100.0) / 1000

# The sites with Hx, Hy that recorded with `site` at the minimum distance or
# farther, nearest first. A site without a position is not among them
_far_remotes(p::TKProc, site::SurveySite) = remote_sites(p.survey, site; min_distance = 1000 * _proc_min_km(p))

_all_remotes_label(p::TKProc) = @sprintf("all ≥ %g m", round(1000 * _proc_min_km(p); sigdigits = 4))
_distance_label(name, km) = km < 1 ? @sprintf("%s · %.0f m", name, 1000km) : @sprintf("%s · %.1f km", name, km)

function _relabel_all_remotes!(p::TKProc)
    opts = p.remote_menu.options[]
    k = findfirst(o -> o[2] === :all, opts)
    k === nothing && return p
    sel = p.remote_menu.i_selected[]
    new = copy(opts)
    new[k] = (_all_remotes_label(p), :all)
    was = p.updating
    p.updating = true
    try
        p.remote_menu.options[] = new
        p.remote_menu.i_selected[] = sel
    finally
        p.updating = was
    end
    return p
end

const PROC_FILE_ONLY = Dict(:min_distance => "100", :min_period => "0", :max_period => "inf", :span => "all", :rotate => "0")

function _read_file_only!(p::TKProc, site::SurveySite)
    file = procset_file(site.path)
    saved = isfile(file) ? (try read_procset(file) catch; Dict{Symbol, String}() end) : Dict{Symbol, String}()
    p.file_only = Dict(k => get(saved, k, d) for (k, d) in PROC_FILE_ONLY)
    return p
end

#---------- the set-up of a site: procset.dat -----

# The keys of procset.dat in the order of the file, with the comment of each
const PROCSET_KEYS = (
    (:base, "the site that gives Hx, Hy; - = the site itself"),
    (:remote, "the reference site; - = none, all = all sites at min_distance or farther"),
    (:min_distance, "metres: remote all takes the sites this far or farther"),
    (:rate, "all, or a rate in Hz"),
    (:section_length, "samples in each data section (Chave & Thomson 2004)"),
    (:section_overlap, "fraction of a section that the next one shares"),
    (:method, "ct2004 (two-stage bounded influence) or eb1986"),
    (:bands_per_decade, "frequency bands per decade"),
    (:ar_terms, "terms of the AR pre-whitening filter"),
    (:top_nyquist, "highest frequency × Nyquist, or auto"),
    (:min_frequency_index, "lowest frequency index of a section, last level"),
    (:min_sections, "data sections that a decimation level needs"),
    (:decimation_levels, "most decimation levels"),
    (:huber_a, "a of the Huber weights (eq. 7)"),
    (:jackknife_blocks, "blocks of sections that the jackknife deletes"),
    (:leverage_weights, "on (bounded influence, eq. 26) or off"),
    (:min_period, "shortest period (s)"),
    (:max_period, "longest period (s), or inf"),
    (:span, "time span UTC, or all"),
    (:rotate, "rotation of the plots, degrees east of north"),
)

"""
    procset_file(site_dir) -> String

The set-up file of a site: `procset.dat` in its directory. For the card folder
of a logger (`DATA0247`) or a `corrected` copy, it is in the folder of the site
above it. The file need not be there.
"""
function procset_file(site_dir::AbstractString)
    d = rstrip(String(site_dir), ['/', '\\'])
    _is_card_dir(basename(d)) && (d = dirname(d))
    return joinpath(d, "procset.dat")
end

"""
    read_procset(path) -> Dict{Symbol, String}
    write_procset(path, settings; site = "")

The processing set-up of one site, as TKProc saves it (Save set-up) and reads
it when you select the site. One `key value` line each; `#` starts a comment:

    base        -
    remote      10
    section_length  256
    method          eb1986
    bands_per_decade 6
    min_period  8
    max_period  3000
    span        2026-10-20 14:00 to 2026-10-21 07:20

The keys: base, remote, min_distance, rate, section_length, section_overlap,
method, bands_per_decade, ar_terms, top_nyquist, min_frequency_index,
min_sections, decimation_levels, huber_a, jackknife_blocks, leverage_weights,
min_period, max_period, span and rotate. Their names are the terms of
Chave & Thomson (2004). A key that the file leaves out keeps
the default of a new TKProc window, and so does every key of a site without
the file. `write_procset` writes the keys of `settings` in that order.
"""
function read_procset(path::AbstractString)
    known = Set(first.(PROCSET_KEYS))
    out = Dict{Symbol, String}()
    for (i, raw) in enumerate(readlines(path))
        line = strip(replace(raw, r"#.*$" => ""))
        isempty(line) && continue
        m = match(r"^(\S+)\s+(.+)$", line)
        m === nothing && error("$path, line $i: write \"key value\", e.g. \"section_length 256\"")
        key = Symbol(lowercase(m.captures[1]))
        key in known || error("$path, line $i: no key $(m.captures[1]); the keys are $(join(first.(PROCSET_KEYS), ", "))")
        out[key] = String(strip(m.captures[2]))
    end
    return out
end

function write_procset(path::AbstractString, settings::AbstractDict; site::AbstractString = "")
    mkpath(dirname(abspath(path)))
    open(path, "w") do io
        println(io, "# Processing set-up", isempty(site) ? "" : " of site $site", " (Timekeepers.jl, TKProc Save set-up).")
        println(io, "# TKProc reads it when you select the site. Edit the values or delete the file.")
        for (k, note) in PROCSET_KEYS
            haskey(settings, k) || continue
            println(io, rpad(string(k), 20), rpad(string(settings[k]), 38), " # ", note)
        end
    end
    return path
end

"""
    parse_span(s) -> Union{Tuple{DateTime, DateTime}, Nothing}

A time span from the two date-times in `s`, as `"2026-10-20 14:00 to
2026-10-21 07:20"` (also with `T`, with seconds, or with `-`, `–` or `/`
between them). `nothing` for `"all"`, an empty text, or a text without two
date-times in order.
"""
function parse_span(s::AbstractString)
    lowercase(strip(s)) in ("all", "") && return nothing
    ms = collect(eachmatch(r"(\d{4})-(\d{2})-(\d{2})[ T](\d{1,2}):(\d{2})(?::(\d{2}))?", s))
    length(ms) == 2 || return nothing
    dt(m) = DateTime(parse.(Int, m.captures[1:5])..., m.captures[6] === nothing ? 0 : parse(Int, m.captures[6]))
    a, b = try
        dt(ms[1]), dt(ms[2])
    catch
        return nothing
    end
    return a < b ? (a, b) : nothing
end

_fmt_span_or_nothing(s::AbstractString) = (sp = parse_span(s); sp === nothing ? nothing : _fmt_span(sp))

_fmt_time(t) = Dates.format(t, second(t) == 0 ? "yyyy-mm-dd HH:MM" : "yyyy-mm-dd HH:MM:SS")
_fmt_span(sp) = sp === nothing ? "all" : _fmt_time(sp[1]) * " to " * _fmt_time(sp[2])

# The set-up on the screen as the text values of procset.dat
function _procset_values(p::TKProc)
    box(k) = String(strip(p.boxes[k].stored_string[]))
    v = Dict{Symbol, String}(k => box(k) for k in keys(p.boxes))
    merge!(v, p.file_only)
    sp = get(v, :span, "all")
    v[:span] = lowercase(strip(sp)) in ("all", "") ? "all" : something(_fmt_span_or_nothing(sp), String(strip(sp)))
    base, remote, rate = p.base_menu.selection[], p.remote_menu.selection[], p.rate_menu.selection[]
    v[:base] = base === nothing ? "-" : string(base)
    v[:remote] = remote === nothing ? "-" : remote === :all ? "all" : remote isa AbstractVector ? join(remote, "+") : string(remote)
    v[:rate] = (rate === nothing || rate === :all) ? "all" : string(Float64(rate))
    v[:method] = string(something(p.method_menu.selection[], :ct2004))
    v[:leverage_weights] = p.leverage_toggle.active[] ? "on" : "off"
    return v
end

# Set a text box as if it were typed in
function _set_box!(b::Textbox, v::AbstractString)
    b.displayed_string[] = v
    b.stored_string[] = v
    return b
end

# The option of a menu whose value matches the text `v` of procset.dat
function _procset_option(menu, key::Symbol, v::AbstractString)
    w = lowercase(strip(v))
    hit(o) = key in (:base, :remote) ? (w == "-" ? o === nothing : w == "all" ? o === :all : o == v) :
               key === :rate ? (w == "all" ? o === :all : o isa Real && (r = tryparse(Float64, w)) !== nothing && isapprox(o, r; rtol = 1e-6)) :
               o isa Symbol && string(o) == w
    return findfirst(o -> hit(o[2]), menu.options[])
end

# Fill the menus and the boxes from procset.dat of `site`, and the options
# that it does not give with the defaults of the window. Without the file, the
# base and remote sites stay as _proc_focus! selected them. The result: the
# file ("" if none), the values that are not in a menu or not valid, and the
# error of a file that cannot be read
function _load_procset!(p::TKProc, site::SurveySite)
    file = procset_file(site.path)
    saved, err = Dict{Symbol, String}(), nothing
    if isfile(file)
        try
            saved = read_procset(file)
        catch e
            err = sprint(showerror, e)
        end
    else
        file = ""
    end
    v = merge(p.defaults, saved)
    p.file_only = Dict(k => get(v, k, d) for (k, d) in PROC_FILE_ONLY)
    missing = String[]
    p.updating = true
    try
        for (key, menu) in ((:base, p.base_menu), (:remote, p.remote_menu), (:rate, p.rate_menu), (:method, p.method_menu))
            haskey(v, key) || continue
            k = _procset_option(menu, key, v[key])
            # a base or remote site that the lists leave out (too near, no
            # position) is added when it is a site of the survey with Hx, Hy
            if k === nothing && key in (:base, :remote)
                names = String[strip(x) for x in split(v[key], r"[+,]") if !isempty(strip(x))]
                js = [_proc_site_index(p.survey, n) for n in names]
                if !isempty(js) && all(j -> j !== nothing && j != p.focus && has_magnetic(p.survey.sites[j]), js) &&
                   (key === :remote || length(js) == 1)
                    value = length(js) == 1 ? p.survey.sites[js[1]].name : [p.survey.sites[j].name for j in js]
                    menu.options[] = vcat(menu.options[], Any[(join(names, " + ") * " · set-up", value)])
                    k = length(menu.options[])
                end
            end
            k === nothing ? push!(missing, "$key $(v[key])") : (menu.i_selected[] = k)
        end
        haskey(v, :leverage_weights) && (p.leverage_toggle.active[] = lowercase(v[:leverage_weights]) in ("on", "true", "1", "yes"))
        for (key, b) in p.boxes
            haskey(v, key) || continue
            ok = try
                b.validator[](v[key])
            catch
                false
            end
            ok ? _set_box!(b, v[key]) : push!(missing, "$key $(v[key])")
        end
    finally
        p.updating = false
    end
    return (file = file, missing = missing, error = err)
end

# Save set-up: write the set-up on the screen to procset.dat of the site
function _save_procset!(p::TKProc)
    (p.busy || p.focus == 0) && return p
    site = p.survey.sites[p.focus]
    file = procset_file(site.path)
    _read_file_only!(p, site)
    try
        write_procset(file, _procset_values(p); site = site.name)
    catch err
        return _proc_status!(p, "Could not write $file: $(sprint(showerror, err))"; error = true)
    end
    _show_site_tables!(p, (file = file, missing = String[], error = nothing))
    return _proc_status!(p, "Saved the set-up of $(site.name) to $file. TKProc loads it when you select $(site.name)."; ok = true)
end

# The line next to Save set-up: where the set-up came from, and the directions
# of the channels in d/orientations.dat (read only: flip a channel in TKApp)
function _show_site_tables!(p::TKProc, saved)
    if p.focus == 0
        p.orient_label.text[] = ""
        return p
    end
    site = p.survey.sites[p.focus]
    setup = isempty(saved.file) ? "no procset.dat (defaults)" : "set-up from procset.dat"
    f = orientation_file(site.path)
    orient = if f === nothing
        "no d/orientations.dat: nominal directions"
    else
        try
            _, rows = read_orientations(f)
            row = get(rows, site.name, nothing)
            row === nothing ? "not in d/orientations.dat: nominal directions" :
                join(("$(ORIENTATION_NAMES[c]) $(row[c])" for c in ORIENTATION_COLUMNS), " · ")
        catch err
            "$(basename(f)): $(sprint(showerror, err))"
        end
    end
    p.orient_label.text[] = setup * " · " * orient
    return p
end

# The keywords of estimate_tf from the option row
# (the rate comes from its menu). The values go through procset.dat's text:
# process_site makes the same estimate from the saved file
function _proc_options(p::TKProc)
    o = procset_options(_procset_values(p))
    return (; (k => getfield(o, k) for k in keys(o) if k !== :rate)...)
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
    _read_file_only!(p, site)
    _relabel_all_remotes!(p)
    base = p.base_menu.selection[]
    remote = p.remote_menu.selection[]
    if remote === :all
        remote = filter(!=(base), _far_remotes(p, site))
        isempty(remote) && return _proc_status!(p, "No site with Hx, Hy recorded with $(site.name) at " *
                                                   "$(_all_remotes_label(p)[5:end]) or farther: lower min_distance in procset.dat, or select a remote site."; error = true)
    end
    base !== nothing && (remote == base || (remote isa AbstractVector && base in remote)) &&
        return _proc_status!(p, "$(base) is the base site of $(site.name): it cannot be a remote site too."; error = true)
    rate = something(p.rate_menu.selection[], :all)
    p.busy = true
    p.process_button.label[] = "Processing…"
    # a new estimate starts from an empty screen: the old one of the site
    # does not stay in view if this one fails
    delete!(p.results, site.name)
    _draw_proc!(p)
    return @async _run_estimate!(p, site.name; base, remote, rate)
end

# Clear: forget every estimate of the window and empty the plots and the
# status line. The loaded site, its menus and the options
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
        setup = _procset_values(p)
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
        file = procset_file(p.survey[name].path)
        tf.metadata[:procset] = setup
        p.results[name] = tf
        _draw_proc!(p)
        done = isempty(tf.periods) ? @sprintf("Done » no periods between min_period and max_period, in %.1f s.", time() - t0) :
               @sprintf("Done » %d periods, %.3g–%.3g s, in %.1f s.", length(tf.periods), first(tf.periods),
                        last(tf.periods), time() - t0)
        try
            write_procset(file, setup; site = name)
            tf.metadata[:procset_file] = file
            _show_site_tables!(p, (file = file, missing = String[], error = nothing))
            _proc_status!(p, done * " The set-up is in $file."; ok = true)
        catch err
            _proc_status!(p, done * " Could not write $file: $(sprint(showerror, err))"; error = true)
        end
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
    site = p.survey.sites[p.focus]
    isempty(tf.periods) && return _proc_status!(p, "$name has no periods to export."; error = true)
    try
        _read_file_only!(p, site)
        setup = copy(get(tf.metadata, :procset, _procset_values(p)))
        setup[:rotate] = p.file_only[:rotate]
        tf.metadata[:procset] = setup
        written = export_site(p.survey, tf; full_tensor = p.full_toggle.active[], errors = p.bars_toggle.active[],
                              wrap_phase = p.wrap_toggle.active[], rotate = _proc_rotation(p))
        _proc_status!(p, "Wrote " * join(basename.(written), ", ") * " to $(dirname(first(written)))."; ok = true)
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
    export_tf(dir, tf::TransferFunction; full_tensor = false, errors = true, wrap_phase = false,
              rotate = 0) -> Vector{String}

Write the estimate of one site into `dir` as `<name>.edi` ([`write_edi`](@ref))
and `<name>.png` ([`plot_tf`](@ref), with or without Zxx and Zyy and the
error bars). `<name>` is [`tf_filename`](@ref): the site, its base and remote
sites and the values of the options. The INFO block of the EDI file holds
each sensor with its calibration file or its dipole length and where the
length came from. The errors are the estimated errors, with no floor. With
`rotate` (degrees east of north), the plot shows the tensor turned by that
angle; the EDI file stays in geographic axes. The function returns the paths.
"""
function export_tf(dir::AbstractString, tf::TransferFunction; full_tensor::Bool = false, errors::Bool = true,
                   wrap_phase::Bool = false, rotate::Real = 0.0)
    mkpath(dir)
    name = tf_filename(tf)
    edi = write_edi(joinpath(dir, name * ".edi"), tf)
    png = joinpath(dir, name * (rotate == 0 ? "" : "-rot" * replace(@sprintf("%g", rotate), "." => "p")) * ".png")
    save(png, plot_tf(rotate == 0 ? tf : rotate_tf(tf, rotate); full_tensor, errors, wrap_phase); px_per_unit = 2)
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
