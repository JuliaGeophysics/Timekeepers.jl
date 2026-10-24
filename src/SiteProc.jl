# SiteProc.jl - the processing of a site from its procset.dat, in TKProc and without it.
# Author: @pankajkmishra
#
# TKProc keeps the set-up of a site in procset.dat (Save set-up). The same
# file, read here, gives the same estimate without the window:
#
#   survey = scan_survey(dir)
#   tf = process_site(survey, "3")          # reads 3/procset.dat
#   export_site(survey, tf)                  # writes 3/proc/
#
# Both TKProc (Process, Export…) and process_site go through procset_options
# and procset_references. Thus, the buttons and the function make the same
# estimate from the same file. export_site writes into the proc/ folder of the
# site: the EDI, the transfer function and a check of the data as PNG figures
# for a paper, and procnotes.md (how it was made, and how to make it again in
# TKProc).

"""
    remote_sites(survey, site; min_distance = 100.0, min_overlap = 0.75) -> Vector{String}

The sites with Hx and Hy at `min_distance` metres or farther that recorded
with `site` for `min_overlap` (a fraction) or more of its recording, nearest
first: "all ≥ …" of the Remote menu of TKProc, and `remote all` in
procset.dat. An estimate with several remote sites uses only the time that
all of them recorded: a site that recorded a short part of the time would cut
the data of the others, so it is left out (select it by itself to use it). A
site without a position is never in the list.
"""
function remote_sites(s::Survey, site; min_distance::Real = 100.0, min_overlap::Real = 0.75)
    return [c.site for c in _concurrent_remotes(s, _target_site(s, site))
            if 1000 * c.distance_km >= min_distance - 1.0e-6 && c.overlap_fraction >= min_overlap]
end

_shared_rate(a::SurveySite, b::SurveySite) =
    any(fs -> _gives_rate(a, fs) && _gives_rate(b, fs), vcat(site_rates(a), site_rates(b)))
function _concurrent_remotes(s::Survey, t::SurveySite)
    refs = site_references(s, t; rate = :all, remote_km = 0.0).remote
    return sort([c for c in refs if _shared_rate(t, s[c.site])]; by = c -> c.distance_km)
end

"""
    procset_options(settings) -> NamedTuple

The keywords of [`estimate_tf`](@ref) from the text values of procset.dat
([`read_procset`](@ref)). A key that `settings` leaves out takes the default
of a new TKProc window. `base`, `remote`, `min_distance` and `rotate` are not
keywords of the estimate: refer to [`procset_references`](@ref).
"""
function procset_options(v::AbstractDict)
    s(k, d) = (x = get(v, k, ""); isempty(strip(x)) ? d : strip(x))
    num(T, k, d) = something(tryparse(T, s(k, string(d))), d)
    top = s(:top_nyquist, "auto")
    rate = s(:rate, "all")
    return (; rate = lowercase(rate) == "all" ? :all : something(tryparse(Float64, rate), :all),
            window = num(Int, :section_length, 256), overlap = num(Float64, :section_overlap, 0.5),
            bands_per_decade = num(Int, :bands_per_decade, 8), prewhiten = num(Int, :ar_terms, 3),
            huber = num(Float64, :huber_a, 1.5), jackknife_groups = num(Int, :jackknife_blocks, 50),
            nyquist_fraction = something(tryparse(Float64, top), :auto), min_harmonic = num(Int, :min_frequency_index, 4),
            min_windows = num(Int, :min_sections, 8), max_levels = num(Int, :decimation_levels, 12),
            min_period = num(Float64, :min_period, 0.0),
            max_period = lowercase(s(:max_period, "inf")) == "inf" ? Inf : num(Float64, :max_period, Inf),
            span = _procset_span(s(:span, "all")),
            method = Symbol(lowercase(s(:method, "ct2004"))),
            leverage = lowercase(s(:leverage_weights, "on")) in ("on", "true", "1", "yes"))
end

function _procset_span(text)
    lowercase(strip(text)) in ("all", "") && return nothing
    sp = parse_span(text)
    sp === nothing && error("span \"$(text)\": write two times in order, e.g. 2026-10-20 14:00 to 2026-10-21 07:20, or all")
    return sp
end

"""
    procset_references(survey, site, settings) -> (base, remote)

The base site and the remote sites of procset.dat: `-` is none (for the base,
the site itself), `all` is [`remote_sites`](@ref) at `min_distance` metres
(100 if the file has none), and other values are site names (remote sites
joined with `+` or `,`).
"""
function procset_references(s::Survey, site, v::AbstractDict)
    t = _target_site(s, site)
    b = strip(get(v, :base, "-"))
    base = b in ("-", "", t.name) ? nothing : String(b)
    r = strip(get(v, :remote, "-"))
    dmin = something(tryparse(Float64, strip(get(v, :min_distance, "100"))), 100.0)
    remote = r in ("-", "") ? nothing :
             lowercase(r) == "all" ? filter(!=(base), remote_sites(s, t; min_distance = dmin)) :
             String[strip(x) for x in split(r, r"[+,]") if !isempty(strip(x))]
    remote isa Vector && base !== nothing && base in remote &&
        error("$(base) is the base site of $(t.name): it cannot be a remote site too")
    remote isa Vector && isempty(remote) &&
        error("No site with Hx, Hy recorded with $(t.name) at $(dmin) m or farther: lower min_distance or name a remote site")
    return (base = base, remote = remote)
end

"""
    process_site(survey, site; procset = procset_file(site.path), progress = nothing) -> TransferFunction

The transfer function of `site` (a name, an index or a
[`SurveySite`](@ref)) with its set-up in procset.dat: the same estimate as
Process in TKProc with that set-up. `procset` can also be a `Dict` of the
values. The result keeps the values in `tf.metadata[:procset]` and the file in
`tf.metadata[:procset_file]`.

See also [`export_site`](@ref), [`read_procset`](@ref).
"""
function process_site(s::Survey, site; procset = nothing, progress = nothing)
    t = _target_site(s, site)
    file = procset === nothing ? procset_file(t.path) : procset isa AbstractString ? String(procset) : ""
    v = procset isa AbstractDict ? Dict{Symbol, String}(procset) :
        isfile(file) ? read_procset(file) : error("$(t.name) has no $(file): save its set-up in TKProc first")
    o = procset_options(v)
    refs = procset_references(s, t, v)
    tf = estimate_tf(s, t.name; base = refs.base, remote = refs.remote, progress,
                     (k => getfield(o, k) for k in keys(o))...)
    tf.metadata[:procset] = v
    tf.metadata[:procset_file] = file
    return tf
end


const PAPER_AXIS = (spinewidth = 2.0, xtickwidth = 2.0, ytickwidth = 2.0, xticksize = 7, yticksize = 7,
                    xgridvisible = false, ygridvisible = false, titlefont = :regular, xlabelfont = :regular,
                    ylabelfont = :regular, backgroundcolor = :white)
const PAPER_FONTSIZE = 20
const PAPER_PX = 4

function _paper_axes!(axes)
    for ax in axes, (k, v) in pairs(PAPER_AXIS)
        setproperty!(ax, k, v)
    end
    return axes
end

_panel_label!(pos, k::Integer) = Label(pos, "($('a' + k - 1))"; fontsize = PAPER_FONTSIZE + 2, font = :regular,
                                        padding = (0, 10, 6, 0), halign = :left, tellwidth = false)

"""
    paper_tf(tf; full_tensor = false, errors = true, wrap_phase = false, rotate = 0.0,
             size = (1100, 720)) -> Figure

The figure of [`plot_tf`](@ref) for a paper: thick axes without a grid, panel
labels (a)-(d), the set-up in the title. `rotate` turns `tf` first (degrees).
Save it as a PNG with `save(path, fig; px_per_unit = 4)`.
"""
function paper_tf(tf::TransferFunction; full_tensor::Bool = false, errors::Bool = true, wrap_phase::Bool = false,
                  rotate::Real = 0.0, size = (1100, 720))
    t = rotate == 0 ? tf : rotate_tf(tf, rotate)
    fig = Figure(; size = size, fontsize = PAPER_FONTSIZE, figure_padding = (12, 24, 8, 8), backgroundcolor = :white)
    title = "Site $(tf.site) · " * _setup_text(tf) * (rotate == 0 ? "" : @sprintf(" · axes turned %+.1f°", rotate))
    Label(fig[1, 1], title; font = :regular, fontsize = PAPER_FONTSIZE, halign = :left, tellwidth = false)
    grid = GridLayout(fig[2, 1])
    axes = _paper_axes!(_tf_axes!(grid))
    mk(c) = MarkerElement(; PROC_MARKER..., color = c)
    zs = full_tensor ? [(2, "Zxy"), (3, "Zyx"), (1, "Zxx"), (4, "Zyy")] : [(2, "Zxy"), (3, "Zyx")]
    leg(pos, items, labels) = Legend(pos, items, labels; orientation = :horizontal, framevisible = false,
                                     labelsize = PAPER_FONTSIZE - 3, patchsize = (14, 12), halign = :left,
                                     tellwidth = false, padding = (0, 0, 0, 0))
    leg(grid[1, 1, Top()], [mk(PROC_ZCOLOURS[i]) for (i, _) in zs], last.(zs))
    leg(grid[1, 2, Top()], [mk(PROC_RE), mk(PROC_IM)], ["Re T", "Im T"])
    axes[1].ylabel = "ρa (Ω·m)"
    axes[2].ylabel = "φ (°)"
    _draw_tf!(axes, t; full_tensor, bars = errors, wrap_phase)
    for (k, (i, j)) in enumerate(((1, 1), (2, 1), (1, 2), (2, 2)))
        _panel_label!(grid[i, j, TopLeft()], k)
    end
    colgap!(grid, 40)
    return fig
end

"""
    paper_qc(tf; size = (1100, 820)) -> Figure

The check of an estimate for a paper: (a) the relative error of Zxy and Zyx,
(b) the coherence of Ex, Ey (and Hz) with the inputs and of the inputs with
the reference, (c) the windows in each band, (d) the data that each rate
used.
"""
function paper_qc(tf::TransferFunction; size = (1100, 820))
    fig = Figure(; size = size, fontsize = PAPER_FONTSIZE, figure_padding = (12, 24, 8, 8), backgroundcolor = :white)
    Label(fig[1, 1:2], "Site $(tf.site) · " * _setup_text(tf) * " · data check"; font = :regular, halign = :left, tellwidth = false)
    T = tf.periods
    xl = isempty(T) ? (1.0, 1.0e3) : (first(T) / 1.5, last(T) * 1.5)
    ax(pos; kw...) = Axis(pos; xscale = log10, xlabel = "Period (s)", PAPER_AXIS..., kw...)
    nice = [0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1, 2, 5, 10, 20, 50, 100, 200, 500, 1000, 2000, 5000]
    a1 = ax(fig[2, 1]; yscale = log10, ylabel = "|δZ| / |Z|", yticks = (nice, [string(x) for x in nice]))
    a2 = ax(fig[2, 2]; ylabel = "Coherence", limits = (nothing, (0, 1.02)))
    a3 = ax(fig[3, 1]; yscale = log10, ylabel = "Windows", yticks = (nice, [string(round(Int, x)) for x in nice]))
    for (i, j, c, l) in ((1, 2, PROC_ZCOLOURS[2], "Zxy"), (2, 1, PROC_ZCOLOURS[3], "Zyx"))
        r = sqrt.(tf.Z_var[i, j, :]) ./ abs.(tf.Z[i, j, :])
        ok = findall(x -> isfinite(x) && x > 0, r)
        isempty(ok) || scatterlines!(a1, T[ok], r[ok]; color = c, linewidth = 2, PROC_MARKER..., label = l)
    end
    hlines!(a1, [0.1]; color = :grey60, linestyle = :dash, linewidth = 1.5)
    axislegend(a1; position = :lt, framevisible = false, labelsize = PAPER_FONTSIZE - 4)
    for (k, c, l) in ((1, PROC_ZCOLOURS[2], "Ex"), (2, PROC_ZCOLOURS[3], "Ey"), (3, PROC_ZCOLOURS[1], "Hz"))
        y = tf.coherence[k, :]
        ok = findall(isfinite, y)
        isempty(ok) || scatterlines!(a2, T[ok], y[ok]; color = c, linewidth = 2, PROC_MARKER..., label = l)
    end
    for (k, l, st) in ((1, "Hx·ref", :dash), (2, "Hy·ref", :dot))
        y = tf.ref_coherence[k, :]
        ok = findall(isfinite, y)
        isempty(ok) || lines!(a2, T[ok], y[ok]; color = :black, linewidth = 2, linestyle = st, label = l)
    end
    axislegend(a2; position = :lb, framevisible = false, labelsize = PAPER_FONTSIZE - 4, nbanks = 2)
    nw = Float64.(tf.n_windows)
    ok = findall(>(0), nw)
    isempty(ok) || scatterlines!(a3, T[ok], nw[ok]; color = :black, linewidth = 2, PROC_MARKER..., markercolor = :grey80)
    for a in (a1, a2, a3)
        a.xticks = _decade_ticks(xl...)
        xlims!(a, xl...)
    end
    a4 = Axis(fig[3, 2]; xlabel = "UTC", yticksvisible = false, yticklabelsvisible = false, PAPER_AXIS...)
    md = tf.metadata
    rates = sort!(collect(get(md, :rates, Dict{Float64, Any}())); by = first, rev = true)
    if !isempty(rates)
        t0 = minimum(i[:start] for (_, i) in rates)
        t1 = maximum(i[:stop] for (_, i) in rates)
        h(t) = Dates.value(t - t0) / 3.6e6
        for (k, (fs, info)) in enumerate(rates)
            barplot!(a4, [k], [h(info[:stop])]; fillto = h(info[:start]), direction = :x,
                     color = :grey55, strokecolor = :black, strokewidth = 1.5, gap = 0.5)
            text!(a4, h(info[:start]), k + 0.38; text = @sprintf("%s · %.1f h in %d spans", _fs_label(fs), info[:hours], info[:spans]),
                  fontsize = PAPER_FONTSIZE - 5, align = (:left, :bottom))
        end
        a4.xlabel = "Hours from " * Dates.format(t0, "yyyy-mm-dd HH:MM") * " UTC"
        ylims!(a4, 0.3, length(rates) + 1.2)
        xlims!(a4, -0.5, h(t1) + 0.5)
    end
    for (k, (i, j)) in enumerate(((2, 1), (2, 2), (3, 1), (3, 2)))
        _panel_label!(fig[i, j, TopLeft()], k)
    end
    colgap!(fig.layout, 40)
    return fig
end


"""
    export_site(survey, tf; dir = joinpath(site.path, "proc"), full_tensor = false, errors = true,
                wrap_phase = false, rotate = nothing) -> Vector{String}

Write the estimate `tf` of a site of `survey` into `dir`, the `proc` folder of
the site (made new: the files of an earlier export there are removed):

- `<site>.edi` -- the transfer function in geographic axes;
- `<site>_tf.png` -- [`paper_tf`](@ref), turned by `rotate` degrees (the
  `rotate` of its procset.dat if `nothing`);
- `<site>_qc.png` -- [`paper_qc`](@ref);
- `procnotes.md` -- how the estimate was made ([`write_tf_report`](@ref)),
  with the set-up of procset.dat and the steps that make it again in TKProc.

The function returns the paths.
"""
function export_site(s::Survey, tf::TransferFunction; dir = nothing, full_tensor::Bool = false, errors::Bool = true,
                     wrap_phase::Bool = false, rotate = nothing)
    t = _target_site(s, tf.site)
    v = get(tf.metadata, :procset, Dict{Symbol, String}())
    θ = rotate === nothing ? something(tryparse(Float64, get(v, :rotate, "0")), 0.0) : Float64(rotate)
    out = dir === nothing ? joinpath(procset_file(t.path) |> dirname, "proc") : String(dir)
    mkpath(out)
    if dir === nothing
        for f in readdir(out; join = true)
            isfile(f) && (endswith(f, ".edi") || endswith(f, "_tf.png") || endswith(f, "_qc.png") ||
                          basename(f) == "procnotes.md") && rm(f)
        end
    end
    name = tf.site
    edi = write_edi(joinpath(out, name * ".edi"), tf)
    png = joinpath(out, name * "_tf.png")
    save(png, paper_tf(tf; full_tensor, errors, wrap_phase, rotate = θ); px_per_unit = PAPER_PX)
    qc = joinpath(out, name * "_qc.png")
    save(qc, paper_qc(tf); px_per_unit = PAPER_PX)
    notes = joinpath(out, "procnotes.md")
    files = [edi, png, qc]
    write_tf_report(notes, tf; files)
    _append_procnotes!(notes, s, t, tf, v, θ)
    return vcat(files, notes)
end

function _append_procnotes!(path, s::Survey, t::SurveySite, tf::TransferFunction, v, θ)
    open(path, "a") do io
        println(io, "\n## Set-up (procset.dat)\n")
        if isempty(v)
            println(io, "No procset.dat: the estimate came from the options above.\n")
        else
            println(io, "```")
            for (k, note) in PROCSET_KEYS
                haskey(v, k) && println(io, rpad(string(k), 20), rpad(v[k], 38), " # ", note)
            end
            println(io, "```\n")
            r = lowercase(strip(get(v, :remote, "-")))
            r == "all" && println(io, "`remote all` with `min_distance ", get(v, :min_distance, "100"),
                                  "` gave the remote sites ", join(tf.remote, ", "), " (nearest first).\n")
        end
        file = procset_file(t.path)
        saved = isfile(file) ? (try read_procset(file) catch; Dict{Symbol, String}() end) : Dict{Symbol, String}()
        diff = [k for (k, _) in PROCSET_KEYS if k !== :rotate && (haskey(v, k) || haskey(saved, k)) &&
                strip(get(saved, k, "")) != strip(get(v, k, ""))]
        val(k) = haskey(v, k) ? strip(v[k]) : "(the default)"
        label(k) = k === :remote && val(k) == "-" ? "Remote **none**" :
                   k === :base && val(k) == "-" ? "Base **$(t.name)**" : "$(k) **$(val(k))**"
        println(io, "## Make it again\n")
        println(io, "In TKProc:\n")
        println(io, "1. `run_tkproc(\"", replace(s.root, '\\' => '/'), "\")`, or Survey… and select that directory.")
        println(io, "2. Select site **", t.name, "**. TKProc reads its set-up from `", basename(dirname(file)), "/procset.dat`.")
        isempty(diff) || println(io, "   This is a variant of that set-up: set ", join(label.(diff), ", "), " in procset.dat before Process.")
        println(io, "3. Press **Process**", isempty(diff) ? ", then **Export…**: it writes this folder again." : ".", "\n")
        println(io, "Without the window:\n")
        if isempty(diff)
            println(io, "```julia\nusing Timekeepers\nsurvey = scan_survey(\"", replace(s.root, '\\' => '/'), "\")\n",
                    "tf = process_site(survey, \"", t.name, "\")\nexport_site(survey, tf)\n```")
        else
            changes = join((haskey(v, k) ? "v[:$(k)] = \"$(v[k])\"" : "delete!(v, :$(k))" for k in diff), "; ")
            println(io, "```julia\nusing Timekeepers\nsurvey = scan_survey(\"", replace(s.root, '\\' => '/'), "\")\n",
                    "v = read_procset(procset_file(survey[\"", t.name, "\"].path)); ", changes, "\n",
                    "tf = process_site(survey, \"", t.name, "\"; procset = v)\n",
                    "export_site(survey, tf; dir = \"", replace(dirname(path), '\\' => '/'), "\")\n```")
        end
        θ == 0 || @printf(io, "\nThe figure is in axes turned %+.1f° (`rotate`); the EDI is in geographic axes.\n", θ)
    end
    return path
end
