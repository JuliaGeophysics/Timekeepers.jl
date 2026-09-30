# Explorer.jl - the interactive GLMakie application.
# Author: @pankajkmishra
#
# Defines TKApp and the whole UI: stacked per-component time series with
# optional PSD panels, a switch per channel, a scrolling time window,
# drag-to-select with mask/unmask, a hover cursor that reads frequencies off
# the spectra panels, and loading or writing single runs and whole sites - a
# Metronix site one sampling rate at a time, behind a rate menu - with a
# loading window that follows the longer jobs.
#
# Three things keep it responsive on long records. Plotted series are
# decimated to a fixed bucket count by min/max per bucket, drawn into buffers
# the plot Observables already own so panning allocates almost nothing.
# Spectral recomputes are debounced behind a timer and reuse cached
# SpectralWorkspaces keyed by their configuration. A change of view builds
# only the panels it adds, and the spectrum panels after the first frame.

const TK_BLACK = RGBf(0.05, 0.05, 0.07)
const TK_BLUE = RGBf(0.114, 0.306, 0.847)
const TK_GREY = RGBf(0.42, 0.45, 0.50)
const TK_MUTED = RGBf(0.62, 0.66, 0.71)
const TK_FRAME = RGBf(0.55, 0.58, 0.63)
const TK_PANEL_BG = RGBf(0.992, 0.993, 0.995)
const TK_MASK_FILL = RGBAf(0.40, 0.43, 0.50, 0.18)
const TK_SEL_FILL = RGBAf(0.114, 0.306, 0.847, 0.18)
const TK_SEL_EDGE = RGBAf(0.114, 0.306, 0.847, 0.70)

# Spectral cursor. Teal for the pin, so it stands apart from the grey traces and
# the grey masked-interval fill.
const TK_CURSOR_LINE = RGBAf(0.30, 0.32, 0.34, 0.70)
const TK_PIN_LINE = RGBAf(0.18, 0.69, 0.78, 0.95)
const TK_PIN_GUIDE = RGBAf(0.18, 0.69, 0.78, 0.40)

# Half-width, in screen pixels, of the neighbourhood the cursor searches for a
# local maximum. A radius in pixels rather than in bins is the only rule that
# behaves the same across a log frequency axis: near DC two bins are tens of
# pixels apart, near Nyquist hundreds of bins fall inside one pixel.
const TK_SNAP_PIXELS = 6.0

# How far a candidate peak must stand above the plain nearest bin before the
# cursor jumps to it: a power ratio of about 1 dB. Without this the readout
# jitters bin to bin while the mouse crosses a flat noise floor.
const TK_SNAP_GAIN = 1.26

# Upper bound on the harmonic guides drawn above a pinned fundamental.
const TK_MAX_HARMONICS = 64

# Widest power range a PSD panel will show. Instrument anti-alias filters put a
# stopband ten or more decades below the passband right under Nyquist; scaling to
# include it squashes everything worth reading into the top of the panel.
const TK_PSD_MAX_DECADES = 8.0

const TK_LOGO_SKY   = RGBf(0.76, 0.84, 0.87)
const TK_LOGO_TEAL  = RGBf(0.18, 0.69, 0.78)
const TK_LOGO_SLATE = RGBf(0.30, 0.32, 0.34)

const TK_CTRL_ACCENT = RGBf(0.13, 0.55, 0.60)
const TK_CTRL_TRACK  = RGBAf(0.62, 0.66, 0.72, 0.35)

function _logo_menu(parent; kwargs...)
    return Menu(parent;
        cell_color_inactive_even = RGBAf(TK_LOGO_SKY.r, TK_LOGO_SKY.g, TK_LOGO_SKY.b, 0.35),
        cell_color_inactive_odd = RGBAf(TK_LOGO_SKY.r, TK_LOGO_SKY.g, TK_LOGO_SKY.b, 0.55),
        cell_color_active = RGBAf(TK_LOGO_TEAL.r, TK_LOGO_TEAL.g, TK_LOGO_TEAL.b, 0.85),
        cell_color_hover = RGBAf(TK_LOGO_TEAL.r, TK_LOGO_TEAL.g, TK_LOGO_TEAL.b, 0.30),
        selection_cell_color_inactive = RGBAf(TK_LOGO_SKY.r, TK_LOGO_SKY.g, TK_LOGO_SKY.b, 0.45),
        dropdown_arrow_color = TK_LOGO_SLATE,
        textcolor = TK_BLACK,
        kwargs...,
    )
end

"""
    _decade_ticks(max_ticks) -> (vmin, vmax) -> (values, labels)

Tick function for a log axis that only ticks whole decades, labelled 10ᵏ.
`LogTicks(WilkinsonTicks(n))` falls back to fractional exponents when a panel
spans few decades, which reads as 10^2.5 or 10^0.0. Past `max_ticks` decades it
keeps every second, third... one, on multiples so the labels stay regular. A
panel inside a single decade gets 1-2-5 steps instead, and one narrower than
those gets plain numbers.
"""
function _decade_ticks(max_ticks::Integer)
    return function (vmin, vmax)
        (vmin > 0 && vmax > vmin) || return (Float64[], Any[])
        lo, hi = log10(vmin), log10(vmax)
        ks = collect(ceil(Int, lo):floor(Int, hi))
        if length(ks) >= 2
            step = cld(length(ks), max_ticks)
            step > 1 && filter!(k -> mod(k, step) == 0, ks)
            return (exp10.(Float64.(ks)), [rich("10", superscript(string(k))) for k in ks])
        end
        steps = [(m, k) for k in floor(Int, lo):ceil(Int, hi) for m in (1, 2, 5)
                 if vmin <= m * exp10(k) <= vmax]
        if length(steps) < 2
            # narrower still: plain numbers, as on a linear axis
            vals = Makie.get_tickvalues(WilkinsonTicks(3), vmin, vmax)
            return (vals, Makie.get_ticklabels(Makie.automatic, vals))
        end
        return ([m * exp10(k) for (m, k) in steps],
                [k == 0 ? rich(string(m)) : rich("$(m)×10", superscript(string(k))) for (m, k) in steps])
    end
end

# The visible span is typed as a whole number of one of these units; "All"
# shows the whole record and ignores the number.
const WINDOW_UNITS = [
    ("seconds", 1.0),
    ("minutes", 60.0),
    ("hours",   3600.0),
    ("days",    86400.0),
    ("All",     Inf),
]

# A window length as typed: a positive whole number, nothing else.
_parse_window_count(s) = (n = tryparse(Int, strip(String(s))); n !== nothing && n > 0 ? n : nothing)

_window_span(count::Integer, unit_seconds::Real) = isfinite(unit_seconds) ? count * Float64(unit_seconds) : Inf

const VIEW_OPTIONS = [
    ("Time", :time),
    ("Spectra", :spectra),
    ("Time | Spectra", :time_spectra),
]

"""
    _shows_time(mode) / _shows_spectra(mode) -> Bool

Which panels a view mode draws. Code that depends on the view asks these rather
than comparing modes, so every mode is handled the same way.
"""
_shows_time(mode::Symbol) = mode === :time || mode === :time_spectra
_shows_spectra(mode::Symbol) = mode === :spectra || mode === :time_spectra

"""
    TKApp(; size = (1600, 900))
    TKApp(path::AbstractString; kwargs...)
    TKApp(ta::TimeArray; size = (1600, 900), source_format = :lemi424, source_path = "")

Build the interactive Timekeepers explorer over a time series, without opening
a window. Given a `path` the data is loaded first — a single file, or a site
directory whose runs are combined with gaps filled. Given nothing, the app
starts empty and waits for **Load**.

The app owns the data, a [`TimekeeperMask`](@ref) and the GLMakie figure. Read
the results back out with [`cleaned_timearray`](@ref), [`good_segments`](@ref),
[`sample_weights`](@ref), [`write_cleaned`](@ref) and [`write_mask`](@ref),
each of which accepts a `TKApp` directly.

Use [`run_tkapp`](@ref) to build and display in one step; `display(app)` opens
the window without blocking.

Requires a desktop session with OpenGL 3.3 or newer.
"""
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
    line_x_scratch::Vector{Float64}                    # redundant x output for channels 2..n
    window_seconds::Observable{Float64}
    window_start::Observable{Float64}
    slider::Any
    window_menu::Any                                   # unit of the visible span
    window_box::Any                                    # whole number of those units
    selection::Observable{Tuple{Float64, Float64}}
    selection_visible::Observable{Bool}
    mask_lows::Observable{Vector{Float64}}
    mask_highs::Observable{Vector{Float64}}
    source_format::Symbol
    source_path::String
    view_mode::Observable{Symbol}
    view_menu::Any
    psd_axes::Vector{Axis}
    psd_freqs::Vector{Observable{Vector{Float64}}}
    psd_values::Vector{Observable{Vector{Float64}}}
    psd_header::Observable{String}                     # the parameter line under the plots
    help_visible::Observable{Bool}                     # info badge toggled the glossary on
    spectral_workspaces::Dict{Tuple{Int, Float64, Int, Symbol, Symbol}, TKSpectralWorkspace}
    spectra_timer::Base.RefValue{Union{Nothing, Timer}}
    # Spectral cursor. The frequency Observables are empty when nothing is shown,
    # which is only safe because every plot drawn from them opts out of the axis
    # autolimits (VLines reports extrema of its own input as data limits, and
    # extrema throws on an empty vector). They live on the app rather than per
    # axis so one shared value lights up every channel panel at once.
    cursor_text::Vector{Observable{String}}            # inline readout, one per spectral panel
    cursor_panel::Base.RefValue{Int}                   # panel the cursor sits in, 0 = none
    cursor_freq::Observable{Vector{Float64}}           # hovered frequency [Hz], 0 or 1 element
    pin_freq::Observable{Vector{Float64}}              # pinned fundamental [Hz], 0 or 1 element
    pin_harmonics::Observable{Vector{Float64}}         # 2f0, 3f0, ... up to Nyquist
    pin_text::Observable{String}                       # label for the pinned fundamental
    # Channel switches, one per column. They survive a change of view and are
    # reset when a record is loaded; the checkboxes are rebuilt with the axes.
    channel_on::Vector{Bool}
    channel_boxes::Vector{Checkbox}
    psd_col::Int                                       # grid column of the spectra, 0 = none
    spectra_placeholder::Any                           # "Computing spectra…" while they are built
    view_generation::Int                               # bumped per view change; stale builds give up
    # Sampling rates of the loaded site, as the rate menu lists them. A Metronix
    # site can hold several; anything else has one. Only the rate on screen is in
    # memory. The others keep their masked intervals, a few timestamps each, so a
    # rate picked again gets its mask back and Write can still cut every rate.
    site_rates::Vector{Float64}
    rate_index::Int
    rate_intervals::Dict{Float64, Vector{Tuple{DateTime, DateTime}}}
    rate_menu::Any
    rate_menu_updating::Bool                           # true while the menu is refilled, not picked
end

function _component_color(name)
    s = lowercase(String(name))
    return startswith(s, "e") ? TK_BLUE : TK_BLACK
end

const _DISPLAY_LABELS = Dict(
    "bx" => "Bx", "by" => "By", "bz" => "Bz",
    "e1" => "Ex", "e2" => "Ey", "e3" => "E3", "e4" => "E4",
)

function _display_label(name)
    s = lowercase(String(name))
    return get(_DISPLAY_LABELS, s, String(name))
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
    names = collect(LEMI424_DEFAULT_COMPONENTS)
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
        ta, fmt = _load_lemi_xyz(path), :lemi_xyz
    else
        fmt = _detect_format(path)
        ta = fmt === :geomag ? load_geomag(path) : load_lemi424(path)
    end
    return _fill_time_gaps(ta), fmt
end

"""
    _is_metronix_run_path(path) -> Bool

Whether `path` names one Metronix run: an `.xml` or `.ats` file inside a
directory of `.ats` files, or such a directory itself. A file dialog cannot
select a directory, so picking the run's XML - or any of its channels - is how
a user selects it.
"""
function _is_metronix_run_path(path::AbstractString)
    isdir(path) && return _is_metronix_dir(path)
    isfile(path) || return false
    lowercase(splitext(path)[2]) in (".ats", ".xml") || return false
    return _is_metronix_dir(dirname(abspath(path)))
end

"""
    _load_metronix_run(path) -> (TimeArray, Symbol)

One Metronix run - the `.ats` channels of one run number at one rate plus
their XML - as a single TimeArray. `path` is the run's `.xml`, one of its
`.ats` files, or a `meas_*` directory holding only that run. `:site_dir` points
at the directory above the `meas_*` directory so a later Write can reuse the
site writer, and `:metronix_runs` narrows that write to this one run.
"""
function _load_metronix_run(path::AbstractString)
    run = read_metronix(path)
    filled = _fill_time_gaps(to_timearray(run; axis = :datetime))
    md = _ta_meta(filled)
    if md isa AbstractDict
        fs = sampling_rate(run)
        md[:source_format] = :metronix
        md[:site_dir] = dirname(run.metadata[:meas_dir])
        md[:metronix_runs] = [run.metadata[:metronix_run_id]]
        md[:sample_rate] = fs
        md[:metronix_rate] = fs
    end
    @info "Loaded Metronix run" run = basename(run.metadata[:metronix_run_id]) rate = _format_fs(sampling_rate(run))
    return filled, :metronix
end

"""
    _load_run_any(path) -> (TimeArray, Symbol)

Load one run from a data file or a Metronix run. Everything the Load Run button
can hand over goes through here.
"""
function _load_run_any(path::AbstractString)
    _is_metronix_run_path(path) && return _load_metronix_run(path)
    return _load_data_file(path)
end

const _SITE_DATA_EXTS = (".txt", ".dat", ".lem", ".xyz")

function _list_data_files(dir::AbstractString)
    site_name = _site_name_from_dir(dir)
    skip = Set(lowercase(site_name * ext) for ext in _SITE_DATA_EXTS)
    files = String[]
    for name in readdir(dir; sort = true)
        full = joinpath(dir, name)
        isfile(full) || continue
        ext = lowercase(splitext(name)[2])
        ext in _SITE_DATA_EXTS || continue
        lowercase(name) in skip && continue
        push!(files, full)
    end
    return files
end

function _site_name_from_dir(dir::AbstractString)
    s = rstrip(String(dir), ['/', '\\'])
    name = basename(s)
    isempty(name) && (name = basename(dirname(s)))
    isempty(name) ? "site" : name
end

function _combine_aux_columns(tas, total::Int, t_start::DateTime, step_ms::Int)
    template = nothing
    for ta in tas
        md = _ta_meta(ta)
        aux = md isa AbstractDict ? get(md, :aux_columns, nothing) : nothing
        if aux isa AbstractDict && !isempty(aux)
            template = aux
            break
        end
    end
    template === nothing && return nothing

    combined = Dict{Symbol, AbstractVector}()
    for (k, vec) in template
        if eltype(vec) <: AbstractString
            default = k === :lat_hemisphere ? "N" : k === :lon_hemisphere ? "E" : ""
            combined[k] = fill(default, total)
        else
            combined[k] = fill(NaN, total)
        end
    end

    for ta in tas
        md = _ta_meta(ta)
        aux = md isa AbstractDict ? get(md, :aux_columns, nothing) : nothing
        aux isa AbstractDict || continue
        ta_times = _ensure_datetime(_ta_timestamps(ta))
        for (k, src_vec) in aux
            haskey(combined, k) || continue
            dst_vec = combined[k]
            @inbounds for i in eachindex(ta_times)
                idx = Int(Dates.value(ta_times[i] - t_start) ÷ step_ms) + 1
                1 <= idx <= total || continue
                i <= length(src_vec) || continue
                dst_vec[idx] = src_vec[i]
            end
        end
    end
    return combined
end

function _combine_site_timearrays(tas::Vector{<:TimeArray}, site::AbstractString)
    isempty(tas) && error("No TimeArrays to combine")
    base = tas[1]
    names = _symbolize.(_ta_colnames(base))
    n_cols = length(names)
    fs = _sample_rate_from_timearray(base)
    fs > 0 || error("Cannot combine: sample rate must be positive (got $fs)")

    for i in 2:length(tas)
        fsi = _sample_rate_from_timearray(tas[i])
        isapprox(fsi, fs; rtol = 1e-6) ||
            @warn "Mixed sample rates across files; resampling to base grid" file_index = i base_fs = fs file_fs = fsi
    end

    step_ms = max(round(Int, 1000 / fs), 1)
    t_start = first(_ensure_datetime(_ta_timestamps(base)))
    t_end = last(_ensure_datetime(_ta_timestamps(base)))
    for ta in tas
        ts = _ensure_datetime(_ta_timestamps(ta))
        t_start = min(t_start, first(ts))
        t_end = max(t_end, last(ts))
    end
    total = Int(Dates.value(t_end - t_start) ÷ step_ms) + 1
    new_times = [t_start + Millisecond(step_ms * (i - 1)) for i in 1:total]
    new_vals = fill(NaN, total, n_cols)

    overlaps = 0
    for ta in tas
        ta_names = _symbolize.(_ta_colnames(ta))
        ta_vals = _ta_values(ta)
        ta_times = _ensure_datetime(_ta_timestamps(ta))
        col_map = Pair{Int, Int}[]
        for (j, name) in enumerate(ta_names)
            target = findfirst(==(name), names)
            target === nothing && continue
            push!(col_map, j => target)
        end
        @inbounds for (i, t) in enumerate(ta_times)
            idx = Int(Dates.value(t - t_start) ÷ step_ms) + 1
            1 <= idx <= total || continue
            for (src, dst) in col_map
                v = ta_vals[i, src]
                isfinite(v) || continue
                isfinite(new_vals[idx, dst]) && (overlaps += 1)
                new_vals[idx, dst] = v
            end
        end
    end
    overlaps > 0 &&
        @info "Site overlap: $overlaps sample-channels overlapped (later file wins)"

    base_meta = _ta_meta(base)
    meta = base_meta isa AbstractDict ? Dict{Symbol, Any}(base_meta) : Dict{Symbol, Any}()
    meta[:site] = String(site)
    meta[:start_time] = first(new_times)
    meta[:end_time] = last(new_times)
    meta[:n_samples] = total
    meta[:sample_rate] = fs
    meta[:source_file] = "<combined site: $(length(tas)) files>"
    meta[:n_files] = length(tas)
    combined_aux = _combine_aux_columns(tas, total, t_start, step_ms)
    combined_aux === nothing ? delete!(meta, :aux_columns) : (meta[:aux_columns] = combined_aux)
    return TimeArray(new_times, new_vals, names, meta)
end

"""
    _concat_runs_timearrays(tas, site) -> TimeArray

Join runs of one rate end to end in time order, without a common grid. Above
1 kHz a millisecond grid cannot hold the samples, so the gap between two runs
stays a jump in the time axis; the plot breaks its line there.
"""
function _concat_runs_timearrays(tas::Vector{<:TimeArray}, site::AbstractString)
    isempty(tas) && error("No TimeArrays to combine")
    order = sortperm([first(_ensure_datetime(_ta_timestamps(ta))) for ta in tas])
    tas = tas[order]
    names = _symbolize.(_ta_colnames(tas[1]))
    for ta in tas
        _symbolize.(_ta_colnames(ta)) == names ||
            error("Cannot join runs with different channels: $(names) vs $(_symbolize.(_ta_colnames(ta)))")
    end
    times = reduce(vcat, (_ensure_datetime(_ta_timestamps(ta)) for ta in tas))
    vals = reduce(vcat, (_ta_values(ta) for ta in tas))
    base_meta = _ta_meta(tas[1])
    meta = base_meta isa AbstractDict ? Dict{Symbol, Any}(base_meta) : Dict{Symbol, Any}()
    meta[:site] = String(site)
    meta[:start_time] = first(times)
    meta[:end_time] = last(times)
    meta[:n_samples] = length(times)
    meta[:source_file] = "<combined site: $(length(tas)) files>"
    meta[:n_files] = length(tas)
    delete!(meta, :aux_columns)
    return TimeArray(times, vals, names, meta; unchecked = true)
end

const TERM_BG = RGBAf(0.95, 0.95, 0.96, 1.0)
const TERM_FG = RGBf(0.10, 0.10, 0.12)
const TERM_DIM = RGBf(0.55, 0.55, 0.60)

# The loading window: its bar, and the marker in front of each log entry.
const PROGRESS_TRACK = RGBAf(0.55, 0.58, 0.63, 0.22)
const PROGRESS_PANEL = RGBf(1.0, 1.0, 1.0)
const PROGRESS_EDGE  = RGBAf(0.55, 0.58, 0.63, 0.35)
const PROGRESS_OK    = RGBf(0.16, 0.56, 0.33)
const PROGRESS_WARN  = RGBf(0.80, 0.52, 0.05)
const PROGRESS_ERROR = RGBf(0.80, 0.20, 0.18)

"""
    ProgressConsole

The loading window behind Load Site and the Metronix write: a progress bar, a
line naming the current step with a counter and a clock, and a log.

Workers report from any thread through [`_progress_step!`](@ref),
[`_progress_note!`](@ref) and [`_progress_finish!`](@ref). Those only touch the
plain fields under `lock`; [`_flush_progress!`](@ref), on the task that owns
the window, copies them into the Observables the window draws.
"""
mutable struct ProgressConsole
    screen::Any
    lock::ReentrantLock
    dirty::Threads.Atomic{Bool}
    t0::Float64
    entries::Vector{Tuple{Float64, Symbol, String}}  # (seconds since open, kind, text)
    max_lines::Int
    activity::String
    step::Int
    total::Int                                       # 0 while the length of the job is unknown
    state::Symbol                                    # :busy, then :ok, :warn or :error
    finished_at::Float64
    log_rows::Vector{NTuple{4, Observable}}          # per row: time, marker colour, text, text colour
    activity_obs::Observable{String}
    counter_obs::Observable{String}
    bar_obs::Observable{Tuple{Float64, Float64}}     # (start, width) as fractions of the track
    state_obs::Observable{Symbol}
end

"""
    _progress_note!(console, kind, msg)

Add `msg` to the log. `kind` picks the marker colour: `:info` and `:step` for
work in progress, `:ok`, `:warn`, `:error`, or `:detail` for an unmarked line
that belongs to the entry above it. Every kind but `:detail` also becomes the
activity line under the bar. A `nothing` console ignores it.
"""
_progress_note!(::Nothing, _kind::Symbol, _msg::AbstractString) = nothing
function _progress_note!(p::ProgressConsole, kind::Symbol, msg::AbstractString)
    lock(p.lock) do
        t = time() - p.t0
        lines = split(String(msg), '\n')
        for (k, line) in enumerate(lines), (j, piece) in enumerate(_wrap_words(line, PROGRESS_WRAP_CHARS))
            # the continuation of a wrapped line reads as a detail: no marker, no time
            push!(p.entries, (j == 1 ? t : NaN, k == 1 && j == 1 ? kind : :detail, piece))
        end
        length(p.entries) > p.max_lines &&
            deleteat!(p.entries, 1:(length(p.entries) - p.max_lines))
        kind === :detail || (p.activity = String(first(lines)))
    end
    p.dirty[] = true
    yield()
    return nothing
end

_progress_println!(p, msg::AbstractString) = _progress_note!(p, :detail, msg)

"""
    _progress_step!(console, i, n, msg)

Log step `i` of `n` and move the bar: it shows `i - 1` of `n` done, scaled to
90% so the work after the last step - sorting, combining - still has room.
"""
_progress_step!(::Nothing, _i::Integer, _n::Integer, _msg::AbstractString) = nothing
function _progress_step!(p::ProgressConsole, i::Integer, n::Integer, msg::AbstractString)
    lock(p.lock) do
        p.step, p.total = i, n
    end
    _progress_note!(p, :step, msg)
end

"""
    _progress_finish!(console, kind, msg)

Close the job as `:ok`, `:warn` or `:error`: the bar fills in that colour, the
clock stops, and `msg` becomes the last log entry.
"""
_progress_finish!(::Nothing, _kind::Symbol, _msg::AbstractString) = nothing
function _progress_finish!(p::ProgressConsole, kind::Symbol, msg::AbstractString)
    lock(p.lock) do
        p.state = kind
        p.finished_at = time() - p.t0
    end
    _progress_note!(p, kind, msg)
    _flush_progress!(p)
    return nothing
end

function _progress_bar_span(state::Symbol, step::Int, total::Int, now::Float64)
    state === :busy || return (0.0, 1.0)
    total > 0 && return (0.0, clamp(0.9 * (step - 1) / total, 0.02, 1.0))
    # length unknown: a segment sweeping across the track
    w = 0.25
    x = mod(now * 0.7, 1 + w) - w
    lo, hi = clamp(x, 0.0, 1.0), clamp(x + w, 0.0, 1.0)
    return (lo, hi - lo)
end

_progress_marker_color(kind::Symbol) =
    kind === :ok ? PROGRESS_OK :
    kind === :warn ? PROGRESS_WARN :
    kind === :error ? PROGRESS_ERROR :
    kind === :step ? TK_BLUE : TERM_DIM

_progress_state_color(state::Symbol) = state === :busy ? TK_BLUE : _progress_marker_color(state)

# Makie wraps plain text only, and wrapping inside a row would break the grid,
# so long lines are split here, at spaces where there are any.
const PROGRESS_WRAP_CHARS = 78

function _wrap_words(line::AbstractString, width::Int)
    length(line) <= width && return [String(line)]
    out = String[]
    current = ""
    for word in split(line, ' ')
        while length(word) > width                     # a path with no spaces
            isempty(current) || (push!(out, current); current = "")
            push!(out, first(word, width))
            word = chop(word; head = width, tail = 0)
        end
        candidate = isempty(current) ? String(word) : current * " " * word
        if length(candidate) > width
            push!(out, current)
            current = String(word)
        else
            current = candidate
        end
    end
    isempty(current) || push!(out, current)
    return out
end

const PROGRESS_CLEAR = RGBAf(0, 0, 0, 0)

# Rows fill from the bottom, so the newest entry always sits on the last row.
function _progress_fill_rows!(rows, entries)
    offset = length(rows) - length(entries)
    for (r, (time_obs, marker_obs, text_obs, color_obs)) in enumerate(rows)
        k = r - offset
        t, kind, text = k >= 1 ? entries[k] : (NaN, :detail, "")
        time_str = isnan(t) ? "" : @sprintf("%.1f s", t)
        marker = kind === :detail ? PROGRESS_CLEAR : RGBAf(_progress_marker_color(kind))
        color = RGBAf(kind === :detail ? TERM_DIM : kind === :error ? PROGRESS_ERROR : TERM_FG)
        time_obs[] == time_str || (time_obs[] = time_str)
        marker_obs[] == marker || (marker_obs[] = marker)
        text_obs[] == text || (text_obs[] = text)
        color_obs[] == color || (color_obs[] = color)
    end
    return nothing
end

_flush_progress!(::Nothing) = nothing
function _flush_progress!(p::ProgressConsole)
    # the bar and the clock move on every pump tick, new entries or not
    state, step, total, activity, t_end, now, dirty = lock(p.lock) do
        (p.state, p.step, p.total, p.activity, p.finished_at, time() - p.t0, p.dirty[])
    end
    elapsed = state === :busy ? now : t_end
    counter = (total > 0 ? "$(min(step, total)) / $(total)    " : "") * @sprintf("%.1f s", elapsed)
    span = _progress_bar_span(state, step, total, now)
    p.bar_obs[] == span || (p.bar_obs[] = span)
    p.counter_obs[] == counter || (p.counter_obs[] = counter)
    p.state_obs[] === state || (p.state_obs[] = state)
    dirty || return nothing
    p.dirty[] = false
    entries = lock(() -> copy(p.entries), p.lock)
    _progress_fill_rows!(p.log_rows, entries)
    p.activity_obs[] = activity
    return nothing
end

function _run_with_progress_pump(f, console::Union{ProgressConsole, Nothing})
    worker = Threads.@spawn f()
    while !istaskdone(worker)
        _flush_progress!(console)
        sleep(0.05)
    end
    _flush_progress!(console)
    try
        return fetch(worker)
    catch err
        err isa TaskFailedException ? throw(err.task.exception) : rethrow()
    end
end

function _load_site_directory(dir::AbstractString;
                              progress::Union{ProgressConsole, Nothing} = nothing)
    isdir(dir) || error("Not a directory: $dir")
    files = _list_data_files(dir)
    isempty(files) && error("No data files ($(join(_SITE_DATA_EXTS, ", "))) found in: $dir")

    site_name = _site_name_from_dir(dir)
    _progress_note!(progress, :info, "Found $(length(files)) data file" *
                                      (length(files) == 1 ? "" : "s") * " in $(site_name)")
    _progress_println!(progress, dir)

    loaded = Tuple{TimeArray, Symbol}[]
    skipped = 0
    for (i, path) in enumerate(files)
        _progress_step!(progress, i, length(files), "Reading $(basename(path))")
        try
            ta, fmt = _load_data_file(path)
            push!(loaded, (ta, fmt))
        catch err
            @warn "Skipping file (could not load)" path exception = err
            _progress_note!(progress, :warn, "Skipped $(basename(path)): it could not be parsed")
            skipped += 1
        end
    end
    isempty(loaded) && error("Could not load any data files from: $dir")

    formats = unique(t[2] for t in loaded)
    length(formats) > 1 &&
        @warn "Mixed file formats in site; treating combined output as $(first(formats))" formats
    fmt = first(formats)

    _progress_note!(progress, :info, "Sorting $(length(loaded)) runs by start time")
    sort!(loaded; by = x -> first(_ensure_datetime(_ta_timestamps(x[1]))))

    _progress_note!(progress, :info, "Combining $(length(loaded)) runs onto one time grid")
    combined = _combine_site_timearrays([t[1] for t in loaded], site_name)
    meta = _ta_meta(combined)
    @info "Loaded site" site = site_name n_files = length(loaded) skipped span = "$(meta[:start_time]) → $(meta[:end_time])"

    _progress_note!(progress, skipped > 0 ? :warn : :ok,
                    "Combined $(length(loaded)) runs" * (skipped > 0 ? ", skipped $skipped" : ""))
    _progress_println!(progress, "$(meta[:n_samples]) samples, $(meta[:start_time]) to $(meta[:end_time])")
    return combined, fmt
end

"""
    load_metronix_site(site_dir; rate = nothing) -> TimeArray

Read the runs of one sampling rate from a Metronix site - or from a single
`meas_*` directory, which can hold runs at several rates - in time order, as
one record. Up to 1 kHz the runs share one time grid with the gaps between
them filled with `NaN`; above it they are joined end to end, the gaps left as
jumps in time. A site holding a single rate needs no `rate`; for a mixed-rate
site it is required, and [`metronix_site_rates`](@ref) lists the choices
without reading any samples.

A raw site of `meas_*` directories is first separated by rate with
[`split_metronix_site`](@ref) - `DF002` into `DF002.TK/128`, `DF002.TK/4096`,
... - and the runs are read from that copy.
"""
function load_metronix_site(site_dir::AbstractString; rate::Union{Nothing, Real} = nothing)
    rates = metronix_site_rates(site_dir)
    isempty(rates) && error("No Metronix meas_ directories found in: $site_dir")
    if rate === nothing
        length(rates) == 1 ||
            error("$site_dir holds $(length(rates)) sampling rates ($(join(_format_fs.(rates), ", "))); pass `rate`")
        rate = only(rates)
    end
    return first(_load_metronix_site(site_dir; rate = rate))
end

function _load_metronix_site(dir::AbstractString;
                             progress::Union{ProgressConsole, Nothing} = nothing,
                             rate::Union{Nothing, Real} = nothing)
    # A raw site is read from its copy separated by rate, <site>.TK, made
    # here unless an earlier load already made it.
    if _is_raw_metronix_site(dir)
        name = _site_name_from_dir(dir)
        split_name = name * _TK_SUFFIX
        _progress_note!(progress, :info, "Checking whether $(name) is split by sampling rate")
        if metronix_site_is_split(dir)
            _progress_note!(progress, :ok, "$(name) is already split by sampling rate in $(split_name)")
            dir = _norm_path(dir) * _TK_SUFFIX
        else
            _progress_note!(progress, :info, "Splitting $(name) by sampling rate into $(split_name)")
            dir = split_metronix_site(dir; on_run = (i, n, id) ->
                _progress_step!(progress, i, n, "Copying $(basename(id))"))
            _progress_note!(progress, :ok, "Split $(name) by sampling rate into $(split_name)")
        end
    elseif _metronix_layout_root(dir) != _norm_path(dir) || any(_is_rate_dirname, readdir(dir))
        _progress_note!(progress, :info, "$(_site_name_from_dir(dir)) is already split by sampling rate")
    end
    runs = metronix_site_runs(dir)
    isempty(runs) && error("No Metronix meas_ directories found in: $dir")
    rates = sort(collect(keys(runs)))
    # One TimeArray carries one sample rate, so a multi-rate site has to be
    # narrowed to one. The GUI asks first and passes the answer in; a silent
    # caller falls back to the lowest rate with a warning.
    chosen = if rate === nothing
        length(rates) > 1 && @warn "Metronix site has multiple sampling rates; loading the lowest. " *
            "Pass `rate` to pick another." rates
        first(rates)
    else
        key = round(Float64(rate); digits = 6)
        haskey(runs, key) ||
            error("No Metronix runs at $(_format_fs(key)) in: $dir (have $(join(_format_fs.(rates), ", ")))")
        key
    end
    run_ids = runs[chosen]                               # in start-time order
    site_name = _site_name_from_dir(dir)

    _progress_note!(progress, :info, "Found $(length(run_ids)) Metronix run" *
                                      (length(run_ids) == 1 ? "" : "s") *
                                      " at $(_format_fs(chosen)) in $(site_name)")
    _progress_println!(progress, String(dir))

    exact_fs = Float64(chosen)
    tas = TimeArray[]
    for (i, id) in enumerate(run_ids)
        _progress_step!(progress, i, length(run_ids), "Reading $(basename(id))")
        run = read_metronix(id)
        i == 1 && (exact_fs = sampling_rate(run))
        push!(tas, to_timearray(run; axis = :datetime))
    end
    ta = length(tas) == 1 ? tas[1] :
         exact_fs > 1000 ? _concat_runs_timearrays(tas, site_name) :
         _combine_site_timearrays(tas, site_name)
    filled = _fill_time_gaps(ta)
    md = _ta_meta(filled)
    if md isa AbstractDict
        md[:source_format] = :metronix
        md[:site_dir] = _metronix_site_root(dir)
        md[:metronix_source_dir] = _norm_path(dir)     # where the rate menu reads from
        md[:metronix_runs] = run_ids
        md[:sample_rate] = exact_fs
        md[:metronix_rate] = exact_fs
    end
    _progress_note!(progress, :ok, "Read $(length(run_ids)) run" *
                                    (length(run_ids) == 1 ? "" : "s") * " at $(_format_fs(chosen))")
    return filled, :metronix
end

# The site a Write covers: the parent of a meas_ directory picked on its own,
# the <site>.TK directory above one of its rate directories, otherwise the
# directory itself.
_metronix_site_root(dir::AbstractString) =
    _is_metronix_dir(dir) ? dirname(_norm_path(dir)) : _metronix_layout_root(dir)

function _try_set_transparent_framebuffer(value::Bool)
    try
        GLMakie.GLFW.WindowHint(GLMakie.GLFW.TRANSPARENT_FRAMEBUFFER, value)
    catch
    end
    return nothing
end

function _center_and_float_window!(screen)
    screen === nothing && return nothing
    try
        glwin = screen.glscreen
        try
            GLMakie.GLFW.SetWindowAttrib(glwin, GLMakie.GLFW.FLOATING, true)
        catch
        end
        mon = GLMakie.GLFW.GetPrimaryMonitor()
        vmode = GLMakie.GLFW.GetVideoMode(mon)
        w, h = GLMakie.GLFW.GetWindowSize(glwin)
        x = (Int(vmode.width) - Int(w)) ÷ 2
        y = (Int(vmode.height) - Int(h)) ÷ 2
        GLMakie.GLFW.SetWindowPos(glwin, max(x, 0), max(y, 0))
    catch
    end
    return nothing
end

"""
    _ask_choice(prompt, labels) -> Union{Nothing, Int}

Open a small floating window with `prompt` above one button per label, and wait
for a pick; returns its index, or `nothing` when the window is closed. Like the
loading window it carries no title of its own.

It waits by polling, the same way [`_run_with_progress_pump`](@ref) does, so it
belongs on a background task - the load handlers already run inside `@async`.
On the render task the poll loop would starve the very window it waits on.
"""
function _ask_choice(prompt::AbstractString, labels::Vector{String})
    isempty(labels) && return nothing
    _try_set_transparent_framebuffer(true)
    fig = Figure(;
        size = (460, 96 + 44 * length(labels)),
        backgroundcolor = TERM_BG,
        figure_padding = (20, 20, 18, 20),
    )
    Label(fig[1, 1], String(prompt);
          color = TERM_FG, halign = :left, justification = :left, tellwidth = false)
    choices = GridLayout(fig[2, 1])
    # `width = nothing` stretches each choice across the window
    buttons = [Button(choices[i, 1]; label = labels[i], width = nothing)
               for i in eachindex(labels)]
    length(labels) > 1 && rowgap!(choices, 8)
    rowgap!(fig.layout, 14)

    # The click lands on the render task and the wait runs on ours, so the
    # handoff crosses threads and needs an atomic.
    picked = Threads.Atomic{Int}(0)
    for (i, b) in enumerate(buttons)
        on(b.clicks) do _
            picked[] = i
        end
    end

    screen = nothing
    try
        screen = display(GLMakie.Screen(; title = "Timekeepers",
                                          visible = true,
                                          focus_on_show = true),
                         fig)
        _center_and_float_window!(screen)
    catch err
        @warn "Could not open the chooser window" exception = err
        _try_set_transparent_framebuffer(false)
        return nothing
    end
    _try_set_transparent_framebuffer(false)

    while picked[] == 0 && isopen(screen)
        sleep(0.05)
    end
    idx = picked[]
    try
        close(screen)
    catch
    end
    return idx == 0 ? nothing : idx
end

"""
    _prompt_metronix_rate(dir) -> (Bool, Union{Nothing, Float64})

Ask which sampling rate to import from a Metronix site holding more than one;
only one is read at a time.

`(true, nothing)` when there is nothing to ask - a plain directory, a single
`meas_` run, or a site holding one rate. `(true, rate)` for a choice.
`(false, nothing)` when the window is closed, which cancels the import.
"""
function _prompt_metronix_rate(dir::AbstractString)
    is_metronix_site(dir) || return (true, nothing)
    runs = metronix_site_runs(dir)
    rates = sort(collect(keys(runs)))
    length(rates) > 1 || return (true, nothing)
    labels = ["$(_format_fs(r))  ·  $(length(runs[r])) run" * (length(runs[r]) == 1 ? "" : "s")
              for r in rates]
    _, empty_xmls = _metronix_site_index(dir)
    skipped = isempty(empty_xmls) ? "" :
              "\n$(length(empty_xmls)) XML" * (length(empty_xmls) == 1 ? " describes" : "s describe") *
              " no recorded data and " * (length(empty_xmls) == 1 ? "is" : "are") * " skipped."
    idx = _ask_choice("$(_site_name_from_dir(dir)) holds $(length(rates)) sampling rates.\n" *
                      "Pick the one to load; the rate menu switches later." * skipped, labels)
    idx === nothing && return (false, nothing)
    return (true, rates[idx])
end

"""
    _show_progress_window(; window_size, max_lines) -> ProgressConsole

Open the loading window: a progress bar, the current step with a counter and a
clock beneath it, and a log panel that keeps the last `max_lines` entries. It
carries no title or heading of its own; the steps say what is happening.
"""
function _show_progress_window(; window_size = (720, 440), max_lines::Int = 13)
    _try_set_transparent_framebuffer(true)
    fig = Figure(;
        size = window_size,
        backgroundcolor = TERM_BG,
        figure_padding = (20, 20, 20, 18),
    )
    activity_obs = Observable("Starting…")
    counter_obs = Observable("")
    bar_obs = Observable((0.0, 0.0))
    state_obs = Observable(:busy)

    # Bar: a track and a fill in the same cell. The fill's halign, as a number,
    # puts its left edge at `start` of the track.
    Box(fig[1, 1]; color = PROGRESS_TRACK, strokevisible = false, cornerradius = 3, height = 6)
    Box(fig[1, 1]; strokevisible = false, cornerradius = 3, height = 6,
        color = lift(_progress_state_color, state_obs),
        width = lift(((x, w),) -> Relative(max(w, 1e-3)), bar_obs),
        halign = lift(((x, w),) -> w >= 1 ? 0.0 : x / (1 - w), bar_obs))

    status = GridLayout(fig[2, 1])
    Label(status[1, 1], activity_obs;
          color = lift(st -> st === :error ? PROGRESS_ERROR : TERM_FG, state_obs),
          halign = :left, tellwidth = false)
    Label(status[1, 2], counter_obs; color = TERM_DIM, halign = :right)

    Box(fig[3, 1]; color = PROGRESS_PANEL, strokecolor = PROGRESS_EDGE,
        strokewidth = 1, cornerradius = 8)
    # One row per entry - elapsed time, marker, text - so the columns line up
    log = GridLayout(fig[3, 1]; alignmode = Outside(16, 16, 12, 12),
                     valign = :bottom, tellheight = false)
    log_rows = NTuple{4, Observable}[]
    for r in 1:max_lines
        time_obs, marker_obs = Observable(""), Observable(PROGRESS_CLEAR)
        text_obs, color_obs = Observable(""), Observable(RGBAf(TERM_FG))
        Label(log[r, 1], time_obs; color = TERM_DIM, halign = :right)
        Label(log[r, 2], "•"; color = marker_obs)
        Label(log[r, 3], text_obs; color = color_obs, halign = :left, tellwidth = false)
        push!(log_rows, (time_obs, marker_obs, text_obs, color_obs))
    end
    colsize!(log, 1, Fixed(52))
    colsize!(log, 3, Auto(true, 1.0))
    colgap!(log, 10)
    rowgap!(log, 4)

    rowsize!(fig.layout, 1, Fixed(6))
    rowsize!(fig.layout, 3, Auto(true, 1.0))
    rowgap!(fig.layout, 1, 14)
    rowgap!(fig.layout, 2, 12)

    screen = nothing
    try
        screen = display(GLMakie.Screen(; title = "Timekeepers",
                                          visible = true,
                                          focus_on_show = true),
                         fig)
        _center_and_float_window!(screen)
    catch err
        @warn "Could not open progress window; falling back to log only" exception = err
    end
    _try_set_transparent_framebuffer(false)
    return ProgressConsole(screen, ReentrantLock(), Threads.Atomic{Bool}(false), time(),
                           Tuple{Float64, Symbol, String}[], max_lines, "Starting…", 0, 0,
                           :busy, 0.0, log_rows, activity_obs, counter_obs, bar_obs, state_obs)
end

_close_progress_window!(::Nothing, _delay_s::Real = 0.0) = nothing
function _close_progress_window!(p::ProgressConsole, delay_s::Real = 0.0)
    p.screen === nothing && return
    delay_s > 0 && sleep(delay_s)
    try
        close(p.screen)
    catch
    end
    return
end

function _combined_site_path(dir::AbstractString, fmt::Symbol)
    site_name = _site_name_from_dir(dir)
    ext = _ext_for_format(fmt)
    return joinpath(rstrip(String(dir), ['/', '\\']), site_name * ext)
end

function _write_combined_site!(ta::TimeArray, dir::AbstractString,
                                fmt::Symbol; progress::Union{ProgressConsole, Nothing} = nothing)
    out_path = _combined_site_path(dir, fmt)
    _progress_note!(progress, :info, "Writing $(basename(out_path))")
    _progress_println!(progress, out_path)
    _write_data_file(out_path, ta, fmt)
    _progress_note!(progress, :ok, "Wrote $(basename(out_path))")
    return out_path
end

function _fill_time_gaps(ta::TimeArray)
    times = _ensure_datetime(_ta_timestamps(ta))
    n = length(times)
    n <= 1 && return ta
    fs = _sample_rate_from_timearray(ta)
    fs > 0 || return ta
    # Above 1 kHz several samples share a millisecond stamp, so a millisecond
    # grid would fold them together. Gaps stay as jumps in time instead.
    fs > 1000 && return ta
    step_ms = max(round(Int, 1000 / fs), 1)
    t0 = first(times)
    elapsed_ms = Dates.value(last(times) - t0)
    total = Int(elapsed_ms ÷ step_ms) + 1
    total <= n && return ta
    vals = _ta_values(ta)
    n_cols = size(vals, 2)
    new_vals = fill(NaN, total, n_cols)
    idx_map = Vector{Int}(undef, n)
    for i in 1:n
        idx = Int(Dates.value(times[i] - t0) ÷ step_ms) + 1
        idx_map[i] = idx
        1 <= idx <= total || continue
        @inbounds for j in 1:n_cols
            new_vals[idx, j] = vals[i, j]
        end
    end
    new_times = [t0 + Millisecond(step_ms * (i - 1)) for i in 1:total]
    names = _symbolize.(_ta_colnames(ta))
    meta = _ta_meta(ta)
    new_meta = meta isa AbstractDict ? Dict{Symbol, Any}(meta) : Dict{Symbol, Any}()
    new_meta[:n_samples] = total
    new_meta[:end_time] = last(new_times)
    aux = meta isa AbstractDict ? get(meta, :aux_columns, nothing) : nothing
    if aux isa AbstractDict && !isempty(aux)
        new_aux = Dict{Symbol, AbstractVector}()
        for (k, vec) in aux
            if eltype(vec) <: AbstractString
                default = k === :lat_hemisphere ? "N" : k === :lon_hemisphere ? "E" : ""
                padded = fill(default, total)
            else
                padded = fill(NaN, total)
            end
            @inbounds for i in 1:min(n, length(vec))
                idx = idx_map[i]
                1 <= idx <= total && (padded[idx] = vec[i])
            end
            new_aux[k] = padded
        end
        new_meta[:aux_columns] = new_aux
    end
    return TimeArray(new_times, new_vals, names, new_meta)
end

function _auto_mask_nan!(mask::TimekeeperMask, vals::AbstractMatrix)
    n_rows = size(vals, 1)
    n_rows == length(mask.masked) || return mask
    # A channel with no finite sample at all was not recorded - an electric-only
    # site, say - and would otherwise mask every row of the record.
    cols = [j for j in axes(vals, 2) if any(isfinite, view(vals, :, j))]
    @inbounds for i in 1:n_rows
        bad = false
        for j in cols
            if !isfinite(vals[i, j])
                bad = true
                break
            end
        end
        bad && (mask.masked[i] = true)
    end
    return _refresh_intervals!(mask)
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
    if n > 1
        dt_min_ms = minimum(Dates.value(times[i + 1] - times[i]) for i in 1:(n - 1))
        fs = dt_min_ms > 0 ? 1000.0 / dt_min_ms : 1.0
    else
        fs = 1.0
    end
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
    elseif source_format === :geomag
        return write_geomag(path, ta)
    end
    return write_lemi424(path, ta)
end

function _format_duration_compact(seconds::Real)
    s = max(0, round(Int, seconds))
    h, rem = divrem(s, 3600)
    m, sec = divrem(rem, 60)
    h > 0 && return "$(h)h$(m)m"
    m > 0 && return "$(m)m$(sec)s"
    return "$(sec)s"
end

function _format_fs(fs::Real)
    isfinite(fs) || return "—Hz"
    rounded = round(fs; digits = 4)
    return rounded == floor(rounded) ? "$(Int(rounded))Hz" : "$(rounded)Hz"
end

function _summary_text(ta::TimeArray)
    metadata = _ta_meta(ta)
    site = metadata isa AbstractDict ? get(metadata, :site, "—") : "—"
    n = length(_ta_timestamps(ta))
    fs = _sample_rate_from_timearray(ta)
    duration = (n > 0 && fs > 0) ? n / fs : 0.0
    return "$(site)  ·  $(_format_duration_compact(duration))  ·  $(_format_fs(fs))"
end

const _PLOT_BUCKETS = 2000

# Min and max of `col[a:c]`, split into clean and masked samples, as two points
# at the ends of the range.
function _push_minmax!(xs, ys_clean, ys_masked, secs, col, masked, a::Int, c::Int)
    clean_min = Inf
    clean_max = -Inf
    masked_min = Inf
    masked_max = -Inf
    @inbounds for i in a:c
        v = col[i]
        isfinite(v) || continue
        if masked[i]
            v < masked_min && (masked_min = v)
            v > masked_max && (masked_max = v)
        else
            v < clean_min && (clean_min = v)
            v > clean_max && (clean_max = v)
        end
    end
    push!(xs, secs[a])
    push!(ys_clean, isfinite(clean_min) ? Float32(clean_min) : NaN32)
    push!(ys_masked, isfinite(masked_min) ? Float32(masked_min) : NaN32)
    push!(xs, secs[c])
    push!(ys_clean, isfinite(clean_max) ? Float32(clean_max) : NaN32)
    push!(ys_masked, isfinite(masked_max) ? Float32(masked_max) : NaN32)
    return nothing
end

function _push_break!(xs, ys_clean, ys_masked, x)
    push!(xs, x)
    push!(ys_clean, NaN32)
    push!(ys_masked, NaN32)
    return nothing
end

"""
    _decimate_minmax!(xs, ys_clean, ys_masked, secs, col, masked, idx_lo, idx_hi, n_buckets; gap_s = Inf)

Reduce `col[idx_lo:idx_hi]` to a min/max pair per bucket for drawing. Where
consecutive samples lie more than `gap_s` seconds apart - two runs joined
without a grid - the line is broken rather than drawn across the gap. The
breaks depend on `secs` alone, so every channel gets the same `xs`.
"""
function _decimate_minmax!(
    xs::Vector{Float64},
    ys_clean::Vector{Float32},
    ys_masked::Vector{Float32},
    secs::Vector{Float64},
    col::AbstractVector{<:Real},
    masked::BitVector,
    idx_lo::Int,
    idx_hi::Int,
    n_buckets::Int;
    gap_s::Float64 = Inf,
)
    empty!(xs)
    empty!(ys_clean)
    empty!(ys_masked)
    n_window = idx_hi - idx_lo + 1
    n_window <= 0 && return 0

    if n_window <= 2 * n_buckets
        @inbounds for i in idx_lo:idx_hi
            i > idx_lo && secs[i] - secs[i - 1] > gap_s && _push_break!(xs, ys_clean, ys_masked, secs[i - 1])
            push!(xs, secs[i])
            v = Float32(col[i])
            if masked[i]
                push!(ys_clean, NaN32)
                push!(ys_masked, v)
            else
                push!(ys_clean, v)
                push!(ys_masked, NaN32)
            end
        end
        return length(xs)
    end

    @inbounds for b in 1:n_buckets
        a = idx_lo + ((b - 1) * n_window) ÷ n_buckets
        c = idx_lo + (b * n_window) ÷ n_buckets - 1
        c > idx_hi && (c = idx_hi)
        a > c && continue
        # a jump from the previous bucket into this one
        a > idx_lo && secs[a] - secs[a - 1] > gap_s && _push_break!(xs, ys_clean, ys_masked, secs[a - 1])
        seg = a
        if isfinite(gap_s)
            for i in a:(c - 1)
                if secs[i + 1] - secs[i] > gap_s
                    _push_minmax!(xs, ys_clean, ys_masked, secs, col, masked, seg, i)
                    _push_break!(xs, ys_clean, ys_masked, secs[i])
                    seg = i + 1
                end
            end
        end
        _push_minmax!(xs, ys_clean, ys_masked, secs, col, masked, seg, c)
    end
    return length(xs)
end

# Samples further apart than this are a gap between runs, not a sample step.
_plot_gap_seconds(app) = max(1.0, 2 / _sample_rate_from_timearray(app.data))

function _refresh_visible_lines!(app::TKApp)
    isempty(app.axes) && return app
    secs = app.time_seconds
    n = length(secs)
    n == 0 && return app
    n_channels = length(app.line_clean)
    n_channels == 0 && return app
    masked = app.mask.masked
    @assert length(masked) == n "Mask length $(length(masked)) != sample count $n"

    x_lo, x_hi = _visible_x_window(app)
    idx_lo = clamp(searchsortedfirst(secs, x_lo), 1, n)
    idx_hi = clamp(searchsortedlast(secs, x_hi), 1, n)
    if idx_lo > idx_hi
        empty!(app.line_x[])
        notify(app.line_x)
        for j in 1:n_channels
            empty!(app.line_clean[j][])
            empty!(app.line_masked[j][])
            notify(app.line_clean[j])
            notify(app.line_masked[j])
        end
        return app
    end

    gap_s = _plot_gap_seconds(app)
    col1 = @view app.raw_values[:, 1]
    _decimate_minmax!(app.line_x[], app.line_clean[1][], app.line_masked[1][],
        secs, col1, masked, idx_lo, idx_hi, _PLOT_BUCKETS; gap_s = gap_s)
    notify(app.line_x)                                 # buffers refilled in place
    notify(app.line_clean[1])
    notify(app.line_masked[1])

    for j in 2:n_channels
        col_j = @view app.raw_values[:, j]
        _decimate_minmax!(app.line_x_scratch, app.line_clean[j][], app.line_masked[j][],
            secs, col_j, masked, idx_lo, idx_hi, _PLOT_BUCKETS; gap_s = gap_s)
        notify(app.line_clean[j])                      # x is identical across channels
        notify(app.line_masked[j])
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
    _refresh_visible_lines!(app)
    _recompute_spectra!(app)
    _refresh_status!(app)
    return app
end

function _format_dt(app::TKApp, secs::Float64)
    return Dates.format(app.origin + Millisecond(round(Int, secs * 1000)), "yyyy-mm-dd HH:MM:SS")
end

function _refresh_status!(app::TKApp)
    n_masked = masked_samples(app.mask)
    n_intervals = length(app.mask.intervals)
    # Only the hints that apply: the drag verbs need a trace to drag on, which
    # the Spectra view does not draw.
    hints = String[]
    if app.selection_visible[]
        lo, hi = app.selection[]
        push!(hints, "Selection " * _format_dt(app, lo) * "  →  " * _format_dt(app, hi))
    elseif _shows_time(app.view_mode[])
        push!(hints, "Left-drag to select  ·  Right-click to mask the selection  ·  Right-drag = pan  ·  Scroll = zoom y")
    end
    if _shows_spectra(app.view_mode[])
        push!(hints, isempty(app.pin_freq[]) ?
            "Spectral panel: hover reads f  ·  Left-click pins f0 and its harmonics" :
            "Spectral panel: hover reads f  ·  Right-click clears the pin")
    end
    sel_text = join(hints, "  ·  ")
    app.status_label.text[] = "$(n_masked) masked samples in $(n_intervals) intervals    ·    $(sel_text)"
    return app
end

function _ready_status_text(app::TKApp)
    if isempty(app.source_path)
        return "Timekeepers ready - use Load to open a data file"
    end
    return "Timekeepers ready - $(basename(app.source_path))"
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
        lo, hi = t0 <= t1 ? (t0, t1) : (t1, t0)
        for (a, b) in app.mask.intervals
            if a <= hi && b >= lo
                lo = min(lo, a)
                hi = max(hi, b)
            end
        end
        unmask_interval!(app.mask, lo, hi)
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
        app.channel_on[j] || continue
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
    _refresh_visible_lines!(app)
    _autoscale_y!(app, x_lo, x_hi)
    return
end

function _visible_x_window(app::TKApp)
    ws = app.window_seconds[]
    span = app.span_seconds
    visible = isfinite(ws) ? min(ws, max(span, 1.0)) : max(span, 1.0)
    x_lo = app.window_start[]
    if span > visible
        x_lo = clamp(x_lo, 0.0, span - visible)
    else
        x_lo = 0.0
    end
    return (x_lo, x_lo + visible)
end

"""
    _visible_good_index_segments(app, x_lo, x_hi) -> (stretches, n_runs)

The unmasked stretches of the visible window, as index ranges, for the spectra
to estimate over. A stretch ends at a masked sample and also at a jump in time
between two runs joined end to end (above 1 kHz), so no FFT segment spans two
recordings - the gap is treated exactly like a mask. `n_runs` counts the runs
the window touches.
"""
function _visible_good_index_segments(app::TKApp, x_lo::Float64, x_hi::Float64)
    secs = app.time_seconds
    masked = app.mask.masked
    isempty(secs) && return (Tuple{Int, Int}[], 0)
    idx_lo = searchsortedfirst(secs, x_lo)
    idx_hi = searchsortedlast(secs, x_hi)
    idx_lo = clamp(idx_lo, 1, length(secs))
    idx_hi = clamp(idx_hi, 1, length(secs))
    idx_lo > idx_hi && return (Tuple{Int, Int}[], 0)
    gap_s = _plot_gap_seconds(app)
    segs = Tuple{Int, Int}[]
    n_runs = 1
    active = false
    start = idx_lo
    @inbounds for i in idx_lo:idx_hi
        if i > idx_lo && secs[i] - secs[i - 1] > gap_s
            n_runs += 1
            active && push!(segs, (start, i - 1))
            active = false
        end
        if !masked[i] && !active
            active = true
            start = i
        elseif masked[i] && active
            push!(segs, (start, i - 1))
            active = false
        end
    end
    active && push!(segs, (start, idx_hi))
    return segs, n_runs
end

# FFT segments Welch's method takes from stretches of these lengths: each
# stretch gives one per half-segment step, and a stretch shorter than nfft none.
_welch_segment_count(lengths, nfft::Integer) =
    sum((L >= nfft ? (L - nfft) ÷ (nfft ÷ 2) + 1 : 0 for L in lengths); init = 0)

function _current_nfft(app::TKApp)
    fs = _sample_rate_from_timearray(app.data)
    ws = app.window_seconds[]
    effective = isfinite(ws) ? ws : max(app.span_seconds, 1.0)
    return _auto_nfft(effective, fs), fs
end

"""
    _spectral_workspace!(app, nfft, fs; noverlap, window, detrend) -> TKSpectralWorkspace

Fetch the cached workspace for one spectral configuration, building it on first
use. Takes the app, the transform length, the sample rate and the window
parameters; returns a concretely typed workspace, so callers dispatch statically
rather than through an `Any` cache.
"""
function _spectral_workspace!(app::TKApp, nfft::Integer, fs::Real; noverlap::Integer = nfft ÷ 2,
    window::Symbol = :hann, detrend::Symbol = :mean)
    key = (Int(nfft), Float64(fs), Int(noverlap), window, detrend)
    return get!(app.spectral_workspaces, key) do
        SpectralWorkspace(nfft, fs; noverlap = noverlap, window = window, detrend = detrend)
    end
end

_fmt_hz(f::Real) = f >= 0.01 ? @sprintf("%.4f Hz", f) : @sprintf("%.2e Hz", f)

"""
    _fmt_hz_short(f) -> String

Frequency for the inline cursor readout, which sits over the trace where every
character costs space: four significant digits, trailing zeros dropped, giving
`0.1 Hz` where the panel header prints `0.1000 Hz`.
"""
_fmt_hz_short(f::Real) = @sprintf("%.4g Hz", f)

function _fmt_dur(t::Real)
    t < 60 && return @sprintf("%.1f s", t)
    t < 3600 && return @sprintf("%.1f min", t / 60)
    return @sprintf("%.1f h", t / 3600)
end

"""
    _fmt_delay(fs, f) -> (String, Bool)

Delay in samples a comb filter needs at rate `fs` to notch `f`, and whether that
delay lands on a whole sample. Only whole-sample delays are realisable without
interpolating the record, so a caller that shows this should flag the ones that
miss.

Nothing draws it today; the labels are tighter for leaving it out. It stays
because it is the one number a comb design needs from a pinned frequency, and
the derivation is easy to reach for when a caller wants it back.
"""
function _fmt_delay(fs::Real, f::Real)
    (f > 0 && fs > 0) || return ("-", false)
    d = fs / f
    dr = round(d)
    exact = dr >= 1 && abs(d - dr) <= 1.0e-3 * dr
    return (exact ? @sprintf("%d smp", Int(dr)) : @sprintf("%.2f smp", d), exact)
end

"""
    _psd_cursor_text(f, p) -> String

Inline readout for a PSD panel: frequency, its period, and the PSD there, kept
to one bracketed group so it reads as an annotation on the trace. The comb delay
that goes with the frequency belongs to the pin, where it holds still long
enough to write down.
"""
function _psd_cursor_text(f::Real, p::Real)
    f > 0 || return ""
    return @sprintf("[%s / %.4gs; PSD=%.3e]", _fmt_hz_short(f), 1 / f, p)
end

"""
    _pin_label_text(f) -> String

Label for the pinned fundamental, shown once in the top spectral panel. The
fundamental is a global quantity, so repeating it in every panel would be noise.
Same bracketed shape as [`_psd_cursor_text`](@ref), and nothing more: the pin
and its harmonic guides carry their own colour, which names the label already.
"""
function _pin_label_text(f::Real)
    f > 0 || return ""
    return @sprintf("[%s / %.4gs]", _fmt_hz_short(f), 1 / f)
end

"""
    _nearest_index(xs, x) -> Int

Index of the element of the sorted vector `xs` closest to `x`; 0 when `xs` is
empty.
"""
function _nearest_index(xs::AbstractVector{<:Real}, x::Real)
    n = length(xs)
    n == 0 && return 0
    hi = searchsortedfirst(xs, x)
    hi <= 1 && return 1
    hi > n && return n
    return (xs[hi] - x) < (x - xs[hi - 1]) ? hi : hi - 1
end

"""
    _peak_window(xs, ys, x, x_lo, x_hi) -> (near, best)

Two candidates for a cursor sitting at `x` on the sorted axis `xs`: `near` is the
plain nearest sample, `best` the largest finite `ys` between `x_lo` and `x_hi`.
The caller picks between them, because what counts as clearly a peak is a
property of the values being searched. Both are 0 when empty.
"""
function _peak_window(xs::AbstractVector{<:Real}, ys::AbstractVector{<:Real},
        x::Real, x_lo::Real, x_hi::Real)
    n = min(length(xs), length(ys))
    n == 0 && return (0, 0)
    near = _nearest_index(view(xs, 1:n), x)
    lo = clamp(searchsortedfirst(xs, x_lo), 1, n)
    hi = clamp(searchsortedlast(xs, x_hi), 1, n)
    best = near
    best_v = ys[near]
    @inbounds for i in lo:hi
        v = ys[i]
        if isfinite(v) && (!isfinite(best_v) || v > best_v)
            best_v = v
            best = i
        end
    end
    return (near, best)
end

"""
    _snap_halfwidth(ax, dim, px) -> Float64

Half-width of a snap window `px` screen pixels wide, in the axis' *transformed*
units along dimension `dim`. On a log axis that is a half-width in decades, so
the corresponding window in data units is multiplicative. Returns 0 before the
axis has a viewport, which degrades the cursor to plain nearest-bin.
"""
function _snap_halfwidth(ax::Axis, dim::Integer, px::Real)
    lims = ax.finallimits[]
    scale = dim == 1 ? ax.xscale[] : ax.yscale[]
    lo = scale(minimum(lims)[dim])
    hi = scale(maximum(lims)[dim])
    w = Makie.widths(ax.scene.viewport[])[dim]
    (w <= 0 || !isfinite(lo) || !isfinite(hi) || hi <= lo) && return 0.0
    return px * (hi - lo) / w
end

function _format_psd_header(nfft::Integer, fs::Real)
    df = fs / nfft
    return "nfft = $(nfft)   ·   df = $(_fmt_hz(df))   ·   f_Nyq = $(_fmt_hz(fs / 2))   ·   seg = $(_fmt_dur(nfft / fs))"
end

"""
    _spectra_help_text() -> String

What the parameter line means, shown in its place while the info badge is on.
The full account of how the spectra are computed is in the Spectral Views page
of the docs.
"""
_spectra_help_text() =
    "nfft: FFT length in samples  ·  df = fs/nfft: spacing between frequency bins" *
    "  ·  f_Nyq = fs/2: highest resolvable frequency  ·  seg = nfft/fs: length of one segment" *
    "  ·  how the spectra are computed: docs, Spectral Views"

_spectra_info_idle() =
    "Time view  ·  use the View menu to add Spectra panels for frequency content"

function _spectra_details_text(nfft::Integer, fs::Real; n_windows::Integer, n_stretches::Integer,
                               n_runs::Integer = 1)
    metrics = _format_psd_header(nfft, fs)
    avg = "averaged over $(n_windows) segment" * (n_windows == 1 ? "" : "s") *
          " from $(n_stretches) unmasked stretch" * (n_stretches == 1 ? "" : "es") *
          (n_runs > 1 ? " in $(n_runs) runs" : "")
    return "PSD · Welch's method     |     x: frequency [Hz], log  ·  y: PSD [amplitude^2/Hz], log" *
           "     |     Hann window, 50% overlap, mean-detrended; $(avg)     |     " * metrics
end

"""
    _autoscale_psd!(ax, psd)

Frame a PSD trace on its log y axis. The lower bound is capped at
`TK_PSD_MAX_DECADES` below the peak: an instrument anti-alias filter drops ten
or more decades in the last few bins before Nyquist, and scaling to include that
cliff squashes the whole usable band into the top sliver of the panel.
"""
function _autoscale_psd!(ax::Axis, psd::Vector{Float64})
    isempty(psd) && return
    ymin = Inf
    ymax = -Inf
    @inbounds for v in psd
        if isfinite(v) && v > 0
            v < ymin && (ymin = v)
            v > ymax && (ymax = v)
        end
    end
    (isfinite(ymin) && isfinite(ymax) && ymin < ymax) || return
    floor_v = ymax * exp10(-TK_PSD_MAX_DECADES)
    ylims!(ax, max(ymin, floor_v) * 0.5, ymax * 2.0)
    return
end

function _compute_psd_for_window!(app::TKApp)
    _shows_spectra(app.view_mode[]) || return app
    isempty(app.psd_axes) && return app
    x_lo, x_hi = _visible_x_window(app)
    segs, n_runs = _visible_good_index_segments(app, x_lo, x_hi)
    nfft, fs = _current_nfft(app)
    workspace = _spectral_workspace!(app, nfft, fs; noverlap = nfft ÷ 2)
    seg_lengths = Int[b - a + 1 for (a, b) in segs]
    too_short = !isempty(segs) && all(L -> L < nfft, seg_lengths)
    if too_short
        @warn "All visible good segments shorter than nfft" nfft maxlen = maximum(seg_lengths)
    end
    n_channels = length(app.psd_axes)
    n_used = 0
    for j in 1:n_channels
        if !app.channel_on[j]
            app.psd_freqs[j][] = Float64[]
            app.psd_values[j][] = Float64[]
            continue
        end
        seg_views = [view(app.raw_values, a:b, j) for (a, b) in segs]
        freqs, psd, nseg = _welch_psd_segments(seg_views, fs;
            nfft = nfft, noverlap = nfft ÷ 2, workspace = workspace)
        n_used = nseg
        if isempty(freqs)
            app.psd_freqs[j][] = Float64[]
            app.psd_values[j][] = Float64[]
        else
            f_plot = freqs[2:end]
            p_plot = psd[2:end]
            app.psd_freqs[j][] = f_plot
            app.psd_values[j][] = p_plot
            _autoscale_psd!(app.psd_axes[j], p_plot)
            if !isempty(f_plot)
                xlims!(app.psd_axes[j], first(f_plot), last(f_plot))
            end
        end
    end
    n_windows = _welch_segment_count(seg_lengths, nfft)
    app.psd_header[] = n_windows == 0 ?
        "No spectra: no unmasked stretch in this window is as long as one segment " *
        "(nfft = $(nfft), $(_fmt_dur(nfft / fs))). Widen the Window or unmask something." :
        _spectra_details_text(nfft, fs; n_windows = n_windows, n_stretches = n_used, n_runs = n_runs)
    return app
end

function _recompute_spectra!(app::TKApp)
    if _shows_spectra(app.view_mode[])
        _compute_psd_for_window!(app)
    else
        app.psd_header[] = _spectra_info_idle()
        return app
    end
    _refresh_status!(app)
    return app
end

function _schedule_spectra_recompute!(app::TKApp; delay_seconds::Real = 0.25)
    _shows_spectra(app.view_mode[]) || return app
    pending = app.spectra_timer[]
    if pending !== nothing
        try
            close(pending)
        catch
        end
    end
    app.spectra_timer[] = Timer(delay_seconds) do _
        try
            _recompute_spectra!(app)
        catch err
            @warn "Spectra recompute failed" exception = err
        end
    end
    return app
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
    elseif et === Makie.MouseEventTypes.rightclick
        if s.app.selection_visible[]
            _apply_selection_mask!(s.app, true)
            return Consume(true)
        end
    end
    return Consume(false)
end

"""
    _set_spectral_pin!(app, f)

Pin `f` as the comb filter fundamental. Fills in guides at 2f, 3f, ... up to
Nyquist so harmonics of the pinned tone can be told apart from unrelated peaks.
The guides are computed once here rather than lifted off the axis limits: lines
outside the current view simply do not render, and the pin then costs nothing
per frame.
"""
function _set_spectral_pin!(app::TKApp, f::Real)
    (isfinite(f) && f > 0) || return app
    _, fs = _current_nfft(app)
    nyq = fs / 2
    harmonics = Float64[]
    k = 2
    while k * f <= nyq && length(harmonics) < TK_MAX_HARMONICS
        push!(harmonics, k * f)
        k += 1
    end
    app.pin_freq[] = Float64[f]
    app.pin_harmonics[] = harmonics
    app.pin_text[] = _pin_label_text(f)
    _refresh_status!(app)
    return app
end

function _clear_spectral_pin!(app::TKApp)
    isempty(app.pin_freq[]) && return app
    app.pin_freq[] = Float64[]
    app.pin_harmonics[] = Float64[]
    app.pin_text[] = ""
    _refresh_status!(app)
    return app
end

"""
    _set_cursor_text!(app, channel, text)

Show `text` in spectral panel `channel` and blank the rest. Panels already blank
are left alone, so a cursor moving inside one panel never notifies the others.
"""
function _set_cursor_text!(app::TKApp, channel::Integer, text::AbstractString)
    for (k, obs) in enumerate(app.cursor_text)
        if k == channel
            obs[] = text
        elseif !isempty(obs[])
            obs[] = ""
        end
    end
    return app
end

"""
    _release_spectral_cursor!(app, channel)

Hide the hover cursor, but only if `channel` still owns it. Each Axis runs its
own mouse state machine, so leaving one panel and entering the next arrive in an
unspecified order; checking ownership makes the result order-independent.
"""
function _release_spectral_cursor!(app::TKApp, channel::Integer)
    app.cursor_panel[] == channel || return app
    app.cursor_panel[] = 0
    app.cursor_freq[] = Float64[]
    _set_cursor_text!(app, 0, "")
    return app
end

"""
    _cursor_data_position(ax, event) -> Union{Nothing, Point2d}

Mouse position in data coordinates. `MouseEvent.data` is in the axis' own
*transformed* space, so on these log-scaled spectral axes it carries log10(f)
rather than f; undo the transform the way Makie's rectangle zoom does. Returns
`nothing` when the result cannot be trusted, which also catches a future Makie
handing us data that is already untransformed.
"""
function _cursor_data_position(ax::Axis, event::Makie.MouseEvent)
    itf = Makie.inverse_transform(Makie.transform_func(ax))
    itf === nothing && return nothing
    p = Makie.apply_transform(itf, event.data)
    (isfinite(p[1]) && isfinite(p[2])) || return nothing
    lims = ax.finallimits[]
    lo = minimum(lims)
    hi = maximum(lims)
    # Generous slack: a decade outside the view either way is still plausibly the
    # user's pointer, ten decades out means we misread the coordinate space.
    @inbounds for d in 1:2
        span = hi[d] - lo[d]
        (p[d] < lo[d] - 10 * span || p[d] > hi[d] + 10 * span) && return nothing
    end
    return p
end

"""
    SpectralCursor

Hover readout and pinning for one PSD panel. Carries the channel the panel
belongs to and the bin last reported, so a mouse move that stays inside the same
frequency bin - which is nearly all of them - returns before touching a single
Observable.
"""
mutable struct SpectralCursor
    app::TKApp
    channel::Int
    last_bin::Int
end

function Makie.process_interaction(c::SpectralCursor, event::Makie.MouseEvent, ax::Axis)
    app = c.app
    et = event.type
    # a switched-off panel reads nothing and pins nothing
    app.channel_on[c.channel] || et === Makie.MouseEventTypes.out || return Consume(false)
    if et === Makie.MouseEventTypes.out
        _release_spectral_cursor!(app, c.channel)
        c.last_bin = 0
        return Consume(false)
    elseif et === Makie.MouseEventTypes.rightclick
        _clear_spectral_pin!(app)
        return Consume(true)
    elseif !(et === Makie.MouseEventTypes.over || et === Makie.MouseEventTypes.enter ||
             et === Makie.MouseEventTypes.leftclick)
        return Consume(false)
    end

    p = _cursor_data_position(ax, event)
    p === nothing && return Consume(false)
    f_cursor = p[1]
    f_cursor > 0 || return Consume(false)

    freqs = app.psd_freqs[c.channel][]
    psd = app.psd_values[c.channel][]
    (isempty(freqs) || length(psd) < length(freqs)) && return Consume(false)

    h = _snap_halfwidth(ax, 1, TK_SNAP_PIXELS)          # in decades, the axis is log10
    near, best = _peak_window(freqs, psd, f_cursor, f_cursor * exp10(-h), f_cursor * exp10(h))
    idx = (best != near && psd[best] > TK_SNAP_GAIN * psd[near]) ? best : near

    if et === Makie.MouseEventTypes.leftclick
        # ctrl+leftclick is Makie's own limit reset; leave it alone.
        Makie.ispressed(ax.scene, Makie.Keyboard.left_control) && return Consume(false)
        _set_spectral_pin!(app, freqs[idx])
        return Consume(true)
    end

    idx == c.last_bin && app.cursor_panel[] == c.channel && return Consume(false)
    c.last_bin = idx
    app.cursor_panel[] = c.channel
    app.cursor_freq[] = Float64[freqs[idx]]
    _set_cursor_text!(app, c.channel, _psd_cursor_text(freqs[idx], psd[idx]))
    return Consume(false)
end

function _clear_psd_axes!(app::TKApp)
    for ax in app.psd_axes
        delete!(ax)
    end
    empty!(app.psd_axes)
    empty!(app.psd_freqs)
    empty!(app.psd_values)
    empty!(app.cursor_text)
    app.psd_col = 0
    app.cursor_panel[] = 0
    app.cursor_freq[] = Float64[]
    # The pin deliberately survives: switching between the time and spectra
    # views rebuilds these axes, and carrying a candidate comb frequency
    # across that switch is the whole point of pinning it.
    return app
end

const TK_OFF_BG = RGBf(0.925, 0.93, 0.935)

"""
    _reset_channel_switches!(app, ta)

Switch every channel of a newly loaded record on, except one the logger did not
record - see [`_channel_recorded`](@ref).
"""
function _reset_channel_switches!(app::TKApp, ta::TimeArray)
    vals = _ta_values(ta)
    app.channel_on = [_channel_recorded(view(vals, :, j)) for j in axes(vals, 2)]
    return app
end

"""
    _channel_recorded(column) -> Bool

Whether a channel holds a recording. An input nothing was connected to comes
back as all `NaN`, or as one constant - loggers, and the LEMI-424 and GEOMAG
writers, fill it with zeros - and no real sensor reads exactly constant. A
one-sample record counts as recorded if that sample is finite.
"""
function _channel_recorded(column::AbstractVector)
    length(column) < 2 && return any(isfinite, column)
    ref = NaN
    for v in column
        isfinite(v) || continue
        isnan(ref) ? (ref = v) : (v != ref && return true)
    end
    return false
end

"""
    _apply_channel_switch!(app, j)

Show or grey out channel `j`. Its panels keep their place in the layout: a
switched-off channel draws on a grey background with a muted y axis, and its
plots - traces, mask spans, spectrum, cursor guides - are hidden. The data and
the mask are untouched.
"""
function _apply_channel_switch!(app::TKApp, j::Integer)
    on = app.channel_on[j]
    for ax in (get(app.axes, j, nothing), get(app.psd_axes, j, nothing))
        ax === nothing && continue
        ax.backgroundcolor = on ? TK_PANEL_BG : TK_OFF_BG
        ax.ylabelcolor = on ? TK_BLACK : TK_MUTED
        # with nothing drawn the y ticks would only label placeholder limits
        ax.yticklabelsvisible = on
        ax.yticksvisible = on
        for plt in ax.scene.plots
            plt.visible = on
        end
    end
    return app
end

# Names and units of the loaded channels, as the panel builders label them.
function _channel_labels(app::TKApp)
    names = _ta_colnames(app.data)
    metadata = _ta_meta(app.data)
    units_map = metadata isa AbstractDict ? get(metadata, :units, Dict{Symbol, String}()) : Dict{Symbol, String}()
    return names, units_map
end

function _clear_time_axes!(app::TKApp)
    foreach(delete!, app.axes)
    empty!(app.axes)
    empty!(app.line_clean)
    empty!(app.line_masked)
    return app
end

function _clear_spectra_placeholder!(app::TKApp)
    app.spectra_placeholder === nothing || delete!(app.spectra_placeholder)
    app.spectra_placeholder = nothing
    return app
end

"""
    _build_time_axes!(app)

One trace panel per channel in column 2 of the plot grid, x-linked.
"""
function _build_time_axes!(app::TKApp)
    names, units_map = _channel_labels(app)
    secs = app.time_seconds
    n = length(names)
    time_col = 2
    # the shared x buffer still holds the last traces; the new, empty y buffers
    # have to start from the same length
    app.line_x[] = Float64[]
    for (i, name) in enumerate(names)
        is_last = i == n
        unit_str = get(units_map, _symbolize(name), component_units(_symbolize(name)))
        col = _component_color(name)
        ax = Axis(app.plot_layout[i, time_col];
            ylabel = "$(_display_label(name)) [$(unit_str)]",
            ylabelrotation = pi / 2,
            ylabelpadding = 8.0,
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

        clean_obs = Observable{Vector{Float32}}(Float32[])
        masked_obs = Observable{Vector{Float32}}(Float32[])
        push!(app.line_clean, clean_obs)
        push!(app.line_masked, masked_obs)

        anchor = isempty(secs) ? 0.0 : first(secs)

        mask_lows_padded = lift(ls -> isempty(ls) ? Float64[anchor] : ls, app.mask_lows)
        mask_highs_padded = lift(hs -> isempty(hs) ? Float64[anchor] : hs, app.mask_highs)
        vspan!(ax, mask_lows_padded, mask_highs_padded; color = TK_MASK_FILL)

        sel_lows = lift((vis, sel) -> vis ? Float64[sel[1]] : Float64[anchor], app.selection_visible, app.selection)
        sel_highs = lift((vis, sel) -> vis ? Float64[sel[2]] : Float64[anchor], app.selection_visible, app.selection)
        vspan!(ax, sel_lows, sel_highs; color = TK_SEL_FILL, strokecolor = TK_SEL_EDGE, strokewidth = 0.8)

        if length(secs) == 1
            scatter!(ax, app.line_x, clean_obs; color = col, markersize = 4)
            scatter!(ax, app.line_x, masked_obs; color = TK_MUTED, markersize = 4)
        else
            lines!(ax, app.line_x, clean_obs; color = col, linewidth = 1.4, joinstyle = :round)
            lines!(ax, app.line_x, masked_obs; color = TK_MUTED, linewidth = 1.4, joinstyle = :round)
        end

        deregister_interaction!(ax, :rectanglezoom)
        register_interaction!(ax, :tk_select, DragSelect(app, false, 0.0))

        push!(app.axes, ax)
    end
    length(app.axes) > 1 && linkxaxes!(app.axes...)
    return app
end

"""
    _build_psd_axes!(app, psd_col)

One spectrum panel per channel in column `psd_col` of the plot grid: 3 beside
the traces, 2 when the spectra have the window to themselves.
"""
function _build_psd_axes!(app::TKApp, psd_col::Integer)
    names, units_map = _channel_labels(app)
    n = length(names)
    time_on = psd_col == 3
    for (i, name) in enumerate(names)
        is_last = i == n
        unit_str = get(units_map, _symbolize(name), component_units(_symbolize(name)))
        col = _component_color(name)
        ax_psd = Axis(app.plot_layout[i, psd_col];
            ylabel = time_on ? "" : "$(_display_label(name)) [$(unit_str)]",
            ylabelrotation = pi / 2,
            ylabelpadding = 8.0,
            xscale = log10,
            yscale = log10,
            xticklabelsvisible = is_last,
            xticksvisible = is_last,
            xgridvisible = false,
            ygridvisible = false,
            topspinevisible = false,
            rightspinevisible = false,
            yticks = _decade_ticks(3),
            xticks = _decade_ticks(4),
            spinewidth = 0.9,
            bottomspinecolor = TK_FRAME,
            leftspinecolor = TK_FRAME,
            xtickcolor = TK_FRAME,
            ytickcolor = TK_FRAME,
            xticklabelcolor = TK_BLACK,
            yticklabelcolor = TK_BLACK,
            backgroundcolor = TK_PANEL_BG,
            tellheight = false,
            tellwidth = false,
        )
        if !is_last
            hidexdecorations!(ax_psd; ticks = true, ticklabels = true, grid = false)
        end
        freqs_obs = Observable{Vector{Float64}}(Float64[])
        psd_obs = Observable{Vector{Float64}}(Float64[])
        push!(app.psd_freqs, freqs_obs)
        push!(app.psd_values, psd_obs)
        lines!(ax_psd, freqs_obs, psd_obs; color = col, linewidth = 1.6)

        # Cursor, pin and harmonic guides. All of them opt out of the axis
        # autolimits: VLines reports the extrema of its own input as data
        # limits, which on a log axis would let a pinned low frequency drag
        # the x range open on the next reset - and would throw from extrema
        # the moment the vector is empty.
        vlines!(ax_psd, app.pin_harmonics; color = TK_PIN_GUIDE, linewidth = 0.8,
            linestyle = :dot, xautolimits = false, yautolimits = false, inspectable = false)
        vlines!(ax_psd, app.pin_freq; color = TK_PIN_LINE, linewidth = 1.4,
            xautolimits = false, yautolimits = false, inspectable = false)
        vlines!(ax_psd, app.cursor_freq; color = TK_CURSOR_LINE, linewidth = 0.9,
            linestyle = :dash, xautolimits = false, yautolimits = false, inspectable = false)

        # Anchored in relative space, so the readout keeps its corner when the
        # limits move and never counts towards the autolimits either.
        cursor_obs = Observable("")
        push!(app.cursor_text, cursor_obs)
        text!(ax_psd, Point2f(0.985, 0.96); text = cursor_obs, space = :relative,
            align = (:right, :top), color = TK_BLACK, inspectable = false)
        if i == 1
            text!(ax_psd, Point2f(0.015, 0.96); text = app.pin_text, space = :relative,
                align = (:left, :top), color = TK_PIN_LINE, inspectable = false)
        end

        deregister_interaction!(ax_psd, :rectanglezoom)
        register_interaction!(ax_psd, :tk_cursor, SpectralCursor(app, i, 0))
        push!(app.psd_axes, ax_psd)
    end
    length(app.psd_axes) > 1 && linkxaxes!(app.psd_axes...)
    app.psd_col = psd_col
    return app
end

"""
    _sync_view_axes!(app; with_spectra = true)

Bring the panels in line with `app.view_mode`, touching only what changes:
going from Time to Time | Spectra keeps the trace panels and adds the spectra,
and back again only removes them. Creating panels is what costs time, so a
view change should not rebuild the ones it keeps.

With `with_spectra = false` the spectra column gets a placeholder instead of
its panels, so the new layout can be drawn before they are built.
"""
function _sync_view_axes!(app::TKApp; with_spectra::Bool = true)
    mode = app.view_mode[]
    want_time, want_psd = _shows_time(mode), _shows_spectra(mode)
    # Column 1 holds the channel switches. Beside traces the spectra take a 40%
    # strip in column 3; on their own they move to column 2 and take the rest.
    psd_col = want_time ? 3 : 2
    n = length(_ta_colnames(app.data))

    want_time || _clear_time_axes!(app)
    keep_psd = want_psd && app.psd_col == psd_col && !isempty(app.psd_axes)
    keep_psd || _clear_psd_axes!(app)
    want_psd || _clear_spectra_placeholder!(app)

    want_time && isempty(app.axes) && _build_time_axes!(app)
    if want_psd && isempty(app.psd_axes)
        if with_spectra
            _clear_spectra_placeholder!(app)
            _build_psd_axes!(app, psd_col)
        elseif app.spectra_placeholder === nothing
            app.spectra_placeholder = Label(app.plot_layout[1:max(n, 1), psd_col], "Computing spectra…";
                                            color = TK_GREY, tellwidth = false, tellheight = false)
        end
    end

    n > 1 && rowgap!(app.plot_layout, 8)
    colsize!(app.plot_layout, 1, Auto(true))            # as wide as a checkbox
    if want_time && want_psd
        colsize!(app.plot_layout, 2, Auto(true, 0.6))
        colsize!(app.plot_layout, 3, Auto(true, 0.4))
    else
        try
            trim!(app.plot_layout)
        catch
        end
        colsize!(app.plot_layout, 2, Auto(true, 1.0))
    end
    colgap!(app.plot_layout, 12)
    foreach(j -> _apply_channel_switch!(app, j), 1:n)
    return app
end

"""
    _build_axes!(app, ta)

Full rebuild for a newly loaded record: time base, channel switches and the
panels of the current view. A change of view goes through
[`_switch_view!`](@ref) instead.
"""
function _build_axes!(app::TKApp, ta::TimeArray)
    _clear_time_axes!(app)
    _clear_psd_axes!(app)
    _clear_spectra_placeholder!(app)
    foreach(delete!, app.channel_boxes)
    empty!(app.channel_boxes)

    names = _ta_colnames(ta)
    vals = _ta_values(ta)
    times = _ensure_datetime(_ta_timestamps(ta))
    origin = first(times)
    secs = _seconds_since(origin, times)
    app.origin = origin
    app.time_seconds = secs
    app.span_seconds = isempty(secs) ? 0.0 : last(secs) - first(secs)
    app.raw_values = Matrix{Float64}(vals)
    app.line_x[] = Float64[]

    n = length(names)
    length(app.channel_on) == n || _reset_channel_switches!(app, ta)
    for i in 1:n
        # tellheight = false: the axes, not the switch, decide how tall a row is
        box = Checkbox(app.plot_layout[i, 1]; checked = app.channel_on[i], tellheight = false)
        on(box.checked) do checked
            app.channel_on[i] = checked
            _apply_channel_switch!(app, i)
            # autoscaling skipped the channel while it was off
            checked && _autoscale_y!(app, _visible_x_window(app)...)
            _recompute_spectra!(app)
        end
        push!(app.channel_boxes, box)
    end
    _sync_view_axes!(app)
    return app.axes
end

"""
    _switch_view!(app, mode)

Change the view without keeping the menu waiting. Panels the new view keeps
stay; the traces are redrawn at once; spectrum panels that have to be created
are built on the next tick, after a frame showing a placeholder in their place.
Picking another view before then drops that build.
"""
function _switch_view!(app::TKApp, mode::Symbol)
    app.view_mode[] = mode
    app.view_generation += 1
    generation = app.view_generation
    _sync_view_axes!(app; with_spectra = false)
    _refresh_mask_overlay!(app)
    _update_x_window!(app)
    if _shows_spectra(mode) && isempty(app.psd_axes)
        app.psd_header[] = "Computing spectra…"
        @async try
            sleep(1 / 20)                              # let the frame without spectra show
            app.view_generation == generation || return
            _sync_view_axes!(app)
            _recompute_spectra!(app)
            _refresh_status!(app)
        catch err
            @warn "Could not build the spectra panels" exception = err
        end
    else
        _recompute_spectra!(app)
    end
    return app
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

function _page_window!(app::TKApp, direction::Integer)
    ws = app.window_seconds[]
    span = app.span_seconds
    visible = isfinite(ws) ? min(ws, max(span, 1.0)) : max(span, 1.0)
    max_start = max(0.0, span - visible)
    max_start <= 0 && return app
    new_start = clamp(app.window_start[] + direction * visible, 0.0, max_start)
    set_close_to!(app.slider, new_start)
    return app
end

"""
    _show_record!(app, ta, mask)

Put `ta` on screen with `mask`, or a fresh one with the `NaN` rows masked when
`mask` is `nothing`. The channel switches carry over when the channels are the
same, as they are across the rates of one site.
"""
function _show_record!(app::TKApp, ta::TimeArray, mask::Union{Nothing, TimekeeperMask})
    same_channels = _ta_colnames(ta) == _ta_colnames(app.data)
    app.data = ta
    app.mask = mask === nothing ? TimekeeperMask(ta) : mask
    app.selection_visible[] = false
    # A pin from the previous record may sit above the new Nyquist, so drop it.
    _clear_spectral_pin!(app)
    app.window_start[] = 0.0
    same_channels || _reset_channel_switches!(app, ta)
    _build_axes!(app, ta)
    mask === nothing && _auto_mask_nan!(app.mask, app.raw_values)
    _refresh_slider_range!(app)
    app.summary_label.text[] = _summary_text(ta)
    _refresh_mask_overlay!(app)
    _update_x_window!(app)
    return app
end

"""
    _apply_loaded_data!(app, ta, fmt, source_path; site_rates = Float64[], rate = nothing)

Show a newly loaded record. For a Metronix site `site_rates` lists every rate
it holds and `rate` the one `ta` carries; the rate menu offers them all.
Anything else offers the one rate `ta` has.
"""
function _apply_loaded_data!(app::TKApp, ta::TimeArray, fmt::Symbol, source_path::AbstractString;
                             site_rates::Vector{Float64} = Float64[], rate = nothing)
    app.source_format = fmt
    app.source_path = String(source_path)
    _reset_channel_switches!(app, ta)
    _show_record!(app, ta, nothing)
    empty!(app.rate_intervals)
    if isempty(site_rates)
        app.site_rates = [_sample_rate_from_timearray(ta)]
        app.rate_index = 1
    else
        app.site_rates = site_rates
        app.rate_index = something(findfirst(r -> isapprox(r, something(rate, first(site_rates))), site_rates), 1)
    end
    _refill_rate_menu!(app)
    return app
end

function _refill_rate_menu!(app::TKApp)
    app.rate_menu_updating = true
    try
        app.rate_menu.options[] = [(_format_fs(r), i) for (i, r) in enumerate(app.site_rates)]
        app.rate_menu.i_selected[] = app.rate_index
    finally
        app.rate_menu_updating = false
    end
    return app
end

"""
    _select_rate!(app, i)

Load rate `i` of the site in place of the one on screen; only one rate is held
in memory. The rate on screen leaves its masked intervals behind, and a rate
that had some gets them back. If the load fails, the record on screen stays.
"""
function _select_rate!(app::TKApp, i::Integer)
    i == app.rate_index && return app
    rate = app.site_rates[i]
    site_dir = app.source_path
    console = _show_progress_window()
    app.status_label.text[] = "Loading $(_format_fs(rate))…"
    @async try
        ta, _ = _run_with_progress_pump(console) do
            _load_metronix_site(site_dir; progress = console, rate = rate)
        end
        _progress_note!(console, :info, "Drawing the record")
        _flush_progress!(console)
        app.rate_intervals[app.site_rates[app.rate_index]] = copy(app.mask.intervals)
        app.rate_index = i
        _show_record!(app, ta, nothing)
        restored = pop!(app.rate_intervals, rate, Tuple{DateTime, DateTime}[])
        if !isempty(restored)
            foreach(((a, b),) -> mask_interval!(app.mask, a, b), restored)
            _refresh_mask_overlay!(app)
        end
        app.status_label.text[] = "Showing $(_format_fs(rate))"
        _progress_finish!(console, :ok, "Showing $(_format_fs(rate))" *
            (isempty(restored) ? "" : ", with its $(length(restored)) masked interval" *
                                      (length(restored) == 1 ? "" : "s")))
        _close_progress_window!(console, 1.5)
    catch err
        @warn "Could not load the $(_format_fs(rate)) runs" exception = err
        msg = sprint(showerror, err)
        _progress_finish!(console, :error, "Could not load $(_format_fs(rate))\n" * msg)
        _refill_rate_menu!(app)                        # back to the rate still on screen
        app.status_label.text[] = "Load failed: $(msg)"
    end
    return app
end

function _load_site_any(dir::AbstractString;
                       progress::Union{ProgressConsole, Nothing} = nothing,
                       rate::Union{Nothing, Real} = nothing)
    # A site, or a single meas_ directory - which can itself hold several
    # runs at several rates - loads one rate at a time.
    return is_metronix_site(dir) ? _load_metronix_site(dir; progress = progress, rate = rate) :
           _load_site_directory(dir; progress = progress)
end

"""
    _load_any_path(path) -> (TimeArray, Symbol)

Load whatever `path` points at: a data file, a Metronix `meas_*` run, a site
directory of runs, or a Metronix site.
"""
_load_any_path(path::AbstractString) =
    isdir(path) ? _load_site_any(path) : _load_run_any(path)

function _ext_for_format(fmt::Symbol)
    fmt === :lemi_xyz && return ".xyz"
    return ".txt"
end

function TKApp(
    ta::TimeArray;
    size = (1600, 900),
    source_format::Symbol = :lemi424,
    source_path::AbstractString = "",
)
    GLMakie.activate!(title = "Timekeepers")
    fig = Figure(; size = size, backgroundcolor = :white,
        figure_padding = (14, 14, 8, 8))

    toolbar = GridLayout(fig[1, 1]; tellheight = true)
    rate_menu = _logo_menu(toolbar[1, 1];
        options = [(_format_fs(_sample_rate_from_timearray(ta)), 1)], width = 110)
    summary_label = Label(toolbar[1, 2], _summary_text(ta);
        color = TK_GREY, halign = :left, tellwidth = false)

    actions = GridLayout(toolbar[1, 3]; halign = :right)
    load_btn = Button(actions[1, 1]; label = "Load Run…")
    load_site_btn = Button(actions[1, 2]; label = "Load Site…")
    mask_btn = Button(actions[1, 3]; label = "Mask")
    unmask_btn = Button(actions[1, 4]; label = "Unmask")
    clear_btn = Button(actions[1, 5]; label = "Clear")
    write_btn = Button(actions[1, 6]; label = "Write")
    view_label = Label(actions[1, 7], "View:"; color = TK_GREY)
    view_menu = _logo_menu(actions[1, 8]; options = VIEW_OPTIONS, default = "Time", width = 170)
    window_label = Label(actions[1, 9], "Window:"; color = TK_GREY)
    window_box = Textbox(actions[1, 10]; stored_string = "1", width = 56,
        restriction = isdigit, validator = s -> _parse_window_count(s) !== nothing,
        textcolor = TK_BLACK, halign = :right)
    window_menu = _logo_menu(actions[1, 11]; options = WINDOW_UNITS, default = "All", width = 100)
    colgap!(actions, 6)

    # Size the controls to their content, so buttons never clip their labels;
    # the summary between them takes whatever width is left.
    colsize!(toolbar, 1, Auto(true))
    colsize!(toolbar, 2, Auto(true, 1.0))
    colsize!(toolbar, 3, Auto(true))

    plot_layout = GridLayout(fig[2, 1]; tellheight = false)

    psd_header_obs = Observable(_spectra_info_idle())
    help_visible_obs = Observable(false)

    info_row = GridLayout(fig[3, 1])
    info_btn = Button(info_row[1, 1]; label = "i")
    Label(info_row[1, 2],
        lift((line, help) -> help ? _spectra_help_text() : line, psd_header_obs, help_visible_obs);
        color = TK_BLACK, halign = :left, tellwidth = false)
    colsize!(info_row, 2, Auto(true, 1.0))
    colgap!(info_row, 8)

    slider_grid = GridLayout(fig[4, 1]; tellheight = true)
    Label(slider_grid[1, 1], "Scroll"; color = TK_GREY, halign = :right)
    prev_btn = Button(slider_grid[1, 2]; label = "<")
    slider = Slider(slider_grid[1, 3]; range = 0.0:1.0:0.0, startvalue = 0.0,
        linewidth = 11.0,
        color_inactive = TK_CTRL_TRACK,
        color_active_dimmed = RGBAf(TK_CTRL_ACCENT.r, TK_CTRL_ACCENT.g, TK_CTRL_ACCENT.b, 0.45),
        color_active = TK_CTRL_ACCENT)
    next_btn = Button(slider_grid[1, 4]; label = ">")
    colsize!(slider_grid, 3, Auto(true, 1.0))
    colgap!(slider_grid, 10)

    status_label = Label(fig[5, 1], "";
        color = TK_GREY, halign = :left, tellwidth = false)

    # the plots take whatever height the control rows, each fitted to its
    # contents, leave over
    rowsize!(fig.layout, 2, Auto(true, 1.0))
    rowgap!(fig.layout, 6)

    selection = Observable((0.0, 0.0))
    selection_visible = Observable(false)
    mask_lows = Observable(Float64[])
    mask_highs = Observable(Float64[])
    window_seconds_obs = Observable(Inf)
    window_start_obs = Observable(0.0)
    line_x_obs = Observable(Float64[])
    view_mode_obs = Observable(:time)
    cursor_freq_obs = Observable(Float64[])
    pin_freq_obs = Observable(Float64[])
    pin_harmonics_obs = Observable(Float64[])
    pin_text_obs = Observable("")

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
        Float64[],
        window_seconds_obs,
        window_start_obs,
        slider,
        window_menu,
        window_box,
        selection,
        selection_visible,
        mask_lows,
        mask_highs,
        source_format,
        String(source_path),
        view_mode_obs,
        view_menu,
        Axis[],
        Observable{Vector{Float64}}[],
        Observable{Vector{Float64}}[],
        psd_header_obs,
        help_visible_obs,
        Dict{Tuple{Int, Float64, Int, Symbol, Symbol}, TKSpectralWorkspace}(),
        Ref{Union{Nothing, Timer}}(nothing),
        Observable{String}[],
        Ref(0),
        cursor_freq_obs,
        pin_freq_obs,
        pin_harmonics_obs,
        pin_text_obs,
        Bool[],
        Checkbox[],
        0,
        nothing,
        0,
        [_sample_rate_from_timearray(ta)],
        1,
        Dict{Float64, Vector{Tuple{DateTime, DateTime}}}(),
        rate_menu,
        false,
    )

    _reset_channel_switches!(app, ta)
    _build_axes!(app, ta)
    _auto_mask_nan!(app.mask, app.raw_values)
    _refresh_mask_overlay!(app)
    _refresh_slider_range!(app)
    _update_x_window!(app)
    app.status_label.text[] = _ready_status_text(app)

    on(slider.value) do v
        app.window_start[] = Float64(v)
        _update_x_window!(app)
        _schedule_spectra_recompute!(app)
    end
    on(prev_btn.clicks) do _
        _page_window!(app, -1)
    end
    on(next_btn.clicks) do _
        _page_window!(app, +1)
    end
    # The span is the typed count times the unit; either one changing
    # applies it. The box only accepts digits and takes a new count on Enter.
    apply_window = function (_)
        unit = window_menu.selection[]
        count = _parse_window_count(window_box.stored_string[])
        (unit === nothing || count === nothing) && return
        app.window_seconds[] = _window_span(count, unit)
        _refresh_slider_range!(app)
        _update_x_window!(app)
        _recompute_spectra!(app)
    end
    on(apply_window, window_menu.selection)
    on(apply_window, window_box.stored_string)
    on(info_btn.clicks) do _
        app.help_visible[] = !app.help_visible[]
    end
    on(view_menu.selection) do mode
        mode === nothing && return
        Symbol(mode) === app.view_mode[] && return
        _switch_view!(app, Symbol(mode))
    end
    on(rate_menu.selection) do i
        (i === nothing || app.rate_menu_updating) && return
        _select_rate!(app, i)
    end

    on(load_btn.clicks) do _
        path = ""
        try
            path = pick_file(; filterlist = "txt,dat,lem,xyz,ats,xml")
        catch err
            @warn "Could not open file picker" exception = err
            return
        end
        isempty(path) && return
        # A meas_ directory can hold several Metronix runs: the user opens it
        # and picks a run's .xml (or one of its .ats files), and that run -
        # its channels plus the XML - is what loads.
        target = path
        app.status_label.text[] = "Loading $(basename(target))…"
        @async begin
            try
                ta, fmt = fetch(Threads.@spawn(_load_run_any(target)))
                _apply_loaded_data!(app, ta, fmt, target)
                app.status_label.text[] = _ready_status_text(app)
            catch err
                err isa TaskFailedException && (err = err.task.exception)
                @warn "Could not load $target" exception = err
                app.status_label.text[] = "Load failed: $(sprint(showerror, err))"
            end
        end
    end
    on(load_site_btn.clicks) do _
        dir = ""
        try
            dir = pick_folder()
        catch err
            @warn "Could not open folder picker" exception = err
            return
        end
        isempty(dir) && return
        site_name = _site_name_from_dir(dir)
        app.status_label.text[] = "Loading site $(site_name)…"
        @async begin
            console = nothing
            try
                # one rate at a time: a mixed-rate site asks which to import
                ok, rate = _prompt_metronix_rate(dir)
                if !ok
                    app.status_label.text[] = "Load cancelled"
                    return
                end
                console = _show_progress_window()
                rate === nothing || _progress_note!(console, :info, "Loading the $(_format_fs(rate)) runs")
                ta, fmt = _run_with_progress_pump(console) do
                    _load_site_any(dir; progress = console, rate = rate)
                end
                _progress_note!(console, :info, "Drawing the record")
                _flush_progress!(console)
                # a raw site was read from its <site>.TK copy; switch rates there too
                md = _ta_meta(ta)
                src = md isa AbstractDict ? String(get(md, :metronix_source_dir, dir)) : dir
                _apply_loaded_data!(app, ta, fmt, src;
                    site_rates = rate === nothing ? Float64[] : metronix_site_rates(src), rate = rate)
                app.status_label.text[] = "Loaded site $(site_name)"
                _progress_finish!(console, :ok, "Loaded $(site_name)")
                _close_progress_window!(console, 3.0)
            catch err
                @warn "Site load failed" dir exception = err
                msg = sprint(showerror, err)
                try
                    _progress_finish!(console, :error, "Could not load $(site_name)\n" * msg)
                catch
                end
                try
                    app.status_label.text[] = "Site load failed: $(msg)"
                catch
                end
            end
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
        md = _ta_meta(app.data)
        if md isa AbstractDict && get(md, :source_format, nothing) === :metronix &&
           haskey(md, :site_dir)
            site_dir = String(md[:site_dir])
            # The whole site is written, every rate: the rate on screen cut by
            # its mask, each rate masked before switching away by its own
            # intervals, and the rest copied as they are.
            cuts = Dict{Float64, Vector{Tuple{DateTime, DateTime}}}(
                r => copy(ivs) for (r, ivs) in app.rate_intervals)
            cuts[app.site_rates[app.rate_index]] = copy(app.mask.intervals)
            n_cut = sum(length, values(cuts); init = 0)
            @async begin
                console = _show_progress_window()
                try
                    _progress_note!(console, :info, "Writing all of $(_site_name_from_dir(site_dir))")
                    _progress_println!(console, String(site_dir))
                    _progress_step!(console, 1, 1,
                        "Copying every run, cutting $(n_cut) masked interval" * (n_cut == 1 ? "" : "s"))
                    dest = _run_with_progress_pump(console) do
                        write_metronix_site_masked(site_dir; rate_intervals = cuts)
                    end
                    _progress_finish!(console, :ok, "Wrote $(basename(dest))")
                    app.status_label.text[] = "Wrote $(basename(dest))"
                    _close_progress_window!(console, 3.0)
                catch err
                    @warn "Metronix write failed" exception = err
                    app.status_label.text[] = "Write failed: $(sprint(showerror, err))"
                    try
                        _progress_finish!(console, :error, "Could not write the site\n" * sprint(showerror, err))
                    catch
                    end
                end
            end
            return
        end
        if isdir(app.source_path)
            site_dir = rstrip(app.source_path, ['/', '\\'])
            stem = _site_name_from_dir(site_dir) * "_combined"
            ext = _ext_for_format(app.source_format)
            clean_path = joinpath(site_dir, stem * "_clean" * ext)
            mask_path = joinpath(site_dir, stem * "_mask.csv")
        else
            dir = dirname(app.source_path)
            stem, ext = splitext(basename(app.source_path))
            clean_path = joinpath(dir, stem * "_clean" * ext)
            mask_path = joinpath(dir, stem * "_mask.csv")
        end
        fmt = app.source_format
        try
            cleaned = cleaned_timearray(app; mode = :drop)
            write_mask(mask_path, app)
            app.status_label.text[] = "Writing $(basename(clean_path))…"
            @async begin
                try
                    fetch(Threads.@spawn _write_data_file(clean_path, cleaned, fmt))
                    @info "Wrote cleaned data and mask" clean_path mask_path
                    app.status_label.text[] = "Wrote $(basename(clean_path))  and  $(basename(mask_path))"
                catch err
                    err isa TaskFailedException && (err = err.task.exception)
                    @warn "Could not write outputs" exception = err
                    app.status_label.text[] = "Write failed: $(sprint(showerror, err))"
                end
            end
        catch err
            @warn "Could not write outputs" exception = err
            app.status_label.text[] = "Write failed: $(sprint(showerror, err))"
        end
    end

    return app
end

function TKApp(path::AbstractString; kwargs...)
    ta, fmt = _load_any_path(path)
    return TKApp(ta; source_format = fmt, source_path = path, kwargs...)
end

function TKApp(; kwargs...)
    return TKApp(_interactive_timearray(); kwargs...)
end

function _open_app_screen(app::TKApp; maximize::Bool)
    try
        screen = display(GLMakie.Screen(; title = "Timekeepers", visible = false,
                                          focus_on_show = true), app.figure)
        if maximize
            try
                GLMakie.GLFW.MaximizeWindow(screen.glscreen)
            catch err
                @warn "Could not maximize window" exception = err
            end
        end
        try
            GLMakie.Makie.colorbuffer(screen)
        catch
        end
        try
            GLMakie.GLFW.ShowWindow(screen.glscreen)
        catch
        end
        return screen
    catch err
        @warn "Could not open window hidden; opening directly" exception = err
        screen = display(app.figure)
        maximize && try
            GLMakie.GLFW.MaximizeWindow(screen.glscreen)
        catch
        end
        return screen
    end
end

function Base.display(app::TKApp)
    @info "Opening Timekeepers window"
    screen = _open_app_screen(app; maximize = false)
    app.status_label.text[] = _ready_status_text(app)
    @info "Timekeepers window is open and ready"
    return screen
end

"""
    run_tkapp(; kwargs...)
    run_tkapp(path::AbstractString; kwargs...)
    run_tkapp(ta::TimeArray; kwargs...)
    run_tkapp(app::TKApp)

Open the Timekeepers explorer window and block until it is closed, then return
the [`TKApp`](@ref) so the mask survives the session.

```julia
using Timekeepers
app = run_tkapp("data/LEMI090.txt")
segments = good_segments(app; min_samples = 256)
```

In the window: **Load Run…** / **Load Site…** to open data, the **Window** menu
and slider to scroll the record, left-drag to select an interval, then
**Mask** / **Unmask** / **Clear** to edit it and **Write** to export in the
source format. The **View** menu adds per-channel PSD panels.

Requires a desktop session with OpenGL 3.3 or newer; see [`TKApp`](@ref) to
build the app without displaying it.
"""
function run_tkapp(app::TKApp)
    @info "Opening Timekeepers window"
    screen = _open_app_screen(app; maximize = true)
    app.status_label.text[] = _ready_status_text(app)
    @info "Timekeepers window is open and ready"
    try
        wait(screen)
    catch
    end
    @info "Timekeepers window closed"
    return app
end

function run_tkapp(; kwargs...)
    @info "Starting Timekeepers"
    return run_tkapp(TKApp(; kwargs...))
end

function run_tkapp(path::AbstractString; kwargs...)
    @info "Starting Timekeepers" path
    return run_tkapp(TKApp(path; kwargs...))
end

function run_tkapp(ta::TimeArray; kwargs...)
    @info "Starting Timekeepers" samples = length(_ta_timestamps(ta))
    return run_tkapp(TKApp(ta; kwargs...))
end

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
