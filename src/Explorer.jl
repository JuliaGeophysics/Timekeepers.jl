# Explorer.jl - the interactive GLMakie application.
# Author: @pankajkmishra
#
# Defines TKApp and the whole UI: stacked per-component time series with
# optional PSD or spectrogram panels, a scrolling time window, drag-to-select
# with mask/unmask, and loading or writing single runs and whole sites.
#
# Two things keep it responsive on long records. Plotted series are decimated
# to a fixed bucket count by min/max per bucket, drawn into buffers the plot
# Observables already own so panning allocates almost nothing. Spectral
# recomputes are debounced behind a timer and reuse cached SpectralWorkspaces
# keyed by their configuration.

const TK_BLACK = RGBf(0.05, 0.05, 0.07)
const TK_BLUE = RGBf(0.114, 0.306, 0.847)
const TK_GREY = RGBf(0.42, 0.45, 0.50)
const TK_MUTED = RGBf(0.62, 0.66, 0.71)
const TK_FRAME = RGBf(0.55, 0.58, 0.63)
const TK_GRID = RGBAf(0.55, 0.58, 0.63, 0.18)
const TK_MINOR_GRID = RGBAf(0.55, 0.58, 0.63, 0.08)
const TK_PANEL_BG = RGBf(0.992, 0.993, 0.995)
const TK_MASK_FILL = RGBAf(0.40, 0.43, 0.50, 0.18)
const TK_SEL_FILL = RGBAf(0.114, 0.306, 0.847, 0.18)
const TK_SEL_EDGE = RGBAf(0.114, 0.306, 0.847, 0.70)

const TK_LOGO_SKY   = RGBf(0.76, 0.84, 0.87)
const TK_LOGO_SAGE  = RGBf(0.42, 0.63, 0.57)
const TK_LOGO_TEAL  = RGBf(0.18, 0.69, 0.78)
const TK_LOGO_SLATE = RGBf(0.30, 0.32, 0.34)
const TK_LOGO_CORAL = RGBf(0.90, 0.35, 0.33)

"""
Default colormap for spectral displays. Reversed so low power reads cool and
high power reads warm, the usual convention for a power spectrogram.
"""
const TK_SPECTRAL_COLORMAP = Makie.Reverse(:Spectral)

const TK_CTRL_BTN    = RGBf(0.88, 0.90, 0.93)
const TK_CTRL_ACCENT = RGBf(0.13, 0.55, 0.60)
const TK_CTRL_TRACK  = RGBAf(0.62, 0.66, 0.72, 0.35)

_shade(c::RGBf, amt::Real) = (f = clamp(1 - amt, 0.0, 1.0); RGBf(c.r * f, c.g * f, c.b * f))

function _logo_button(parent, label, color::RGBf; textcolor = :white, fontsize = 11,
                      font = :regular, height = nothing)
    return Button(parent;
        label = label, fontsize = fontsize, font = font, height = height,
        buttoncolor = color,
        buttoncolor_hover = _shade(color, 0.10),
        buttoncolor_active = _shade(color, 0.22),
        labelcolor = textcolor,
        labelcolor_hover = textcolor,
        labelcolor_active = textcolor,
    )
end

function _logo_menu(parent; kwargs...)
    return Menu(parent;
        fontsize = 11,
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

const VIEW_OPTIONS = [
    ("Time", :time),
    ("Time | Spectra", :time_spectra),
    ("Time | Spectrogram", :time_spectrogram),
]

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
    window_menu::Any
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
    psd_header::Any
    spec_axes::Vector{Axis}
    spec_times::Vector{Observable{Vector{Float64}}}
    spec_freqs::Vector{Observable{Vector{Float64}}}
    spec_matrix::Vector{Observable{Matrix{Float64}}}
    spectral_workspaces::Dict{Tuple{Int, Float64, Int, Symbol, Symbol}, TKSpectralWorkspace}
    spectra_timer::Base.RefValue{Union{Nothing, Timer}}
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

const TERM_BG = RGBAf(0.95, 0.95, 0.96, 1.0)
const TERM_FG = RGBf(0.10, 0.10, 0.12)
const TERM_TITLE = RGBf(0.0, 0.0, 0.0)
const TERM_DIM = RGBf(0.55, 0.55, 0.60)
const TERM_FONT = "Consolas"

mutable struct ProgressConsole
    screen::Any
    text_obs::Observable{String}
    lines::Vector{String}
    max_lines::Int
    lock::ReentrantLock
    dirty::Threads.Atomic{Bool}
end

_progress_println!(::Nothing, _msg::AbstractString) = nothing
function _progress_println!(p::ProgressConsole, msg::AbstractString)
    lock(p.lock) do
        for raw in split(String(msg), '\n')
            push!(p.lines, "  " * raw)
        end
        if length(p.lines) > p.max_lines
            deleteat!(p.lines, 1:(length(p.lines) - p.max_lines))
        end
    end
    p.dirty[] = true
    yield()
    return nothing
end

_flush_progress!(::Nothing) = nothing
function _flush_progress!(p::ProgressConsole)
    p.dirty[] || return nothing
    text = lock(() -> join(p.lines, "\n"), p.lock)
    p.dirty[] = false
    p.text_obs[] = text
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
    _progress_println!(progress, "> site      = $(site_name)")
    _progress_println!(progress, "> directory = $(dir)")
    _progress_println!(progress, "> found $(length(files)) file(s):")
    for f in files
        _progress_println!(progress, "    - $(basename(f))")
    end
    _progress_println!(progress, "")
    _progress_println!(progress, "> starting load...")

    loaded = Tuple{TimeArray, Symbol}[]
    skipped = 0
    for (i, path) in enumerate(files)
        _progress_println!(progress, "  [$i/$(length(files))] loading $(basename(path))")
        try
            ta, fmt = _load_data_file(path)
            push!(loaded, (ta, fmt))
        catch err
            @warn "Skipping file (could not load)" path exception = err
            _progress_println!(progress, "      ! skip (parse failed)")
            skipped += 1
        end
    end
    isempty(loaded) && error("Could not load any data files from: $dir")

    formats = unique(t[2] for t in loaded)
    length(formats) > 1 &&
        @warn "Mixed file formats in site; treating combined output as $(first(formats))" formats
    fmt = first(formats)

    _progress_println!(progress, "")
    _progress_println!(progress, "> sorting $(length(loaded)) runs by start time...")
    sort!(loaded; by = x -> first(_ensure_datetime(_ta_timestamps(x[1]))))

    _progress_println!(progress, "> combining $(length(loaded)) runs onto a single time grid...")
    combined = _combine_site_timearrays([t[1] for t in loaded], site_name)
    meta = _ta_meta(combined)
    @info "Loaded site" site = site_name n_files = length(loaded) skipped span = "$(meta[:start_time]) → $(meta[:end_time])"

    _progress_println!(progress, "")
    _progress_println!(progress, "[ok] combined $(length(loaded)) runs" *
                                  (skipped > 0 ? "  (skipped $skipped)" : ""))
    _progress_println!(progress, "     samples = $(meta[:n_samples])")
    _progress_println!(progress, "     span    = $(meta[:start_time]) -> $(meta[:end_time])")
    return combined, fmt
end

function _load_metronix_site(dir::AbstractString;
                             progress::Union{ProgressConsole, Nothing} = nothing)
    runs = metronix_site_runs(dir)
    isempty(runs) && error("No Metronix meas_ directories found in: $dir")
    rates = sort(collect(keys(runs)))
    length(rates) > 1 && @warn "Metronix site has multiple sampling rates; loading the first. " *
        "Split it first with scripts/split_metronix_by_rate.jl" rates
    rate = first(rates)
    meas_dirs = sort(runs[rate])
    site_name = _site_name_from_dir(dir)

    _progress_println!(progress, "> site      = $(site_name)")
    _progress_println!(progress, "> rate      = $(_format_fs(rate))")
    _progress_println!(progress, "> found $(length(meas_dirs)) run(s)")

    exact_fs = Float64(rate)
    tas = TimeArray[]
    for (i, d) in enumerate(meas_dirs)
        _progress_println!(progress, "  [$i/$(length(meas_dirs))] $(basename(d))")
        run = read_metronix(d)
        i == 1 && (exact_fs = sampling_rate(run))
        push!(tas, to_timearray(run; axis = :datetime))
    end
    ta = length(tas) == 1 ? tas[1] : _combine_site_timearrays(tas, site_name)
    filled = _fill_time_gaps(ta)
    md = _ta_meta(filled)
    if md isa AbstractDict
        md[:source_format] = :metronix
        md[:site_dir] = abspath(dir)
        md[:sample_rate] = exact_fs
        md[:metronix_rate] = exact_fs
    end
    _progress_println!(progress, "[ok] loaded $(length(meas_dirs)) run(s) @ $(_format_fs(rate))")
    return filled, :metronix
end

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

function _show_progress_window(title::AbstractString;
                                window_size = (900, 540),
                                max_lines::Int = 30)
    _try_set_transparent_framebuffer(true)
    fig = Figure(;
        size = window_size,
        backgroundcolor = TERM_BG,
        fontsize = 12,
        figure_padding = (16, 14, 12, 14),
    )
    Label(fig[1, 1], "▶  " * String(title);
          fontsize = 13, font = TERM_FONT, color = TERM_TITLE,
          halign = :left, tellwidth = false)
    Box(fig[2, 1]; color = RGBAf(TERM_DIM.r, TERM_DIM.g, TERM_DIM.b, 0.6),
        strokevisible = false)
    text_obs = Observable("")
    Label(fig[3, 1], text_obs;
          fontsize = 12, font = TERM_FONT, color = TERM_FG,
          halign = :left, valign = :bottom, justification = :left,
          tellwidth = false, tellheight = false, word_wrap = true)
    rowsize!(fig.layout, 1, Fixed(22))
    rowsize!(fig.layout, 2, Fixed(2))
    rowsize!(fig.layout, 3, Auto(true, 1.0))
    rowgap!(fig.layout, 6)
    screen = nothing
    try
        screen = display(GLMakie.Screen(; title = String(title),
                                          visible = true,
                                          focus_on_show = true),
                         fig)
        _center_and_float_window!(screen)
    catch err
        @warn "Could not open progress window; falling back to log only" exception = err
    end
    _try_set_transparent_framebuffer(false)
    return ProgressConsole(screen, text_obs, String[], max_lines,
                           ReentrantLock(), Threads.Atomic{Bool}(false))
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
    _progress_println!(progress, "")
    _progress_println!(progress, "> writing combined file:")
    _progress_println!(progress, "    $(out_path)")
    _write_data_file(out_path, ta, fmt)
    _progress_println!(progress, "[ok] wrote $(basename(out_path))")
    return out_path
end

function _fill_time_gaps(ta::TimeArray)
    times = _ensure_datetime(_ta_timestamps(ta))
    n = length(times)
    n <= 1 && return ta
    fs = _sample_rate_from_timearray(ta)
    fs > 0 || return ta
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
    n_cols = size(vals, 2)
    @inbounds for i in 1:n_rows
        bad = false
        for j in 1:n_cols
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

function _decimate_minmax!(
    xs::Vector{Float64},
    ys_clean::Vector{Float32},
    ys_masked::Vector{Float32},
    secs::Vector{Float64},
    col::AbstractVector{<:Real},
    masked::BitVector,
    idx_lo::Int,
    idx_hi::Int,
    n_buckets::Int,
)
    empty!(xs)
    empty!(ys_clean)
    empty!(ys_masked)
    n_window = idx_hi - idx_lo + 1
    n_window <= 0 && return 0

    if n_window <= 2 * n_buckets
        @inbounds for i in idx_lo:idx_hi
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
        clean_min = Inf
        clean_max = -Inf
        masked_min = Inf
        masked_max = -Inf
        for i in a:c
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
    end
    return length(xs)
end

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

    col1 = @view app.raw_values[:, 1]
    _decimate_minmax!(app.line_x[], app.line_clean[1][], app.line_masked[1][],
        secs, col1, masked, idx_lo, idx_hi, _PLOT_BUCKETS)
    notify(app.line_x)                                 # buffers refilled in place
    notify(app.line_clean[1])
    notify(app.line_masked[1])

    for j in 2:n_channels
        col_j = @view app.raw_values[:, j]
        _decimate_minmax!(app.line_x_scratch, app.line_clean[j][], app.line_masked[j][],
            secs, col_j, masked, idx_lo, idx_hi, _PLOT_BUCKETS)
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
    if app.selection_visible[]
        lo, hi = app.selection[]
        sel_text = "Selection " * _format_dt(app, lo) * "  →  " * _format_dt(app, hi)
    else
        sel_text = "Left-drag to select  ·  Right-click to mask the selection  ·  Right-drag = pan  ·  Scroll = zoom y"
    end
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

function _visible_good_index_segments(app::TKApp, x_lo::Float64, x_hi::Float64)
    secs = app.time_seconds
    masked = app.mask.masked
    isempty(secs) && return Tuple{Int, Int}[]
    idx_lo = searchsortedfirst(secs, x_lo)
    idx_hi = searchsortedlast(secs, x_hi)
    idx_lo = clamp(idx_lo, 1, length(secs))
    idx_hi = clamp(idx_hi, 1, length(secs))
    idx_lo > idx_hi && return Tuple{Int, Int}[]
    segs = Tuple{Int, Int}[]
    active = false
    start = idx_lo
    for i in idx_lo:idx_hi
        if !masked[i] && !active
            active = true
            start = i
        elseif masked[i] && active
            push!(segs, (start, i - 1))
            active = false
        end
    end
    active && push!(segs, (start, idx_hi))
    return segs
end

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

function _fmt_dur(t::Real)
    t < 60 && return @sprintf("%.1f s", t)
    t < 3600 && return @sprintf("%.1f min", t / 60)
    return @sprintf("%.1f h", t / 3600)
end

function _format_psd_header(nfft::Integer, fs::Real)
    df = fs / nfft
    return "nfft = $(nfft)   ·   df = $(_fmt_hz(df))   ·   f_Nyq = $(_fmt_hz(fs / 2))   ·   seg = $(_fmt_dur(nfft / fs))"
end

_spectra_info_idle() =
    "Time view  ·  use the View menu to add Spectra or Spectrogram panels for frequency content"

function _spectra_details_text(mode::Symbol, nfft::Integer, fs::Real; n_segments = nothing)
    metrics = _format_psd_header(nfft, fs)
    if mode === :time_spectrogram
        return "Spectrogram · STFT     |     x: time [s]  ·  y: frequency [Hz], log  ·  color: log10 power" *
               "     |     Hann window, 50% overlap, mean-detrended; masked windows left blank     |     " * metrics
    end
    segtxt = n_segments === nothing ? "unmasked segments" :
             "$(n_segments) unmasked segment" * (n_segments == 1 ? "" : "s")
    return "PSD · Welch's method     |     x: frequency [Hz], log  ·  y: PSD [amplitude^2/Hz], log" *
           "     |     Hann window, 50% overlap, mean-detrended; averaged over $(segtxt)     |     " * metrics
end

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
    if isfinite(ymin) && isfinite(ymax) && ymin < ymax
        ylims!(ax, ymin * 0.5, ymax * 2.0)
    end
end

function _compute_psd_for_window!(app::TKApp)
    app.view_mode[] === :time_spectra || return app
    isempty(app.psd_axes) && return app
    x_lo, x_hi = _visible_x_window(app)
    segs = _visible_good_index_segments(app, x_lo, x_hi)
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
    if app.psd_header !== nothing
        header_text = isempty(segs) ?
            "Window too short for nfft = $(nfft)" :
            _spectra_details_text(:time_spectra, nfft, fs; n_segments = n_used)
        app.psd_header.text[] = header_text
    end
    return app
end

function _compute_spectrogram_for_window!(app::TKApp)
    app.view_mode[] === :time_spectrogram || return app
    isempty(app.spec_axes) && return app
    x_lo, x_hi = _visible_x_window(app)
    secs = app.time_seconds
    isempty(secs) && return app
    idx_lo = clamp(searchsortedfirst(secs, x_lo), 1, length(secs))
    idx_hi = clamp(searchsortedlast(secs, x_hi), 1, length(secs))
    idx_lo > idx_hi && return app
    nfft, fs = _current_nfft(app)
    workspace = _spectral_workspace!(app, nfft, fs; noverlap = nfft ÷ 2)
    masked_window = view(app.mask.masked, idx_lo:idx_hi)
    for j in 1:length(app.spec_axes)
        x_view = view(app.raw_values, idx_lo:idx_hi, j)
        freqs, times, spec = _stft_psd(x_view, masked_window, fs;
            nfft = nfft, noverlap = nfft ÷ 2, workspace = workspace)
        if isempty(freqs)
            app.spec_times[j][] = Float64[0.0, 1.0]
            app.spec_freqs[j][] = Float64[1.0, 2.0]
            app.spec_matrix[j][] = fill(NaN, 2, 2)
            continue
        end
        log_spec = Matrix{Float64}(undef, length(times), length(freqs) - 1)
        @inbounds for ti in eachindex(times), fi in 2:length(freqs)
            v = spec[fi, ti]
            log_spec[ti, fi - 1] = (isfinite(v) && v > 0) ? log10(v) : NaN
        end
        app.spec_times[j][] = times .+ x_lo
        app.spec_freqs[j][] = freqs[2:end]
        app.spec_matrix[j][] = log_spec
        finite_vals = filter(isfinite, log_spec)
        if !isempty(finite_vals)
            lo, hi = extrema(finite_vals)
            if lo < hi
                app.spec_axes[j].limits[] = (nothing, nothing)
                xlims!(app.spec_axes[j], first(times) + x_lo, last(times) + x_lo)
                ylims!(app.spec_axes[j], freqs[2], freqs[end])
            end
        end
    end
    if app.psd_header !== nothing
        app.psd_header.text[] = _spectra_details_text(:time_spectrogram, nfft, fs)
    end
    return app
end

function _recompute_spectra!(app::TKApp)
    mode = app.view_mode[]
    if mode === :time_spectra
        _compute_psd_for_window!(app)
    elseif mode === :time_spectrogram
        _compute_spectrogram_for_window!(app)
    else
        app.psd_header !== nothing && (app.psd_header.text[] = _spectra_info_idle())
        return app
    end
    _refresh_status!(app)
    return app
end

function _schedule_spectra_recompute!(app::TKApp; delay_seconds::Real = 0.25)
    app.view_mode[] === :time && return app
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

function _clear_psd_axes!(app::TKApp)
    for ax in app.psd_axes
        delete!(ax)
    end
    empty!(app.psd_axes)
    empty!(app.psd_freqs)
    empty!(app.psd_values)
    for ax in app.spec_axes
        delete!(ax)
    end
    empty!(app.spec_axes)
    empty!(app.spec_times)
    empty!(app.spec_freqs)
    empty!(app.spec_matrix)
    return app
end

function _build_axes!(app::TKApp, ta::TimeArray, existing::Vector{Axis})
    for ax in existing
        delete!(ax)
    end
    empty!(existing)
    empty!(app.line_clean)
    empty!(app.line_masked)
    _clear_psd_axes!(app)

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
    app.line_x[] = Float64[]

    spectra_on = app.view_mode[] === :time_spectra
    spectrogram_on = app.view_mode[] === :time_spectrogram
    right_col_on = spectra_on || spectrogram_on
    n = length(names)
    axes = Axis[]
    psd_axes = Axis[]
    spec_axes = Axis[]
    for (i, name) in enumerate(names)
        is_last = i == n
        unit_str = get(units_map, _symbolize(name), component_units(_symbolize(name)))
        ax = Axis(app.plot_layout[i, 1];
            ylabel = "$(_display_label(name)) [$(unit_str)]",
            ylabelrotation = pi / 2,
            ylabelpadding = 8.0,
            ylabelsize = 12,
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

        col = _component_color(name)
        if length(secs) == 1
            scatter!(ax, app.line_x, clean_obs; color = col, markersize = 4)
            scatter!(ax, app.line_x, masked_obs; color = TK_MUTED, markersize = 4)
        else
            lines!(ax, app.line_x, clean_obs; color = col, linewidth = 1.4, joinstyle = :round)
            lines!(ax, app.line_x, masked_obs; color = TK_MUTED, linewidth = 1.4, joinstyle = :round)
        end

        deregister_interaction!(ax, :rectanglezoom)
        register_interaction!(ax, :tk_select, DragSelect(app, false, 0.0))

        push!(axes, ax)

        if spectra_on
            ax_psd = Axis(app.plot_layout[i, 2];
                xscale = log10,
                yscale = log10,
                yticklabelsize = 10,
                xticklabelsize = 10,
                xticklabelsvisible = is_last,
                xticksvisible = is_last,
                xgridvisible = false,
                ygridvisible = false,
                topspinevisible = false,
                rightspinevisible = false,
                yticks = LogTicks(WilkinsonTicks(3)),
                xticks = LogTicks(WilkinsonTicks(4)),
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
            deregister_interaction!(ax_psd, :rectanglezoom)
            push!(psd_axes, ax_psd)
        end

        if spectrogram_on
            ax_spec = Axis(app.plot_layout[i, 2];
                yscale = log10,
                yticklabelsize = 10,
                xticklabelsize = 10,
                xticklabelsvisible = is_last,
                xticksvisible = is_last,
                xgridvisible = false,
                ygridvisible = false,
                topspinevisible = false,
                rightspinevisible = false,
                yticks = LogTicks(WilkinsonTicks(3)),
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
                hidexdecorations!(ax_spec; ticks = true, ticklabels = true, grid = false)
            end
            t_obs = Observable{Vector{Float64}}(Float64[0.0, 1.0])
            f_obs = Observable{Vector{Float64}}(Float64[1.0, 2.0])
            m_obs = Observable{Matrix{Float64}}(fill(NaN, 2, 2))
            push!(app.spec_times, t_obs)
            push!(app.spec_freqs, f_obs)
            push!(app.spec_matrix, m_obs)
            heatmap!(ax_spec, t_obs, f_obs, m_obs;
                colormap = TK_SPECTRAL_COLORMAP, nan_color = RGBAf(0, 0, 0, 0))
            deregister_interaction!(ax_spec, :rectanglezoom)
            push!(spec_axes, ax_spec)
        end
    end
    if length(axes) > 1
        linkxaxes!(axes...)
        rowgap!(app.plot_layout, 8)
    end
    if spectra_on && length(psd_axes) > 1
        linkxaxes!(psd_axes...)
    end
    if spectrogram_on && length(spec_axes) > 1
        linkxaxes!(spec_axes...)
        linkyaxes!(spec_axes...)
    end
    if right_col_on
        colsize!(app.plot_layout, 1, Auto(true, 0.6))
        colsize!(app.plot_layout, 2, Auto(true, 0.4))
        colgap!(app.plot_layout, 12)
    else
        try
            trim!(app.plot_layout)
        catch
        end
        colsize!(app.plot_layout, 1, Auto(true, 1.0))
    end
    append!(existing, axes)
    append!(app.psd_axes, psd_axes)
    append!(app.spec_axes, spec_axes)
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

function _apply_loaded_data!(app::TKApp, ta::TimeArray, fmt::Symbol, source_path::AbstractString)
    app.data = ta
    app.mask = TimekeeperMask(ta)
    app.source_format = fmt
    app.source_path = String(source_path)
    app.selection_visible[] = false
    app.window_start[] = 0.0
    _build_axes!(app, ta, app.axes)
    _auto_mask_nan!(app.mask, app.raw_values)
    _refresh_slider_range!(app)
    app.summary_label.text[] = _summary_text(ta)
    _refresh_mask_overlay!(app)
    _update_x_window!(app)
    return app
end

function _load_site_any(dir::AbstractString; progress::Union{ProgressConsole, Nothing} = nothing)
    return is_metronix_site(dir) ? _load_metronix_site(dir; progress = progress) :
           _load_site_directory(dir; progress = progress)
end

function _load_into_app!(app::TKApp, path::AbstractString)
    ta, fmt = isdir(path) ? _load_site_any(path) : _load_data_file(path)
    return _apply_loaded_data!(app, ta, fmt, path)
end

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
    fig = Figure(; size = size, backgroundcolor = :white, fontsize = 12,
        figure_padding = (14, 14, 8, 8))

    toolbar = GridLayout(fig[1, 1]; tellheight = true)
    summary_label = Label(toolbar[1, 1], _summary_text(ta);
        fontsize = 12, color = TK_GREY, halign = :left, tellwidth = false)

    actions = GridLayout(toolbar[1, 2]; tellheight = false, halign = :right)
    load_btn = _logo_button(actions[1, 1], "Load Run…", TK_LOGO_SKY; textcolor = TK_BLACK)
    load_site_btn = _logo_button(actions[1, 2], "Load Site…", TK_LOGO_SKY; textcolor = TK_BLACK)
    mask_btn = _logo_button(actions[1, 3], "Mask", TK_LOGO_CORAL)
    unmask_btn = _logo_button(actions[1, 4], "Unmask", TK_LOGO_SAGE)
    clear_btn = _logo_button(actions[1, 5], "Clear", TK_LOGO_SLATE)
    write_btn = _logo_button(actions[1, 6], "Write", TK_LOGO_TEAL)
    view_label = Label(actions[1, 7], "View:"; fontsize = 11, color = TK_GREY)
    view_menu = _logo_menu(actions[1, 8]; options = VIEW_OPTIONS, default = "Time", width = 170)
    window_label = Label(actions[1, 9], "Window:"; fontsize = 11, color = TK_GREY)
    window_menu = _logo_menu(actions[1, 10]; options = WINDOW_OPTIONS, default = "1 hour", width = 110)
    colgap!(actions, 6)

    colsize!(toolbar, 1, Auto(true, 1.0))
    colsize!(toolbar, 2, Auto(false))

    plot_layout = GridLayout(fig[2, 1]; tellheight = false)

    spectra_info = Label(fig[3, 1], _spectra_info_idle();
        fontsize = 11, color = TK_BLACK, halign = :left, tellwidth = false)

    slider_grid = GridLayout(fig[4, 1]; tellheight = true)
    Label(slider_grid[1, 1], "Scroll"; fontsize = 11, color = TK_GREY, halign = :right)
    prev_btn = _logo_button(slider_grid[1, 2], "<", TK_CTRL_BTN; textcolor = TK_BLACK, fontsize = 16, font = :bold, height = 22)
    slider = Slider(slider_grid[1, 3]; range = 0.0:1.0:0.0, startvalue = 0.0,
        linewidth = 11.0,
        color_inactive = TK_CTRL_TRACK,
        color_active_dimmed = RGBAf(TK_CTRL_ACCENT.r, TK_CTRL_ACCENT.g, TK_CTRL_ACCENT.b, 0.45),
        color_active = TK_CTRL_ACCENT)
    next_btn = _logo_button(slider_grid[1, 4], ">", TK_CTRL_BTN; textcolor = TK_BLACK, fontsize = 16, font = :bold, height = 22)
    colsize!(slider_grid, 1, Fixed(50))
    colsize!(slider_grid, 2, Fixed(44))
    colsize!(slider_grid, 3, Auto(true, 1.0))
    colsize!(slider_grid, 4, Fixed(44))
    colgap!(slider_grid, 10)

    status_label = Label(fig[5, 1], "";
        fontsize = 10.5, color = TK_GREY, halign = :left, tellwidth = false)

    rowsize!(fig.layout, 1, Fixed(38))
    rowsize!(fig.layout, 2, Auto(true, 1.0))
    rowsize!(fig.layout, 3, Fixed(20))
    rowsize!(fig.layout, 4, Fixed(34))
    rowsize!(fig.layout, 5, Fixed(20))
    rowgap!(fig.layout, 6)

    selection = Observable((0.0, 0.0))
    selection_visible = Observable(false)
    mask_lows = Observable(Float64[])
    mask_highs = Observable(Float64[])
    window_seconds_obs = Observable(3600.0)
    window_start_obs = Observable(0.0)
    line_x_obs = Observable(Float64[])
    view_mode_obs = Observable(:time)

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
        spectra_info,
        Axis[],
        Observable{Vector{Float64}}[],
        Observable{Vector{Float64}}[],
        Observable{Matrix{Float64}}[],
        Dict{Tuple{Int, Float64, Int, Symbol, Symbol}, TKSpectralWorkspace}(),
        Ref{Union{Nothing, Timer}}(nothing),
    )

    _build_axes!(app, ta, app.axes)
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
    on(window_menu.selection) do secs
        secs === nothing && return
        app.window_seconds[] = Float64(secs)
        _refresh_slider_range!(app)
        _update_x_window!(app)
        _recompute_spectra!(app)
    end
    on(view_menu.selection) do mode
        mode === nothing && return
        mode === app.view_mode[] && return
        app.view_mode[] = Symbol(mode)
        _build_axes!(app, app.data, app.axes)
        _refresh_mask_overlay!(app)
        _update_x_window!(app)
        _recompute_spectra!(app)
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
        app.status_label.text[] = "Loading $(basename(path))…"
        @async begin
            try
                ta, fmt = fetch(Threads.@spawn(_load_data_file(path)))
                _apply_loaded_data!(app, ta, fmt, path)
                app.status_label.text[] = _ready_status_text(app)
            catch err
                err isa TaskFailedException && (err = err.task.exception)
                @warn "Could not load $path" exception = err
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
        console = _show_progress_window("TIMEKEEPERS // LOAD SITE  ::  $(site_name)")
        app.status_label.text[] = "Loading site $(site_name)…"
        @async begin
            try
                ta, fmt = _run_with_progress_pump(console) do
                    _load_site_any(dir; progress = console)
                end
                _progress_println!(console, "")
                _progress_println!(console, "> applying to viewer...")
                _flush_progress!(console)
                _apply_loaded_data!(app, ta, fmt, dir)
                app.status_label.text[] = "Loaded site $(site_name)"
                _progress_println!(console, "")
                _progress_println!(console, "[done] closing in 3 s ...")
                _flush_progress!(console)
                _close_progress_window!(console, 3.0)
            catch err
                @warn "Site load failed" dir exception = err
                msg = sprint(showerror, err)
                try
                    _progress_println!(console, "")
                    _progress_println!(console, "[error]")
                    for line in split(msg, '\n')
                        _progress_println!(console, "  " * line)
                    end
                    _flush_progress!(console)
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
            intervals = copy(app.mask.intervals)
            @async begin
                console = _show_progress_window("TIMEKEEPERS // WRITE METRONIX  ::  $(_site_name_from_dir(site_dir))")
                try
                    _progress_println!(console, "> source = $(site_dir)")
                    _progress_println!(console, "> cuts   = $(length(intervals)) interval(s)")
                    _progress_println!(console, "> writing split meas_ dirs...")
                    dest = _run_with_progress_pump(console) do
                        write_metronix_site_masked(site_dir; intervals = intervals)
                    end
                    _progress_println!(console, "[ok] wrote $(basename(dest))")
                    _flush_progress!(console)
                    app.status_label.text[] = "Wrote $(basename(dest))"
                    _close_progress_window!(console, 3.0)
                catch err
                    @warn "Metronix write failed" exception = err
                    app.status_label.text[] = "Write failed: $(sprint(showerror, err))"
                    try
                        _progress_println!(console, "[error] " * sprint(showerror, err))
                        _flush_progress!(console)
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
    ta, fmt = isdir(path) ? _load_site_any(path) : _load_data_file(path)
    return TKApp(ta; source_format = fmt, source_path = path, kwargs...)
end

function TKApp(; kwargs...)
    return TKApp(_interactive_timearray(); kwargs...)
end

function _open_app_screen(app::TKApp; maximize::Bool)
    try
        screen = display(GLMakie.Screen(; title = "Timekeepers", visible = false,
                                          focus_on_show = true), app.figure)
        _apply_timekeepers_icon!(screen)
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
        _apply_timekeepers_icon!(screen)
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
source format. The **View** menu adds per-channel PSD and spectrogram panels.

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
