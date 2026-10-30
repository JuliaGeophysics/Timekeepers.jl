# Processing.jl - robust estimation of magnetotelluric transfer functions.
# Author: @pankajkmishra
#
# The impedance Z and the tipper T relate the fields at one site in the
# frequency domain:
#   [Ex, Ey] = Z [Hx, Hy]        Hz = T [Hx, Hy]
# The file estimates them from the time series in three set-ups:
# - single site: the site gives the outputs (E, Hz) and the inputs (Hx, Hy)
# - base site: a site near the target gives the inputs. Thus, a site with only
#   electric channels (a telluric site) gets an impedance
# - remote reference: one or more far sites give a reference field. Noise in
#   the local inputs is not in the reference. Thus, the estimate is free of
#   the downward bias of a single site. A base site and remote sites can go
#   together
#
# The steps are:
# 1. Load the runs of each site at one rate, change them to physical units
#    (Calibration.jl), and keep the spans that all the sites recorded
#    together, without gaps and masks.
# 2. A cascade of decimation levels: each level is a factor of 4 below the one
#    above, with an anti-alias FIR filter. Each level gives the frequencies
#    from 1/16 to 1/4 of its rate. The last level also gives its lowest
#    harmonics.
# 3. At each level, an AR filter whitens each channel. Then a Hann window, an
#    FFT, the correction for the AR filter and the sensor, and a rotation to
#    geographic north and east give the Fourier coefficients.
# 4. Frequency bands on a fixed grid of periods (the same for all rates and
#    sites). Each band collects its harmonics in all the windows.
# 5. For each band, a robust regression. Two methods:
#    - CT2004 (the default), after Chave & Thomson (2004): two stages with
#      bounded influence. Stage 1, with a reference: a robust fit of the local
#      inputs to the reference. Its weights remove noise bursts in the local
#      or the reference field, and the fitted inputs are free of local noise.
#      Stage 2: a robust instrumental variable fit of each output, a Huber
#      M-estimate, then a redescending weight, with a bounded influence
#      weight for points of high leverage.
#    - EB1986, after Egbert & Booker (1986): one robust M-estimate with Huber
#      weights on the residuals of each output. The reference enters as it
#      is: the instrument itself with one remote site, a plain least squares
#      projection with more. No stage 1 weights, no redescending weights and
#      no leverage weights.
#    Neither wins on all data. CT2004 resists noise in the local H and
#    windows of high leverage. EB1986 keeps more of the data when the remote
#    site is far or weakly coherent, and the strong natural events that the
#    redescending and leverage weights can remove.
# 6. The errors: a jackknife that deletes blocks of windows, one block at a
#    time, with the final weights. Adjacent harmonics of one window are not
#    independent. The blocks keep them together.
#
# At more than one rate, the function estimates each rate. Then, for each
# period, it keeps the estimate with the smallest relative error

const TF_DECIMATION = 4
const TF_FIR_TAPS = 65
const TF_FIR_BETA = 8.0
const TF_BAND_HARMONICS = 3
const TF_EDGE_COHERENCE = 0.5
const TF_EDGE_MARGIN = 0.95

"""
    TransferFunction

The impedance and the tipper of one site, with their errors, as
[`estimate_tf`](@ref) returns them.

The fields are in geographic coordinates: `x` is north and `y` is east. `Z`
is in (mV/km)/nT and has the time dependence `exp(+iωt)`. Thus, the phase of
`Zxy` over a uniform half-space is +45°.

# Fields
- `site::String` -- the site.
- `mode::Symbol` -- `:single`, `:base`, `:remote` or `:base_remote`.
- `base::String` -- the site that gave the inputs Hx, Hy (`""` for the site
  itself).
- `remote::Vector{String}` -- the reference sites.
- `latitude`, `longitude`, `elevation` -- the position of the site.
- `periods::Vector{Float64}` -- seconds, increasing.
- `Z::Array{ComplexF64, 3}` -- `2 × 2 × nperiods`. `NaN` where there is no
  estimate.
- `Z_var::Array{Float64, 3}` -- the variance `E|δZ|²` of each element.
- `T::Matrix{ComplexF64}` -- `2 × nperiods`: `Tzx`, `Tzy`.
- `T_var::Matrix{Float64}` -- the variance of each tipper element.
- `coherence::Matrix{Float64}` -- `3 × nperiods`: the multiple coherence of
  Ex, Ey and Hz with the inputs, with the robust weights.
- `ref_coherence::Matrix{Float64}` -- `2 × nperiods`: the coherence of the
  inputs Hx, Hy with the reference (`NaN` without a reference).
- `n_windows::Vector{Int}` -- the windows in each band.
- `n_obs::Vector{Int}` -- the Fourier coefficients in each band.
- `sample_rate::Vector{Float64}` -- the rate that gave each period.
- `metadata::Dict{Symbol, Any}` -- the options, the sensors and the data
  spans.

See also [`apparent_resistivity`](@ref), [`impedance_phase`](@ref),
[`write_modem`](@ref).
"""
struct TransferFunction
    site::String
    mode::Symbol
    base::String
    remote::Vector{String}
    latitude::Float64
    longitude::Float64
    elevation::Float64
    periods::Vector{Float64}
    Z::Array{ComplexF64, 3}
    Z_var::Array{Float64, 3}
    T::Matrix{ComplexF64}
    T_var::Matrix{Float64}
    coherence::Matrix{Float64}
    ref_coherence::Matrix{Float64}
    n_windows::Vector{Int}
    n_obs::Vector{Int}
    sample_rate::Vector{Float64}
    metadata::Dict{Symbol, Any}
end

function Base.show(io::IO, tf::TransferFunction)
    span = isempty(tf.periods) ? "no periods" :
           @sprintf("%d periods %.3g–%.3g s", length(tf.periods), first(tf.periods), last(tf.periods))
    refs = tf.mode === :single ? "" :
           tf.mode === :base ? ", base $(tf.base)" :
           tf.mode === :remote ? ", remote $(join(tf.remote, "+"))" :
           ", base $(tf.base), remote $(join(tf.remote, "+"))"
    print(io, "TransferFunction(\"$(tf.site)\", :$(tf.mode)$refs, $span)")
end

const _MU0 = 4π * 1.0e-7

"""
    apparent_resistivity(tf::TransferFunction) -> (rho, rho_err)

The apparent resistivity `ρ = 0.2 T |Z|²` of each element in Ω·m, and its
error, both `2 × 2 × nperiods`. The error comes from the standard deviation
of `Z` along one axis, `σ = √(var/2)`: `δρ = 2ρσ/|Z|`.
"""
function apparent_resistivity(tf::TransferFunction)
    T = reshape(tf.periods, 1, 1, :)
    rho = 0.2 .* T .* abs2.(tf.Z)
    σ = sqrt.(tf.Z_var ./ 2)
    return rho, 2 .* rho .* σ ./ abs.(tf.Z)
end

"""
    impedance_phase(tf::TransferFunction) -> (phi, phi_err)

The phase of each element of `Z` in degrees, in `(-180°, 180°]`, and its error
`δφ = σ/|Z|` in degrees, both `2 × 2 × nperiods`. Refer to
[`apparent_resistivity`](@ref) for `σ`.
"""
function impedance_phase(tf::TransferFunction)
    σ = sqrt.(tf.Z_var ./ 2)
    return rad2deg.(angle.(tf.Z)), rad2deg.(min.(σ ./ abs.(tf.Z), π))
end

"""
    rotate_tf(tf::TransferFunction, angle) -> TransferFunction

The transfer function in axes that are turned clockwise by `angle` degrees
(east of north): `Z' = R Z Rᵀ` and `T' = T Rᵀ`. The variances turn with the
same rule, for independent elements.
"""
function rotate_tf(tf::TransferFunction, angle::Real)
    c, s = cosd(angle), sind(angle)
    R = [c s; -s c]
    Z = similar(tf.Z)
    V = similar(tf.Z_var)
    T = similar(tf.T)
    TV = similar(tf.T_var)
    R2 = R .^ 2
    for k in axes(tf.Z, 3)
        Z[:, :, k] = R * tf.Z[:, :, k] * R'
        V[:, :, k] = R2 * tf.Z_var[:, :, k] * R2'
        T[:, k] = R * tf.T[:, k]
        TV[:, k] = R2 * tf.T_var[:, k]
    end
    md = copy(tf.metadata)
    md[:rotation] = get(tf.metadata, :rotation, 0.0) + angle
    return TransferFunction(tf.site, tf.mode, tf.base, tf.remote, tf.latitude, tf.longitude,
                            tf.elevation, tf.periods, Z, V, T, TV, tf.coherence, tf.ref_coherence,
                            tf.n_windows, tf.n_obs, tf.sample_rate, md)
end

#---------- sites and their data -----

# One site at one rate: its contiguous good segments and its channel
# responses. `k0` is the index of the first sample on a grid of `fs` samples
# from 1970
const _Col = SubArray{Float64, 1, Vector{Float64}, Tuple{UnitRange{Int}}, true}

_subview(v::_Col, r::UnitRange{Int}) = view(parent(v), parentindices(v)[1][r])
_whole(v::Vector{Float64}) = view(v, 1:length(v))

struct _Segment
    k0::Int
    cols::Dict{Symbol, _Col}
end

struct _SiteData
    name::String
    latitude::Float64
    longitude::Float64
    elevation::Float64
    segments::Vector{_Segment}
    responses::Dict{Symbol, _ChannelResponse}
end

_site_label(x::SurveySite) = x.name
_site_label(x::TimekeeperRun) = x.site
_site_label(x::AbstractVector{TimekeeperRun}) = isempty(x) ? "" : first(x).site
_site_label(x::AbstractString) = (s = _as_site(x); s isa SurveySite ? s.name : _site_from_path(x))

# A path is a site directory (scanned for its runs) or one data file
function _as_site(x::AbstractString)
    p = _norm_path(x)
    if isdir(p)
        s = _scan_site(p)
        s === nothing && error("No recordings in $p")
        return s
    end
    isfile(p) || error("No such file or directory: $p")
    return read_timekeeper(p)
end
_as_site(x) = x

_rates_of(x::SurveySite) = site_rates(x)
_rates_of(x::TimekeeperRun) = [_rate_key(sampling_rate(x))]
_rates_of(x::AbstractVector{TimekeeperRun}) = sort!(unique([_rate_key(sampling_rate(r)) for r in x]))

_components_of(x::SurveySite) = site_components(x)
_components_of(x::TimekeeperRun) = components(x)
_components_of(x::AbstractVector{TimekeeperRun}) = unique(c for r in x for c in components(r))

# A rate `r` gives the rate `fs`: the same rate, or an integer multiple that
# a decimation brings down to it (a 10 Hz GEOMAG with 1 Hz LEMI-424 sites)
_decimation_factor(r::Real, fs::Real) = _same_rate(r, fs) ? 1 :
    (r > fs && isapprox(r / fs, round(r / fs); rtol = 1.0e-6)) ? round(Int, r / fs) : 0
_gives_rate(x, fs::Real) = any(r -> _decimation_factor(r, fs) > 0, _rates_of(x))

# The runs of a site at `fs`. A site without runs at `fs` gives its runs at an
# integer multiple of `fs`, decimated
function _runs_at(x::SurveySite, fs, comps)
    runs = TimekeeperRun[]
    native = any(r -> _same_rate(r.sample_rate, fs), x.runs)
    for r in x.runs
        k = _decimation_factor(r.sample_rate, fs)
        (k == 1 || (k > 1 && !native)) || continue
        want = [c for c in comps if c in r.components]
        isempty(want) && continue
        # a run that cannot be read (e.g. a test file whose clock jumped back
        # at the GPS fix) is left out, not the whole site
        run = try
            x.format === :metronix ? read_metronix(r.path; site = x.name, components = want) :
            read_timekeeper(r.path; format = x.format, site = x.name)
        catch err
            @warn "$(x.name): $(basename(r.path)) is left out: $(sprint(showerror, err))" maxlog = 8
            continue
        end
        push!(runs, k == 1 ? run : _decimate_run(run, fs, comps))
    end
    return runs
end
function _runs_at(x::TimekeeperRun, fs, comps)
    k = _decimation_factor(sampling_rate(x), fs)
    return k == 0 ? TimekeeperRun[] : k == 1 ? [x] : [_decimate_run(x, fs, comps)]
end
function _runs_at(x::AbstractVector{TimekeeperRun}, fs, comps)
    native = any(r -> _same_rate(sampling_rate(r), fs), x)
    out = TimekeeperRun[]
    for r in x
        k = _decimation_factor(sampling_rate(r), fs)
        k == 1 ? push!(out, r) : (k > 1 && !native) ? push!(out, _decimate_run(r, fs, comps)) : nothing
    end
    return out
end

# `run` at the lower rate `fs`: a linear-phase low-pass FIR (windowed sinc,
# Blackman, cut-off at 0.8 of the new Nyquist frequency), then the samples on
# the grid of `fs` from 1970. The filter has no delay. A NaN reaches the
# samples within its half length
function _decimate_run(run::TimekeeperRun, fs::Real, comps)
    chans = Dict{Symbol, TimekeeperChannel}()
    for c in comps
        haskey(run.channels, c) || continue
        ch = run.channels[c]
        r = ch.sample_rate
        k = _decimation_factor(r, fs)
        k >= 1 || error("$(run.site): $(c) at $(_fs_label(r)) does not give $(_fs_label(fs))")
        h = 8k                                   # half length of the filter in samples
        fc = 0.8 * (fs / 2) / r                  # cycles per input sample
        w = [(n == 0 ? 2fc : sin(2π * fc * n) / (π * n)) *
             (0.42 + 0.5 * cos(π * n / (h + 1)) + 0.08 * cos(2π * n / (h + 1))) for n in -h:h]
        w ./= sum(w)
        # the first input sample on the output grid
        period_ms = 1000 / fs
        t0 = Dates.value(ch.start - _EPOCH)
        wait_ms = mod(-t0, period_ms)
        i0 = round(Int, wait_ms * r / 1000) + 1
        abs(wait_ms - (i0 - 1) * 1000 / r) > 1 &&
            @warn "$(run.site): the samples of $c are not on the grid of $(_fs_label(fs)); the nearest is used" maxlog = 2
        d = ch.data
        n = length(d)
        idx = i0:k:n
        out = Vector{Float64}(undef, length(idx))
        @inbounds for (j, i) in enumerate(idx)
            s = 0.0
            if i - h < 1 || i + h > n
                s = NaN
            else
                for m in -h:h
                    s += w[m + h + 1] * d[i + m]
                end
            end
            out[j] = s
        end
        start = ch.start + Millisecond(round(Int, (i0 - 1) * 1000 / r))
        header = copy(ch.header)
        header["decimated_from"] = r
        chans[c] = TimekeeperChannel(c, out, Float64(fs), start, ch.units, ch.source_file, header)
    end
    md = copy(run.metadata)
    md[:sample_rate] = Float64(fs)
    md[:decimated_from] = sampling_rate(run)
    return TimekeeperRun(run.site, run.instrument, run.source_format, chans, md)
end

function _position_of(x::SurveySite, runs)
    return (x.latitude, x.longitude, x.elevation)
end

_position_of(x, runs) = _position_of_runs(runs)

function _runs_in_span(runs, span)
    span === nothing && return runs
    function overlaps(r)
        isempty(r.channels) && return false
        ch = first(values(r.channels))
        t1 = ch.start + Millisecond(round(Int, 1000 * length(ch.data) / ch.sample_rate))
        return ch.start < span[2] && t1 > span[1]
    end
    return filter(overlaps, runs)
end

function _position_of_runs(runs)
    for r in runs
        lat = get(r.metadata, :latitude, NaN)
        lon = get(r.metadata, :longitude, NaN)
        lat isa Real && lon isa Real && isfinite(lat) && isfinite(lon) &&
            return (Float64(lat), Float64(lon), Float64(get(r.metadata, :elevation, 0.0)))
        for ch in values(r.channels)
            haskey(ch.header, "ats_header_bytes") || continue
            p = _ats_position(ch.header["ats_header_bytes"])
            isfinite(p[1]) && return p
        end
    end
    return (NaN, NaN, NaN)
end

const _EPOCH = DateTime(1970, 1, 1)

function _grid_index(t::DateTime, fs::Real)
    k = Dates.value(t - _EPOCH) / 1000 * fs
    kr = round(Int, k)
    abs(k - kr) > 0.01 &&
        @warn "A run starts $(round(abs(k - kr); digits = 3)) samples off the common sample grid; it is moved to the nearest sample" maxlog = 2
    return kr
end

# The masked intervals of a site as index ranges on its grid
function _mask_ranges(masks, name, fs)
    iv = masks === nothing ? () : get(masks, name, ())
    return [(_grid_index(a, fs), _grid_index(b, fs)) for (a, b) in iv]
end

function _load_site(x, fs, comps; calibration, dipole, azimuths, masks, span)
    runs = _runs_in_span(_runs_at(x, fs, comps), span)
    name = _site_label(x)
    isempty(runs) && error("$name has no runs at $(_fs_label(fs))" * (span === nothing ? "" : " in the span"))
    lat, lon, elev = span === nothing ? _position_of(x, runs) : _position_of_runs(runs)
    isfinite(lat) || ((lat, lon, elev) = _position_of(x, runs))
    responses = Dict{Symbol, _ChannelResponse}()
    segments = _Segment[]
    cut = _mask_ranges(masks, name, fs)
    span_ranges = span === nothing ? nothing : (_grid_index(span[1], fs), _grid_index(span[2], fs))
    # the responses come from the longest run: an automatic fluxgate direction
    # (d/orientations.dat) needs a long static record, not a short test file
    full = [r for r in runs if all(c -> haskey(r.channels, c), comps)]
    if !isempty(full)
        longest = full[argmax([minimum(length(r.channels[c].data) for c in comps) for r in full])]
        orientation = any(in(ORIENTATION_COLUMNS), comps) ?
                      site_orientation(longest; site_dir = String(get(longest.metadata, :site_dir,
                                                                      dirname(first(values(longest.channels)).source_file)))) : nothing
        for c in comps
            responses[c] = _channel_response(longest, c; calibration, dipole, azimuths, orientation)
        end
    end
    for run in runs
        all(c -> haskey(run.channels, c), comps) || continue
        k0 = _grid_index(minimum(run.channels[c].start for c in comps), fs)
        n = minimum(length(run.channels[c].data) for c in comps)
        bad = falses(n)
        for c in comps
            d = run.channels[c].data
            off = _grid_index(run.channels[c].start, fs) - k0
            off == 0 || error("The channels of a run of $name start at different times")
            @inbounds for i in 1:n
                isfinite(d[i]) || (bad[i] = true)
            end
        end
        for (a, b) in cut
            lo, hi = max(a - k0 + 1, 1), min(b - k0 + 1, n)
            lo <= hi && (bad[lo:hi] .= true)
        end
        if span_ranges !== nothing
            lo, hi = span_ranges[1] - k0 + 1, span_ranges[2] - k0 + 1
            lo > 1 && (bad[1:min(lo - 1, n)] .= true)
            hi < n && (bad[max(hi + 1, 1):n] .= true)
        end
        i = 1
        while i <= n
            if bad[i]
                i += 1
                continue
            end
            j = i
            while j < n && !bad[j + 1]
                j += 1
            end
            cols = Dict{Symbol, _Col}(c => view(run.channels[c].data, i:j) for c in comps)
            push!(segments, _Segment(k0 + i - 1, cols))
            i = j + 1
        end
    end
    isempty(segments) && error("$name has no data with all of $(join(comps, ", ")) at $(_fs_label(fs))")
    sort!(segments; by = s -> s.k0)
    # segments that overlap (the same data in two runs) keep the first copy
    kept = _Segment[segments[1]]
    for s in @view segments[2:end]
        last_end = kept[end].k0 + length(first(values(kept[end].cols))) - 1
        if s.k0 <= last_end
            skip = last_end - s.k0 + 1
            len = length(first(values(s.cols)))
            skip >= len && continue
            s = _Segment(s.k0 + skip, Dict(c => _subview(v, (skip + 1):len) for (c, v) in s.cols))
        end
        push!(kept, s)
    end
    return _SiteData(name, lat, lon, elev, kept, responses)
end

_seg_range(s::_Segment) = (s.k0, s.k0 + length(first(values(s.cols))) - 1)

# The spans that all the sites recorded, as (lo, hi, segment index of each site)
function _common_spans(sites::Vector{_SiteData})
    acc = [(a, b, [i]) for (i, (a, b)) in enumerate(_seg_range.(sites[1].segments))]
    for s in sites[2:end]
        ranges = _seg_range.(s.segments)
        out = Tuple{Int, Int, Vector{Int}}[]
        i = j = 1
        while i <= length(acc) && j <= length(ranges)
            lo = max(acc[i][1], ranges[j][1])
            hi = min(acc[i][2], ranges[j][2])
            hi >= lo && push!(out, (lo, hi, vcat(acc[i][3], j)))
            acc[i][2] < ranges[j][2] ? (i += 1) : (j += 1)
        end
        acc = out
    end
    return acc
end

#---------- the channel plan -----

# The channels of one estimate in the column order of the Fourier matrices:
# outputs (Ex, Ey, Hz of the target), inputs (Hx, Hy), then reference pairs
struct _Plan
    sites::Vector{_SiteData}
    chans::Vector{Tuple{Int, Symbol}}        # (site index, component)
    outputs::Vector{Pair{Symbol, Int}}       # :ex, :ey, :hz => column
    inputs::Vector{Int}
    refs::Vector{Int}
    pairs::Vector{Tuple{Int, Int, Matrix{Float64}}}   # columns turned to north and east
end

function _rotation_pair(r1::_ChannelResponse, r2::_ChannelResponse)
    a1, a2 = r1.azimuth, r2.azimuth
    (isapprox(a1, 0; atol = 1.0e-6) && isapprox(a2, 90; atol = 1.0e-6)) && return nothing
    M = [cosd(a1) sind(a1); cosd(a2) sind(a2)]
    abs(det(M)) < 0.1 && error("Sensors at $(a1)° and $(a2)° are too close to parallel to resolve north and east")
    return inv(M)
end

function _make_plan(sites, target::Int, insite::Int, refsites::Vector{Int}, outputs::Vector{Symbol})
    chans = Tuple{Int, Symbol}[]
    outs = Pair{Symbol, Int}[]
    pairs = Tuple{Int, Int, Matrix{Float64}}[]
    col(s, c) = (push!(chans, (s, c)); length(chans))
    names = Dict(:e1 => :ex, :e2 => :ey, :bz => :hz)
    for c in outputs
        push!(outs, names[c] => col(target, c))
    end
    if :e1 in outputs && :e2 in outputs
        R = _rotation_pair(sites[target].responses[:e1], sites[target].responses[:e2])
        R === nothing || push!(pairs, (outs[1][2], outs[2][2], R))
    end
    function hpair(s)
        i, j = col(s, :bx), col(s, :by)
        R = _rotation_pair(sites[s].responses[:bx], sites[s].responses[:by])
        R === nothing || push!(pairs, (i, j, R))
        return [i, j]
    end
    inputs = hpair(insite)
    refs = isempty(refsites) ? inputs : reduce(vcat, (hpair(s) for s in refsites))
    return _Plan(sites, chans, outs, inputs, refs, pairs)
end

#---------- decimation and prewhitening -----

function _bessel_i0(x::Real)
    s, t, k = 1.0, 1.0, 0
    while true
        k += 1
        t *= (x / (2k))^2
        s += t
        t < 1.0e-16 * s && return s
    end
end

# A low-pass FIR for a decimation by `D`: a sinc with a Kaiser window, cut at
# half of the new Nyquist frequency, with a gain of 1 at zero frequency
function _decimation_fir(D::Integer = TF_DECIMATION; taps::Integer = TF_FIR_TAPS, beta::Real = TF_FIR_BETA)
    M = (taps - 1) / 2
    fc = 0.5 / D
    h = [begin
             x = n - M
             sinc_v = x == 0 ? 2fc : sin(2π * fc * x) / (π * x)
             sinc_v * _bessel_i0(beta * sqrt(max(0.0, 1 - (x / M)^2))) / _bessel_i0(beta)
         end for n in 0:(taps - 1)]
    return h ./ sum(h)
end

function _decimate(x::AbstractVector{Float64}, h::Vector{Float64}, D::Integer)
    L = length(h)
    n = length(x)
    m = n < L ? 0 : (n - L) ÷ D + 1
    y = Vector{Float64}(undef, m)
    @inbounds for i in 1:m
        s = 0.0
        base = (i - 1) * D
        @simd for j in 1:L
            s += h[j] * x[base + j]
        end
        y[i] = s
    end
    return y
end

# The coefficients `a` of an AR model x[t] ≈ Σ a[k] x[t-k], from up to 2^20
# samples of the spans (Levinson-Durbin with the biased autocorrelation). The
# correction for the filter is exact. Thus, the model only has to make the
# spectrum flatter
function _ar_model(spans, c::Integer, order::Integer)
    order <= 0 && return Float64[]
    budget = 1 << 20
    r = zeros(order + 1)
    used = 0
    for sp in spans
        used >= budget && break
        x = sp.cols[c]
        n = min(length(x), budget - used)
        n <= order + 1 && continue
        μ = sum(@view x[1:n]) / n
        for lag in 0:order
            s = 0.0
            @inbounds for t in (lag + 1):n
                s += (x[t] - μ) * (x[t - lag] - μ)
            end
            r[lag + 1] += s
        end
        used += n
    end
    (used == 0 || r[1] <= 0) && return Float64[]
    r ./= used
    a = zeros(order)
    e = r[1]
    for k in 1:order
        acc = r[k + 1]
        for j in 1:(k - 1)
            acc -= a[j] * r[k - j + 1]
        end
        κ = acc / e
        prev = copy(a)
        a[k] = κ
        for j in 1:(k - 1)
            a[j] = prev[j] - κ * prev[k - j]
        end
        e *= (1 - κ^2)
        e <= 0 && break
    end
    return a
end

_ar_response(a::Vector{Float64}, f::Real, fs::Real) =
    1 - sum((a[k] * cis(-2π * f * k / fs) for k in eachindex(a)); init = 0.0im)

#---------- bands -----

struct _Band
    period::Float64
    level::Int
    harmonics::Vector{Int}
end

# The bands on the fixed grid of periods 10^(j/bpd) s, each in the level where
# its centre is between 1/16 and 1/4 of the level rate. The first level goes
# up to `top` times its rate. The last level goes down to its harmonic
# `min_harmonic`
function _make_bands(fs, nfft, nlevels, bpd; min_period = 0.0, max_period = Inf, top = 0.25,
                     min_harmonic = 4)
    bands = _Band[]
    half = 10.0^(1 / (2bpd))
    for level in 0:(nlevels - 1)
        fl = fs / TF_DECIMATION^level
        hi = level == 0 ? top * fl : fl / 4
        lo = level == nlevels - 1 ? min_harmonic * fl / nfft : fl / 16
        jlo = ceil(Int, -bpd * log10(hi) - 1.0e-9)
        jhi = floor(Int, -bpd * log10(lo) + 1.0e-9)
        for j in jlo:jhi
            T = 10.0^(j / bpd)
            f = 1 / T
            (f <= hi * (1 + 1.0e-9) && f > lo * (1 - 1.0e-9)) || continue
            (level < nlevels - 1 && f <= fl / 16 * (1 + 1.0e-9)) && continue
            min_period <= T <= max_period || continue
            # no harmonic of the first level above the top frequency
            kmax = level == 0 ? min(floor(Int, top * nfft + 1.0e-9), nfft ÷ 2 - 1) : nfft ÷ 2 - 1
            k1 = max(ceil(Int, f / half * nfft / fl - 1.0e-9), 1)
            k2 = min(floor(Int, f * half * nfft / fl - 1.0e-9), kmax)
            # a band at the low harmonics of the last level holds one or two
            # harmonics. It takes the nearest ones up to TF_BAND_HARMONICS.
            # Thus, it has enough coefficients when the windows are few
            kc = f * nfft / fl
            kmin = max(min_harmonic, 1)
            k1 = max(k1, kmin)
            while k2 - k1 + 1 < TF_BAND_HARMONICS && (k1 > kmin || k2 < kmax)
                if k1 > kmin && (kc - (k1 - 1) <= (k2 + 1) - kc || k2 >= kmax)
                    k1 -= 1
                else
                    k2 += 1
                end
            end
            k1 <= k2 || continue
            push!(bands, _Band(T, level, collect(k1:k2)))
        end
    end
    return sort!(bands; by = b -> b.period)
end

#---------- Fourier coefficients of one level -----

struct _Span
    t0::Float64                          # seconds from 1970 of the first sample
    cols::Vector{_Col}
end

_span_windows(n, nfft, step, order) = n < nfft + order ? 0 : (n - nfft - order) ÷ step + 1

# The Fourier coefficients of each band at one level. A band holds a matrix
# with one row for each (window, harmonic) and one column for each channel,
# and the window index of each row
function _level_coefficients(plan::_Plan, spans::Vector{_Span}, bands::Vector{_Band}, level, fs_level,
                             nfft, step, order)
    nch = length(plan.chans)
    ks = sort!(unique(reduce(vcat, (b.harmonics for b in bands); init = Int[])))
    pos = Dict(k => i for (i, k) in enumerate(ks))
    nw = sum((_span_windows(length(sp.cols[1]), nfft, step, order) for sp in spans); init = 0)
    # the factor for each channel and harmonic: AR filter and sensor
    mult = Matrix{ComplexF64}(undef, length(ks), nch)
    ar = Vector{Vector{Float64}}(undef, nch)
    for c in 1:nch
        a = _ar_model(spans, c, order)
        ar[c] = vcat(a, zeros(order - length(a)))
        resp = plan.sites[plan.chans[c][1]].responses[plan.chans[c][2]]
        for (i, k) in enumerate(ks)
            f = k * fs_level / nfft
            mult[i, c] = 1 / (_ar_response(ar[c], f, fs_level) * _response(resp, f))
        end
    end
    store = [Matrix{ComplexF32}(undef, nw * length(b.harmonics), nch) for b in bands]
    times = Vector{Float64}(undef, nw)
    taper = _hann_window(nfft)
    tc = [i - (nfft + 1) / 2 for i in 1:nfft]
    stc = sum(abs2, tc)
    buf = Vector{Float64}(undef, nfft)
    spec = Vector{ComplexF64}(undef, nfft ÷ 2 + 1)
    fplan = plan_rfft(buf; flags = FFTW.ESTIMATE)
    coef = Matrix{ComplexF64}(undef, length(ks), nch)
    bpos = [[pos[k] for k in b.harmonics] for b in bands]
    w = 0
    for sp in spans
        n = length(sp.cols[1])
        nwin = _span_windows(n, nfft, step, order)
        for iw in 1:nwin
            s0 = order + (iw - 1) * step            # the window is s0+1 : s0+nfft
            w += 1
            times[w] = sp.t0 + (s0 + nfft / 2) / fs_level
            for c in 1:nch
                x = sp.cols[c]
                a = ar[c]
                μ = 0.0
                @inbounds for i in 1:nfft
                    t = s0 + i
                    v = x[t]
                    for k in 1:order
                        v -= a[k] * x[t - k]
                    end
                    buf[i] = v
                    μ += v
                end
                μ /= nfft
                slope = 0.0
                @inbounds for i in 1:nfft
                    slope += tc[i] * buf[i]
                end
                slope /= stc
                @inbounds for i in 1:nfft
                    buf[i] = (buf[i] - μ - slope * tc[i]) * taper[i]
                end
                mul!(spec, fplan, buf)
                @inbounds for (i, k) in enumerate(ks)
                    coef[i, c] = spec[k + 1] * mult[i, c]
                end
            end
            for (i, j, R) in plan.pairs
                @inbounds for h in eachindex(ks)
                    u, v = coef[h, i], coef[h, j]
                    coef[h, i] = R[1, 1] * u + R[1, 2] * v
                    coef[h, j] = R[2, 1] * u + R[2, 2] * v
                end
            end
            for (b, hp) in enumerate(bpos)
                M = store[b]
                row0 = (w - 1) * length(hp)
                @inbounds for (h, p) in enumerate(hp), c in 1:nch
                    M[row0 + h, c] = ComplexF32(coef[p, c])
                end
            end
        end
    end
    return store, times
end

function _decimate_spans(spans::Vector{_Span}, h, fs_level)
    out = _Span[]
    delay = (length(h) - 1) / 2 / fs_level
    for sp in spans
        length(sp.cols[1]) < length(h) && continue
        push!(out, _Span(sp.t0 + delay, [_whole(_decimate(x, h, TF_DECIMATION)) for x in sp.cols]))
    end
    return out
end

#---------- robust regression -----

_huber(x, k) = x <= k ? 1.0 : k / x

# A redescending weight: 1 near zero, 1/e at x = β, and fast to zero above.
# The form is the one of Chave, Thomson & Ander (1987) and Chave & Thomson
# (2004)
_thomson(x, β) = exp(exp(-β^2) - exp(β * (x - β)))

# The robust scale of complex residuals: for a circular Gaussian with
# E|r|² = σ², the median of |r| is σ·√ln 2
function _robust_scale(r, prior)
    a = Float64[abs(r[i]) for i in eachindex(r) if prior[i] > 0]
    isempty(a) && return 0.0
    return median!(a) / sqrt(log(2))
end

# The (1 - 1/n) quantile of a Gamma(p, 1) variable: the largest normalized
# leverage that n points of p circular Gaussian inputs would give
function _gamma_quantile(p::Integer, n::Real)
    target = 1 / max(n, 2)
    surv(t) = exp(-t) * sum(t^k / factorial(k) for k in 0:(p - 1))
    lo, hi = 0.0, 10.0 + 4p
    while surv(hi) > target
        hi *= 2
    end
    for _ in 1:80
        mid = (lo + hi) / 2
        surv(mid) > target ? (lo = mid) : (hi = mid)
    end
    return (lo + hi) / 2
end

function _solve_weighted(y, X, Xh, w)
    Xw = X .* w
    A = Xh' * Xw
    b = Xh' * (w .* y)
    return A \ b
end

"""
    _robust_fit(y, X, Xh, prior; huber = 1.5, leverage = true) -> (β, w, scale, converged)

A robust instrumental variable fit of `y ≈ X β`, with the instruments `Xh`
(`Xh = X` gives a least squares fit). `prior` are weights from a previous
stage, multiplied into each step. The steps:
1. a weighted fit with `prior`;
2. Huber weights on the residuals scaled by a robust scale, to convergence;
3. a redescending weight with a threshold at the largest residual that a
   circular Gaussian gives for this number of points, and a weight of the
   same form on the normalized leverage of each point.
"""
function _robust_fit(y::AbstractVector, X::AbstractMatrix, Xh::AbstractMatrix, prior::Vector{Float64};
                     huber::Real = 1.5, leverage::Bool = true, redescend::Bool = true, maxit::Integer = 40,
                     tol::Real = 1.0e-5)
    n, p = size(X)
    w = copy(prior)
    β = _solve_weighted(y, X, Xh, w)
    r = similar(y, ComplexF64)
    m = count(>(0), prior)
    βres = sqrt(log(max(m, 3)))
    βlev = sqrt(_gamma_quantile(p, max(m, 3)))
    scale = 0.0
    converged = false
    lev = ones(n)
    for stage in (redescend ? (:huber, :redescend) : (:huber,))
        for it in 1:maxit
            mul!(r, X, β)
            @. r = y - r
            scale = _robust_scale(r, prior)
            scale > 0 || break
            if stage === :redescend && leverage
                Xw = Xh .* w
                C = Hermitian(Xh' * Xw)
                Q = Xh * inv(Matrix(C))
                sw = sum(w)
                @inbounds for i in 1:n
                    g = 0.0
                    for j in 1:p
                        g += real(conj(Xh[i, j]) * Q[i, j])
                    end
                    lev[i] = _thomson(sqrt(max(g * sw, 0.0)), βlev)
                end
            end
            @inbounds for i in 1:n
                x = abs(r[i]) / scale
                wr = stage === :huber ? _huber(x, huber) : _thomson(x, βres) * lev[i]
                w[i] = prior[i] * wr
            end
            βn = _solve_weighted(y, X, Xh, w)
            δ = maximum(abs, βn .- β) / max(maximum(abs, βn), eps())
            β = βn
            if δ < tol
                converged = stage === :redescend || !redescend
                break
            end
            stage === :redescend && it >= 10 && (converged = true; break)
        end
    end
    return β, w, scale, converged
end

# The weighted multiple coherence of y with its prediction
function _coherence(y, pred, w)
    a = sum(w[i] * conj(y[i]) * pred[i] for i in eachindex(y); init = 0.0im)
    b = sum(w[i] * abs2(y[i]) for i in eachindex(y); init = 0.0)
    c = sum(w[i] * abs2(pred[i]) for i in eachindex(y); init = 0.0)
    return b > 0 && c > 0 ? sqrt(clamp(abs2(a) / (b * c), 0.0, 1.0)) : NaN
end

# The sums of one group of rows for the jackknife: Σ w Rᴴ X and Σ w Rᴴ y
function _group_sums(R, X, y, w, groups, G)
    q, p = size(R, 2), size(X, 2)
    SX = [zeros(ComplexF64, q, p) for _ in 1:G]
    Sy = [zeros(ComplexF64, q) for _ in 1:G]
    @inbounds for i in eachindex(groups)
        wi = w[i]
        wi == 0 && continue
        g = groups[i]
        for a in 1:q
            ra = conj(R[i, a]) * wi
            for b in 1:p
                SX[g][a, b] += ra * X[i, b]
            end
            Sy[g][a] += ra * y[i]
        end
    end
    return SX, Sy
end

function _jackknife_cov(estimates::Vector{Vector{ComplexF64}})
    G = length(estimates)
    G < 3 && return nothing
    μ = sum(estimates) / G
    p = length(μ)
    C = zeros(ComplexF64, p, p)
    for e in estimates
        d = e - μ
        C .+= d * d'
    end
    return C .* ((G - 1) / G)
end

#---------- one band -----

struct _BandResult
    rows::Dict{Symbol, Vector{ComplexF64}}       # :ex, :ey, :hz => [b_x, b_y]
    var::Dict{Symbol, Vector{Float64}}
    coh::Dict{Symbol, Float64}
    ref_coh::Vector{Float64}
    n_windows::Int
    n_obs::Int
end

function _estimate_band(plan::_Plan, M::Matrix{ComplexF32}, wid::Vector{Int}, nwin::Int, opts)
    n = size(M, 1)
    X = ComplexF64.(M[:, plan.inputs])
    p = size(X, 2)
    reference = plan.refs != plan.inputs
    G = min(opts.jackknife_groups, nwin)
    groups = [1 + ((wid[i] - 1) * G) ÷ nwin for i in 1:n]
    u = ones(n)
    ref_coh = fill(NaN, p)
    ct = opts.method === :ct2004
    if reference
        R = ComplexF64.(M[:, plan.refs])
        for j in 1:p
            xj = X[:, j]
            if ct
                # stage 1: a robust fit of the input to the reference. Its
                # weights remove noise bursts of the local or the reference H
                βj, wj, _, _ = _robust_fit(xj, R, R, ones(n); huber = opts.huber, leverage = opts.leverage)
                u .= min.(u, wj)
                ref_coh[j] = _coherence(xj, R * βj, wj)
            else
                # the reference as it is: a plain least squares projection
                ref_coh[j] = _coherence(xj, R * ((R' * R) \ (R' * xj)), u)
            end
        end
        S1RX, _ = _group_sums(R, X, X[:, 1], u, groups, G)
        S1RR, _ = _group_sums(R, R, R[:, 1], u, groups, G)
        B = sum(S1RR) \ sum(S1RX)
        Xh = R * B
    else
        R = X
        Xh = X
    end
    rows = Dict{Symbol, Vector{ComplexF64}}()
    vars = Dict{Symbol, Vector{Float64}}()
    cohs = Dict{Symbol, Float64}()
    for (name, col) in plan.outputs
        y = ComplexF64.(M[:, col])
        β, w, _, _ = ct ? _robust_fit(y, X, Xh, u; huber = opts.huber, leverage = opts.leverage) :
                     _robust_fit(y, X, Xh, u; huber = opts.huber, leverage = false, redescend = false)
        pred = X * β
        cohs[name] = _coherence(y, pred, w)
        # the jackknife: delete one block of windows, refit with the final weights
        S2RX, S2Ry = _group_sums(R, X, y, w, groups, G)
        tRX, tRy = sum(S2RX), sum(S2Ry)
        if reference
            tRR, tR1 = sum(S1RR), sum(S1RX)
        end
        est = Vector{ComplexF64}[]
        for g in 1:G
            Bg = reference ? (tRR - S1RR[g]) \ (tR1 - S1RX[g]) : Matrix{ComplexF64}(I, p, p)
            A = Bg' * (tRX - S2RX[g])
            b = Bg' * (tRy - S2Ry[g])
            cond_ok = all(isfinite, A) && abs(det(A)) > 0
            cond_ok && push!(est, A \ b)
        end
        C = length(est) == G ? _jackknife_cov(est) : nothing
        # the sandwich (heteroscedastic) covariance of the weighted estimate
        r = y .- pred
        Ginv = inv(Xh' * (X .* w))
        Mid = Xh' * (Xh .* (w .^ 2 .* abs2.(r)))
        S = Ginv * Mid * Ginv'
        v_jk = C === nothing ? fill(NaN, p) : real.(diag(C))
        v_sw = real.(diag(S))
        # the jackknife. Where it has too few blocks, the sandwich
        v = [isfinite(a) ? a : b for (a, b) in zip(v_jk, v_sw)]
        rows[name] = β
        vars[name] = v
    end
    return _BandResult(rows, vars, cohs, ref_coh, nwin, n)
end

#---------- one rate -----

function _estimate_rate(inputs, fs, opts; progress)
    target, insrc, refsrcs = inputs
    tcomps = _components_of(target)
    outs = Symbol[c for c in (:e1, :e2) if c in tcomps]
    :bz in tcomps && push!(outs, :bz)
    isempty(outs) && error("$(_site_label(target)) has no Ex, Ey or Hz to estimate")
    base_mode = insrc !== nothing
    target_comps = base_mode ? outs : unique(vcat(outs, [:bx, :by]))
    load(x, comps) = _load_site(x, fs, comps; calibration = opts.calibration, dipole = opts.dipole,
                                azimuths = opts.azimuths, masks = opts.masks, span = opts.span)
    progress("Loading $(_site_label(target)) at $(_fs_label(fs)): $(_channel_names(target_comps))")
    sites = _SiteData[load(target, target_comps)]
    insite = 1
    if base_mode
        progress("Loading base $(_site_label(insrc)): Hx Hy")
        push!(sites, load(insrc, [:bx, :by]))
        insite = 2
    end
    refsites = Int[]
    for r in refsrcs
        progress("Loading remote $(_site_label(r)): Hx Hy")
        push!(sites, load(r, [:bx, :by]))
        push!(refsites, length(sites))
    end
    plan = _make_plan(sites, 1, insite, refsites, outs)
    results, run = _run_plan(plan, sites, fs, opts; progress)
    lat, lon, elev = sites[1].latitude, sites[1].longitude, sites[1].elevation
    info = Dict{Symbol, Any}(
        :hours => run.hours, :levels => run.levels, :spans => run.spans, :start => run.start, :stop => run.stop,
        :nyquist_fraction => run.nyquist_fraction, :edge => run.edge,
        :dipoles => Dict(c => 1000 * r.gain for (c, r) in sites[1].responses if r.kind === :dipole),
        # the directions of the channels of the target (and of the inputs)
        :azimuths => Dict(c => r.azimuth for (c, r) in sites[1].responses),
        :input_azimuths => Dict(c => sites[insite].responses[c].azimuth for c in (:bx, :by)),
        :sensors => Dict(sites[s].name * "." * _channel_names([c]) => sites[s].responses[c].note for (s, c) in plan.chans),
        :files => sort!(unique(String[sites[s].responses[c].file for (s, c) in plan.chans if !isempty(sites[s].responses[c].file)])),
    )
    return results, (lat, lon, elev), info
end

# Where the recorder's anti-alias filter ends the signal: the lowest frequency
# above half the Nyquist frequency where the coherence of the outputs (E) with
# the inputs (H) falls below TF_EDGE_COHERENCE. Each output takes the input of
# the best pair at each frequency, because the field changes its direction.
# The coherence needs no calibration: the two channels of a pair pass the same
# filters. At most 2000 windows of 1024 samples, spread over the spans. NaN if
# the coherence does not fall before the Nyquist frequency
function _coherence_edge(plan::_Plan, spans::Vector{_Span}, fs; nfft::Integer = 1024)
    outs = [c for (k, c) in plan.outputs if k in (:ex, :ey)]
    ins = plan.inputs
    isempty(outs) && return NaN
    nf = nfft ÷ 2 + 1
    chans = unique(vcat(outs, ins))
    S = Dict(c => zeros(nf) for c in chans)
    C = Dict((o, i) => zeros(ComplexF64, nf) for o in outs for i in ins)
    starts = [(sp, i0) for sp in spans for i0 in 1:nfft:(length(sp.cols[1]) - nfft + 1)]
    isempty(starts) && return NaN
    stride = max(1, length(starts) ÷ 2000)
    taper = _hann_window(nfft)
    buf = Vector{Float64}(undef, nfft)
    fplan = plan_rfft(buf; flags = FFTW.ESTIMATE)
    F = Dict(c => Vector{ComplexF64}(undef, nf) for c in chans)
    for (sp, i0) in starts[1:stride:end]
        for c in chans
            x = view(sp.cols[c], i0:(i0 + nfft - 1))
            μ = sum(x) / nfft
            @. buf = (x - μ) * taper
            mul!(F[c], fplan, buf)
            S[c] .+= abs2.(F[c])
        end
        for (o, i) in keys(C)
            C[(o, i)] .+= F[o] .* conj.(F[i])
        end
    end
    coh = [minimum(maximum(abs2(C[(o, i)][k]) / (S[o][k] * S[i][k]) for i in ins) for o in outs) for k in 1:nf]
    # a median over five bins: one narrow line or one noisy bin does not move
    # the edge
    smooth = [median(coh[max(1, k - 2):min(nf, k + 2)]) for k in 1:nf]
    k0 = nfft ÷ 4 + 1
    k = findfirst(j -> j >= k0 && !(smooth[j] >= TF_EDGE_COHERENCE), 1:nf)
    return k === nothing ? NaN : (k - 1) * fs / nfft
end

# The estimate of a plan at one rate: the spans that all the sites recorded,
# the decimation levels, the bands and the robust fit of each band
function _run_plan(plan::_Plan, sites::Vector{_SiteData}, fs, opts; progress)
    common = _common_spans(sites)
    nfft = opts.window
    step = max(1, round(Int, nfft * (1 - opts.overlap)))
    order = opts.prewhiten
    spans = _Span[]
    for (lo, hi, segs) in common
        hi - lo + 1 >= nfft + order || continue
        cols = _Col[]
        for (s, comp) in plan.chans
            seg = sites[s].segments[segs[s]]
            push!(cols, _subview(seg.cols[comp], (lo - seg.k0 + 1):(hi - seg.k0 + 1)))
        end
        push!(spans, _Span(lo / fs, cols))
    end
    isempty(spans) && error("The sites never recorded together for one window ($(nfft) samples at $(_fs_label(fs)))")
    hours = sum(length(sp.cols[1]) for sp in spans) / fs / 3600
    first_t = _datetime(spans[1].t0)
    last_t = _datetime(spans[end].t0 + length(spans[end].cols[1]) / fs)
    progress(@sprintf("%.1f h recorded together in %d span%s, %s to %s", hours, length(spans),
                      length(spans) == 1 ? "" : "s", Dates.format(first_t, "yyyy-mm-dd HH:MM"),
                      Dates.format(last_t, "yyyy-mm-dd HH:MM")))
    # the levels that have enough windows
    lengths = [length(sp.cols[1]) for sp in spans]
    L = TF_FIR_TAPS
    nlevels = 0
    while nlevels < opts.max_levels
        sum((_span_windows(n, nfft, step, order) for n in lengths); init = 0) >= opts.min_windows || break
        nlevels += 1
        lengths = [n >= L ? (n - L) ÷ TF_DECIMATION + 1 : 0 for n in lengths]
    end
    nlevels == 0 && error("Too little data at $(_fs_label(fs)) for $(opts.min_windows) windows")
    nyq = opts.nyquist_fraction
    edge = NaN
    if nyq === :auto
        progress("finding the top frequency from the coherence of E and H")
        edge = _coherence_edge(plan, spans, fs)
        # a coherence that does not fall: no filter edge below the Nyquist
        # frequency
        nyq = isfinite(edge) ? clamp(TF_EDGE_MARGIN * edge / (fs / 2), 0.5, 0.9) : 0.9
        progress(isfinite(edge) ? @sprintf("filter edge at %.3g Hz: top frequency %.3g Hz", edge, nyq * fs / 2) :
                                  @sprintf("no filter edge below Nyquist: top frequency %.3g Hz", nyq * fs / 2))
    end
    bands = _make_bands(fs, nfft, nlevels, opts.bands_per_decade;
                        min_period = opts.min_period, max_period = opts.max_period,
                        top = nyq / 2, min_harmonic = opts.min_harmonic)
    isempty(bands) || progress(@sprintf("%d decimation level%s, %d bands from %.3g s to %.3g s", nlevels,
                                        nlevels == 1 ? "" : "s", length(bands), minimum(b.period for b in bands),
                                        maximum(b.period for b in bands)))
    h = _decimation_fir()
    results = Dict{Float64, _BandResult}()
    for level in 0:(nlevels - 1)
        fl = fs / TF_DECIMATION^level
        lb = [b for b in bands if b.level == level]
        if !isempty(lb)
            stage = @sprintf("level %d of %d (%s)", level + 1, nlevels, _fs_label(fl))
            progress("$stage: Fourier coefficients")
            store, times = _level_coefficients(plan, spans, lb, level, fl, nfft, step, order)
            nwin = length(times)
            if nwin >= opts.min_windows
                out = Vector{Union{Nothing, _BandResult}}(nothing, length(lb))
                done = Threads.Atomic{Int}(0)
                progress(@sprintf("%s: %d windows, robust fit of %d bands", stage, nwin, length(lb)))
                Threads.@threads for b in eachindex(lb)
                    nh = length(lb[b].harmonics)
                    wid = [1 + (i - 1) ÷ nh for i in 1:size(store[b], 1)]
                    size(store[b], 1) >= 2 * length(plan.refs) + 4 || continue
                    out[b] = try
                        _estimate_band(plan, store[b], wid, nwin, opts)
                    catch err
                        err isa LinearAlgebra.SingularException || err isa LinearAlgebra.LAPACKException || rethrow()
                        nothing
                    end
                    k = Threads.atomic_add!(done, 1) + 1
                    progress(@sprintf("%s: band %d of %d done (%.3g s)", stage, k, length(lb), lb[b].period))
                end
                for (b, r) in zip(lb, out)
                    r === nothing || (results[b.period] = r)
                end
            end
            store = nothing
        end
        if level < nlevels - 1
            progress(@sprintf("decimating to %s", _fs_label(fl / TF_DECIMATION)))
            spans = _decimate_spans(spans, h, fl)
        end
    end
    return results, (hours = hours, levels = nlevels, spans = length(common), start = first_t, stop = last_t,
                     nyquist_fraction = nyq, edge = edge)
end

#---------- the public entry points -----

"""
    estimate_tf(site; base = nothing, remote = nothing, rate = :all, kwargs...) -> TransferFunction
    estimate_tf(survey::Survey, site; base = nothing, remote = nothing, kwargs...) -> TransferFunction

Estimate the impedance and the tipper of `site`. A site is a site directory
(Metronix, LEMI-424 or GEOMAG), a data file, a [`SurveySite`](@ref), a
[`TimekeeperRun`](@ref) or a vector of runs. With a [`Survey`](@ref), give
the sites by their names; `:auto` then selects the base or the remote site
with the longest overlap (refer to [`site_references`](@ref)).

The set-up:
- no `base`, no `remote` -- single site: the inputs are Hx, Hy of the site;
- `base` -- the inputs are Hx, Hy of the base site. Use it for a site without
  magnetic channels, or with poor ones;
- `remote` -- one site or a vector of sites. Their Hx, Hy are the reference.
  In a vector of runs, the runs of each site name are one remote site.
  With more than one remote site, the reference has all their channels (a
  two-stage least squares fit);
- `base` and `remote` together.

The method is in the header of `Processing.jl`. The keywords:
- `rate` -- the sampling rate in Hz, or `:all` (each rate that all the sites
  recorded, then the best estimate at each period).
- `window = 256` -- the section length: samples in each data section
  (Chave & Thomson 2004), at each decimation level.
- `overlap = 0.5` -- the section overlap.
- `bands_per_decade = 8` -- the bands at the periods `10^(j/8)` s.
- `prewhiten = 3` -- the terms of the least-squares AR pre-whitening filter (0:
  none).
- `min_windows = 8` -- the data sections that a level needs. Fewer gives
  longer periods with larger errors.
- `min_harmonic = 4` -- the lowest frequency index of a section at the last
  level. A lower one gives longer periods, closer to the trend of each section.
- `nyquist_fraction = :auto` -- the highest frequency as a fraction of the
  Nyquist frequency of the record (0.5 is a quarter of the rate). `:auto`
  finds where the anti-alias filter of the recorder ends the signal: the
  lowest frequency above half the Nyquist frequency where the coherence of E
  with H falls below 0.5. The top is 5 % below it, from 0.5 to 0.9. A number
  sets it. A higher one gives shorter periods, closer to the filter.
- `max_levels = 12` -- the most decimation levels.
- `min_period = 0`, `max_period = Inf` -- the periods to estimate.
- `method = :ct2004` -- `:ct2004` (two stages with bounded influence, after
  Chave & Thomson 2004) or `:eb1986` (one robust M-estimate with Huber
  weights, after Egbert & Booker 1986). Step 5 of the header of
  Processing.jl compares them.
- `huber = 1.5` -- the parameter a of the Huber weights (Chave & Thomson
  2004, eq. 7).
- `leverage = true` -- the leverage weights of the bounded influence
  estimator (Chave & Thomson 2004, eq. 26; `:ct2004` only).
- `jackknife_groups = 50` -- the blocks of data sections that the jackknife
  deletes one at a time.
- `calibration` -- a directory (or a vector of directories) with the coil
  calibration files. If you do not give it, the function looks in the site
  directory and up to three directories above it, and in the `s` subdirectory
  of each of those (the calibration directory of a survey is `s/` of its root).
- `dipole` -- `Dict(:e1 => L1, :e2 => L2)` in metres. If you do not give it,
  the lengths come from the row of the site in `d/dipoles.dat` of the survey
  (refer to [`read_dipoles`](@ref)), then from the electrode positions in the
  `.ats` headers. Without either, each electrode is 50 m from the centre.
- `azimuths` -- `Dict(component => degrees east of north)` to replace the
  directions of the sensors.
- `masks` -- `Dict(site name => [(start, stop), …])`: intervals to leave out.
- `span` -- `(start, stop)`: use only this time span.
- `progress` -- a function that receives a line of text at each step. The
  function can call it from more than one thread.
"""
function estimate_tf(site; base = nothing, remote = nothing, rate = :all, window::Integer = 256,
                     overlap::Real = 0.5, bands_per_decade::Integer = 8, prewhiten::Integer = 3,
                     min_windows::Integer = 8, min_harmonic::Integer = 4, nyquist_fraction = :auto,
                     max_levels::Integer = 12, min_period::Real = 0.0,
                     max_period::Real = Inf, method::Symbol = :ct2004, huber::Real = 1.5, leverage::Bool = true,
                     jackknife_groups::Integer = 50,
                     calibration = nothing, dipole = Dict{Symbol, Float64}(),
                     azimuths = Dict{Symbol, Float64}(), masks = nothing, span = nothing,
                     progress = nothing)
    ispow2(window) || error("window must be a power of 2, got $window")
    window >= 64 || error("window must be 64 samples or more")
    0 <= overlap < 1 || error("overlap must be in [0, 1)")
    method in (:ct2004, :eb1986) || error("method must be :ct2004 or :eb1986")
    nyquist_fraction === :auto || 0.1 <= nyquist_fraction <= 0.9 ||
        error("nyquist_fraction must be :auto or in [0.1, 0.9]")
    1 <= min_harmonic < window ÷ 8 || error("min_harmonic must be in [1, $(window ÷ 8 - 1)]")
    min_windows >= 3 || error("min_windows must be 3 or more")
    say = progress === nothing ? (_ -> nothing) : progress
    target = _as_site(site)
    insrc = base === nothing ? nothing : _as_site(base)
    refsrcs = remote === nothing ? Any[] :
              remote isa AbstractVector{TimekeeperRun} ? _runs_by_site(remote) :
              remote isa Union{TimekeeperRun, AbstractString, SurveySite} ?
              Any[_as_site(remote)] : Any[_as_site(r) for r in remote]
    name = _site_label(target)
    insrc === nothing && !(:bx in _components_of(target) && :by in _components_of(target)) &&
        error("$name has no Hx and Hy: give a base site for the inputs")
    insrc === nothing || (:bx in _components_of(insrc) && :by in _components_of(insrc)) ||
        error("The base site $(_site_label(insrc)) has no Hx and Hy")
    for r in refsrcs
        :bx in _components_of(r) && :by in _components_of(r) ||
            error("The remote site $(_site_label(r)) has no Hx and Hy")
    end
    # the rates that every site gives: its own, or a decimation of a faster one
    # (a 10 Hz site with 1 Hz sites is processed at 1 Hz)
    sources = vcat(Any[target], insrc === nothing ? Any[] : Any[insrc], refsrcs)
    rates = sort!(unique(fs for x in sources for fs in _rates_of(x)); rev = true)
    rates = [fs for fs in rates if all(x -> _gives_rate(x, fs), sources)]
    if rate !== :all
        rates = [fs for fs in rates if _same_rate(fs, rate)]
        isempty(rates) && error("Not all the sites recorded at $(_fs_label(rate))")
    end
    isempty(rates) && error("The sites have no sampling rate in common")
    opts = (; window = Int(window), overlap = Float64(overlap), bands_per_decade = Int(bands_per_decade),
            prewhiten = Int(prewhiten), min_windows = Int(min_windows), min_harmonic = Int(min_harmonic),
            nyquist_fraction = nyquist_fraction === :auto ? :auto : Float64(nyquist_fraction), max_levels = Int(max_levels),
            min_period = Float64(min_period), max_period = Float64(max_period), huber = Float64(huber),
            method, leverage, jackknife_groups = Int(jackknife_groups), calibration,
            dipole = Dict{Symbol, Float64}(Symbol(k) => Float64(v) for (k, v) in dipole),
            azimuths = Dict{Symbol, Float64}(Symbol(k) => Float64(v) for (k, v) in azimuths),
            masks, span)
    per_rate = Dict{Float64, Any}()
    pos = (NaN, NaN, NaN)
    infos = Dict{Float64, Any}()
    failures = String[]
    say("$(length(rates)) rate$(length(rates) == 1 ? "" : "s") in common: " * join(_fs_label.(sort(rates; rev = true)), ", "))
    for (k, fs) in enumerate(sort(rates; rev = true))
        length(rates) > 1 && say("Rate $k of $(length(rates)): $(_fs_label(fs))")
        try
            res, p, info = _estimate_rate((target, insrc, refsrcs), fs, opts; progress = say)
            per_rate[fs] = res
            infos[fs] = info
            isfinite(pos[1]) || (pos = p)
        catch err
            length(rates) == 1 && rethrow()
            push!(failures, "$(_fs_label(fs)): $(sprint(showerror, err))")
            @warn "No estimate at $(_fs_label(fs)) for $name" exception = err
        end
    end
    isempty(per_rate) && error("No estimate for $name:\n" * join(failures, "\n"))
    mode = insrc === nothing ? (isempty(refsrcs) ? :single : :remote) :
           (isempty(refsrcs) ? :base : :base_remote)
    md = Dict{Symbol, Any}(:options => opts, :rate => rate, :rates => infos, :failures => failures,
                           :start => minimum(i[:start] for i in values(infos)),
                           :stop => maximum(i[:stop] for i in values(infos)),
                           :dipoles => first(values(infos))[:dipoles])
    tf = _assemble(name, mode, insrc === nothing ? "" : _site_label(insrc),
                   String[_site_label(r) for r in refsrcs], pos, per_rate, md)
    say("Done: $name")
    return tf
end

# The relative error of a band: the median of var/|Z|² of the elements that
# it has, used to select a rate at a period
function _relative_error(r::_BandResult)
    v = Float64[]
    for name in (:ex, :ey)
        haskey(r.rows, name) || continue
        for j in 1:2
            a = abs2(r.rows[name][j])
            a > 0 && push!(v, r.var[name][j] / a)
        end
    end
    return isempty(v) ? Inf : median(v)
end

function _assemble(name, mode, base, remote, pos, per_rate, md)
    periods = sort!(unique(reduce(vcat, (collect(keys(r)) for r in values(per_rate)); init = Float64[])))
    n = length(periods)
    Z = fill(complex(NaN, NaN), 2, 2, n)
    ZV = fill(NaN, 2, 2, n)
    T = fill(complex(NaN, NaN), 2, n)
    TV = fill(NaN, 2, n)
    coh = fill(NaN, 3, n)
    rcoh = fill(NaN, 2, n)
    nw = zeros(Int, n)
    no = zeros(Int, n)
    rates = fill(NaN, n)
    for (k, T0) in enumerate(periods)
        best, bestfs, beste = nothing, NaN, Inf
        for (fs, res) in per_rate
            r = get(res, T0, nothing)
            r === nothing && continue
            e = _relative_error(r)
            if best === nothing || e < beste
                best, bestfs, beste = r, fs, e
            end
        end
        r = best
        for (i, row) in ((1, :ex), (2, :ey))
            haskey(r.rows, row) || continue
            Z[i, :, k] = r.rows[row]
            ZV[i, :, k] = r.var[row]
            coh[i, k] = r.coh[row]
        end
        if haskey(r.rows, :hz)
            T[:, k] = r.rows[:hz]
            TV[:, k] = r.var[:hz]
            coh[3, k] = r.coh[:hz]
        end
        rcoh[:, k] = r.ref_coh
        nw[k], no[k], rates[k] = r.n_windows, r.n_obs, bestfs
    end
    return TransferFunction(name, mode, base, remote, pos..., periods, Z, ZV, T, TV, coh, rcoh, nw, no, rates, md)
end

# Runs as sites: the runs of each site name are one site
function _runs_by_site(runs::AbstractVector{TimekeeperRun})
    names = unique(r.site for r in runs)
    return Any[TimekeeperRun[r for r in runs if r.site == n] for n in names]
end

function _survey_source(s::Survey, x)
    x === nothing && return nothing
    x isa AbstractString && return s[x]
    x isa SurveySite && return x
    x isa Integer && return s.sites[x]
    return x
end

function estimate_tf(s::Survey, site; base = nothing, remote = nothing, kwargs...)
    t = _target_site(s, site)
    refs = nothing
    if base === :auto || remote === :auto
        refs = site_references(s, t)
    end
    base === :auto && (base = isempty(refs.base) ? nothing : first(refs.base).site)
    remote === :auto && (remote = _pick_remote(refs.remote, base isa SurveySite ? base.name : base))
    b = _survey_source(s, base)
    r = remote === nothing ? nothing :
        remote isa Union{AbstractString, SurveySite, Integer} ? _survey_source(s, remote) :
        [_survey_source(s, x) for x in remote]
    return estimate_tf(t; base = b, remote = r, kwargs...)
end

"""
    default_references(survey, site; base_km = 5.0, remote_km = 0.1) -> (base, remote)

The set-up that TKProc selects when it loads `site`: names or `nothing`.
- A site with Hx and Hy: no base site, and the remote site (none if no site
  is far enough).
- A site without them: the base site with the longest overlap, and the
  remote site.
The remote site is the farthest of the sites that recorded with `site` at
least 90 % as long as the best one: its local noise is the least like that
of the site.
"""
function default_references(s::Survey, site; base_km::Real = 5.0, remote_km::Real = 0.1)
    t = _target_site(s, site)
    refs = site_references(s, t; base_km, remote_km)
    base = has_magnetic(t) || isempty(refs.base) ? nothing : first(refs.base).site
    return (base = base, remote = _pick_remote(refs.remote, base))
end

function _pick_remote(candidates, base)
    pool = [c for c in candidates if c.site != base]
    isempty(pool) && return nothing
    best = maximum(c.overlap_hours for c in pool)
    near_best = [c for c in pool if c.overlap_hours >= 0.9best]
    return near_best[argmax([c.distance_km for c in near_best])].site
end
