# Spectra.jl - power spectral density and spectrogram estimation.
# Author: @pankajkmishra
#
# Implements Welch's method (averaged over overlapping segments, or over a set
# of disjoint good segments) and the STFT used by the spectrogram view, which
# blanks any window touching masked samples. All of it runs through a reusable
# SpectralWorkspace holding the window, the FFT plan and its scratch buffers,
# so repeated estimates at one configuration allocate nothing extra.

using FFTW
using LinearAlgebra: mul!

function _hann_window(n::Integer)
    n <= 1 && return ones(Float64, max(n, 0))
    return Float64[0.5 * (1 - cos(2π * (k - 1) / (n - 1))) for k in 1:n]
end

function _make_window(kind::Symbol, n::Integer)
    (kind === :hann || kind === :hanning) && return _hann_window(n)
    kind === :rect && return ones(Float64, n)
    error("Unsupported window: $kind")
end

mutable struct SpectralWorkspace{P}
    nfft::Int
    noverlap::Int
    fs::Float64
    window::Symbol
    detrend::Symbol
    weights::Vector{Float64}
    norm_psd::Float64
    freqs::Vector{Float64}
    buf_x::Vector{Float64}
    fft_x::Vector{ComplexF64}
    plan::P
end

"""
Concrete plan type produced by every [`SpectralWorkspace`](@ref). Real-FFT plan
types do not depend on transform length, so one workspace type covers every
`nfft` and lets caches hold workspaces concretely instead of as `Any`.
"""
const SpectralPlanType = typeof(plan_rfft(Vector{Float64}(undef, 2); flags = FFTW.ESTIMATE))

"""
Concretely typed [`SpectralWorkspace`](@ref), for cache and container fields.
"""
const TKSpectralWorkspace = SpectralWorkspace{SpectralPlanType}

function SpectralWorkspace(
    nfft::Integer,
    fs::Real;
    noverlap::Integer = nfft ÷ 2,
    window::Symbol = :hann,
    detrend::Symbol = :mean,
)
    nfft_i = Int(nfft)
    noverlap_i = Int(noverlap)
    nfft_i > 0 || error("nfft must be positive")
    0 <= noverlap_i < nfft_i || error("noverlap must be >= 0 and < nfft")
    fs_f = Float64(fs)
    fs_f > 0 || error("Sample rate must be positive, got $fs")
    weights = _make_window(window, nfft_i)
    sum_w_sq = sum(abs2, weights)
    norm_psd = 1.0 / (fs_f * sum_w_sq)
    freqs = Float64[(i - 1) * fs_f / nfft_i for i in 1:(nfft_i ÷ 2 + 1)]
    n_freqs = length(freqs)
    buf_x = Vector{Float64}(undef, nfft_i)
    fft_x = Vector{ComplexF64}(undef, n_freqs)
    plan = plan_rfft(buf_x; flags = FFTW.ESTIMATE)
    return SpectralWorkspace{typeof(plan)}(
        nfft_i,
        noverlap_i,
        fs_f,
        window,
        detrend,
        weights,
        norm_psd,
        freqs,
        buf_x,
        fft_x,
        plan,
    )
end

function _validate_workspace(
    ws::SpectralWorkspace,
    fs::Real,
    nfft::Integer,
    noverlap::Integer,
    window::Symbol,
    detrend::Symbol,
)
    ws.nfft == Int(nfft) || error("SpectralWorkspace nfft mismatch")
    ws.noverlap == Int(noverlap) || error("SpectralWorkspace noverlap mismatch")
    ws.fs == Float64(fs) || error("SpectralWorkspace sample-rate mismatch")
    ws.window === window || error("SpectralWorkspace window mismatch")
    ws.detrend === detrend || error("SpectralWorkspace detrend mismatch")
    return ws
end

function _detrend_segment!(seg::AbstractVector{<:Real}, kind::Symbol)
    if kind === :mean
        m = mean(seg)
        @inbounds for i in eachindex(seg)
            seg[i] = seg[i] - m
        end
    elseif kind === :none
        return seg
    else
        error("Unsupported detrend: $kind")
    end
    return seg
end

function _windowed_rfft!(
    out::Vector{ComplexF64},
    buf::Vector{Float64},
    x::AbstractVector{<:Real},
    start::Integer,
    ws::SpectralWorkspace,
)
    @inbounds for i in 1:ws.nfft
        buf[i] = Float64(x[start + i - 1])
    end
    _detrend_segment!(buf, ws.detrend)
    @inbounds for i in 1:ws.nfft
        buf[i] *= ws.weights[i]
    end
    mul!(out, ws.plan, buf)
    return out
end

function _welch_accumulate!(
    pxx::Vector{Float64},
    x::AbstractVector{<:Real},
    ws::SpectralWorkspace,
)
    n = length(x)
    nfft = ws.nfft
    n < nfft && return 0
    step = nfft - ws.noverlap

    n_freqs = length(ws.freqs)

    k = 0
    start = 1
    while start + nfft - 1 <= n
        Xf = _windowed_rfft!(ws.fft_x, ws.buf_x, x, start, ws)
        @inbounds for i in 1:n_freqs
            pxx[i] += abs2(Xf[i]) * ws.norm_psd
        end

        k += 1
        start += step
    end
    return k
end

function _welch_accumulate!(
    pxx::Vector{Float64},
    x::AbstractVector{<:Real},
    fs::Real,
    nfft::Integer,
    noverlap::Integer,
    window::Symbol,
    detrend::Symbol,
)
    ws = SpectralWorkspace(nfft, fs; noverlap = noverlap, window = window, detrend = detrend)
    return _welch_accumulate!(pxx, x, ws)
end

function _finalize_psd!(pxx::Vector{Float64}, nfft::Integer, k::Integer)
    k == 0 && return pxx
    inv_k = 1.0 / k
    n_freqs = nfft ÷ 2 + 1
    @inbounds for i in 1:n_freqs
        pxx[i] *= inv_k
    end
    @inbounds for i in 2:(n_freqs - 1)
        pxx[i] *= 2.0
    end
    if !iseven(nfft)
        @inbounds pxx[end] *= 2.0
    end
    return pxx
end

function _freq_axis(nfft::Integer, fs::Real)
    n_freqs = nfft ÷ 2 + 1
    return Float64[(i - 1) * fs / nfft for i in 1:n_freqs]
end

function _welch_psd(
    x::AbstractVector{<:Real},
    fs::Real;
    nfft::Integer,
    noverlap::Integer = nfft ÷ 2,
    window::Symbol = :hann,
    detrend::Symbol = :mean,
    workspace = nothing,
)
    ws = workspace === nothing ?
        SpectralWorkspace(nfft, fs; noverlap = noverlap, window = window, detrend = detrend) :
        _validate_workspace(workspace, fs, nfft, noverlap, window, detrend)
    n_freqs = nfft ÷ 2 + 1
    pxx = zeros(Float64, n_freqs)
    k = _welch_accumulate!(pxx, x, ws)
    k == 0 && return (Float64[], Float64[])
    _finalize_psd!(pxx, nfft, k)
    return (ws.freqs, pxx)
end

function _welch_psd_segments(
    segments::AbstractVector,
    fs::Real;
    nfft::Integer,
    noverlap::Integer = nfft ÷ 2,
    window::Symbol = :hann,
    detrend::Symbol = :mean,
    workspace = nothing,
)
    ws = workspace === nothing ?
        SpectralWorkspace(nfft, fs; noverlap = noverlap, window = window, detrend = detrend) :
        _validate_workspace(workspace, fs, nfft, noverlap, window, detrend)
    n_freqs = nfft ÷ 2 + 1
    acc = zeros(Float64, n_freqs)
    total_k = 0
    n_used = 0
    for seg in segments
        length(seg) < nfft && continue
        k = _welch_accumulate!(acc, seg, ws)
        k == 0 && continue
        total_k += k
        n_used += 1
    end
    total_k == 0 && return (Float64[], Float64[], 0)
    _finalize_psd!(acc, nfft, total_k)
    return (ws.freqs, acc, n_used)
end

function _stft_psd(
    x::AbstractVector{<:Real},
    masked::AbstractVector{Bool},
    fs::Real;
    nfft::Integer,
    noverlap::Integer = nfft ÷ 2,
    window::Symbol = :hann,
    detrend::Symbol = :mean,
    workspace = nothing,
)
    n = length(x)
    n == length(masked) || error("x and masked length mismatch")
    step = nfft - noverlap
    step <= 0 && error("noverlap must be < nfft")
    if n < nfft
        return (Float64[], Float64[], zeros(Float64, 0, 0))
    end
    n_freqs = nfft ÷ 2 + 1
    n_times = (n - nfft) ÷ step + 1
    ws = workspace === nothing ?
        SpectralWorkspace(nfft, fs; noverlap = noverlap, window = window, detrend = detrend) :
        _validate_workspace(workspace, fs, nfft, noverlap, window, detrend)
    spec = Matrix{Float64}(undef, n_freqs, n_times)
    times = Vector{Float64}(undef, n_times)
    for t in 1:n_times
        start = (t - 1) * step + 1
        bad = false
        @inbounds for i in 0:(nfft - 1)
            if masked[start + i]
                bad = true
                break
            end
        end
        times[t] = (start - 1 + nfft / 2) / fs
        if bad
            @inbounds for i in 1:n_freqs
                spec[i, t] = NaN
            end
            continue
        end
        Xf = _windowed_rfft!(ws.fft_x, ws.buf_x, x, start, ws)
        @inbounds for i in 1:n_freqs
            p = abs2(Xf[i]) * ws.norm_psd
            if i > 1 && i < n_freqs
                p *= 2.0
            end
            spec[i, t] = p
        end
    end
    if !iseven(nfft)
        @inbounds for t in 1:n_times
            if isfinite(spec[end, t])
                spec[end, t] *= 2.0
            end
        end
    end
    return (ws.freqs, times, spec)
end

function _auto_nfft(window_seconds::Real, fs::Real)
    isfinite(window_seconds) || return 8192
    window_samples = round(Int, window_seconds * fs)
    window_samples <= 0 && return 256
    base = window_samples ÷ 8
    base <= 0 && return 256
    p = 2^floor(Int, log2(base))
    return clamp(p, 256, 8192)
end
