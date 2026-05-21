using FFTW

function _hann_window(n::Integer)
    n <= 1 && return ones(Float64, max(n, 0))
    return Float64[0.5 * (1 - cos(2π * (k - 1) / (n - 1))) for k in 1:n]
end

function _make_window(kind::Symbol, n::Integer)
    (kind === :hann || kind === :hanning) && return _hann_window(n)
    kind === :rect && return ones(Float64, n)
    error("Unsupported window: $kind")
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

function _welch_accumulate!(
    pxx::Vector{Float64},
    pyy::Vector{Float64},
    pxy::Vector{ComplexF64},
    x::AbstractVector{<:Real},
    y::Union{Nothing, AbstractVector{<:Real}},
    fs::Real,
    nfft::Integer,
    noverlap::Integer,
    window::Symbol,
    detrend::Symbol,
)
    n = length(x)
    n < nfft && return 0
    step = nfft - noverlap
    step <= 0 && error("noverlap must be < nfft")
    use_y = y !== nothing
    use_y && (length(y) == n || error("x and y must have the same length"))

    w = _make_window(window, nfft)
    sum_w_sq = sum(abs2, w)
    norm_psd = 1.0 / (fs * sum_w_sq)

    n_freqs = nfft ÷ 2 + 1
    buf_x = Vector{Float64}(undef, nfft)
    buf_y = use_y ? Vector{Float64}(undef, nfft) : Float64[]

    k = 0
    start = 1
    while start + nfft - 1 <= n
        @inbounds for i in 1:nfft
            buf_x[i] = Float64(x[start + i - 1])
        end
        _detrend_segment!(buf_x, detrend)
        @inbounds for i in 1:nfft
            buf_x[i] *= w[i]
        end
        Xf = rfft(buf_x)

        if use_y
            @inbounds for i in 1:nfft
                buf_y[i] = Float64(y[start + i - 1])
            end
            _detrend_segment!(buf_y, detrend)
            @inbounds for i in 1:nfft
                buf_y[i] *= w[i]
            end
            Yf = rfft(buf_y)
            @inbounds for i in 1:n_freqs
                pxx[i] += abs2(Xf[i]) * norm_psd
                pyy[i] += abs2(Yf[i]) * norm_psd
                pxy[i] += (Xf[i] * conj(Yf[i])) * norm_psd
            end
        else
            @inbounds for i in 1:n_freqs
                pxx[i] += abs2(Xf[i]) * norm_psd
            end
        end

        k += 1
        start += step
    end
    return k
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

function _finalize_cross!(pxy::Vector{ComplexF64}, nfft::Integer, k::Integer)
    k == 0 && return pxy
    inv_k = 1.0 / k
    n_freqs = nfft ÷ 2 + 1
    @inbounds for i in 1:n_freqs
        pxy[i] *= inv_k
    end
    @inbounds for i in 2:(n_freqs - 1)
        pxy[i] *= 2.0
    end
    if !iseven(nfft)
        @inbounds pxy[end] *= 2.0
    end
    return pxy
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
)
    n_freqs = nfft ÷ 2 + 1
    pxx = zeros(Float64, n_freqs)
    dummy_y = zeros(Float64, 0)
    dummy_xy = zeros(ComplexF64, 0)
    k = _welch_accumulate!(pxx, dummy_y, dummy_xy, x, nothing, fs,
        nfft, noverlap, window, detrend)
    k == 0 && return (Float64[], Float64[])
    _finalize_psd!(pxx, nfft, k)
    return (_freq_axis(nfft, fs), pxx)
end

function _welch_psd_segments(
    segments::AbstractVector,
    fs::Real;
    nfft::Integer,
    noverlap::Integer = nfft ÷ 2,
    window::Symbol = :hann,
    detrend::Symbol = :mean,
)
    n_freqs = nfft ÷ 2 + 1
    acc = zeros(Float64, n_freqs)
    dummy_y = zeros(Float64, 0)
    dummy_xy = zeros(ComplexF64, 0)
    total_k = 0
    n_used = 0
    for seg in segments
        length(seg) < nfft && continue
        k = _welch_accumulate!(acc, dummy_y, dummy_xy, seg, nothing, fs,
            nfft, noverlap, window, detrend)
        k == 0 && continue
        total_k += k
        n_used += 1
    end
    total_k == 0 && return (Float64[], Float64[], 0)
    _finalize_psd!(acc, nfft, total_k)
    return (_freq_axis(nfft, fs), acc, n_used)
end

function _welch_coherence(
    x::AbstractVector{<:Real},
    y::AbstractVector{<:Real},
    fs::Real;
    nfft::Integer,
    noverlap::Integer = nfft ÷ 2,
    window::Symbol = :hann,
    detrend::Symbol = :mean,
)
    n_freqs = nfft ÷ 2 + 1
    pxx = zeros(Float64, n_freqs)
    pyy = zeros(Float64, n_freqs)
    pxy = zeros(ComplexF64, n_freqs)
    k = _welch_accumulate!(pxx, pyy, pxy, x, y, fs,
        nfft, noverlap, window, detrend)
    k == 0 && return (Float64[], Float64[])
    inv_k = 1.0 / k
    @inbounds for i in 1:n_freqs
        pxx[i] *= inv_k
        pyy[i] *= inv_k
        pxy[i] *= inv_k
    end
    gamma_sq = Vector{Float64}(undef, n_freqs)
    @inbounds for i in 1:n_freqs
        denom = pxx[i] * pyy[i]
        gamma_sq[i] = denom > 0 ? min(abs2(pxy[i]) / denom, 1.0) : 0.0
    end
    return (_freq_axis(nfft, fs), gamma_sq)
end

function _welch_coherence_segments(
    segs_x::AbstractVector,
    segs_y::AbstractVector,
    fs::Real;
    nfft::Integer,
    noverlap::Integer = nfft ÷ 2,
    window::Symbol = :hann,
    detrend::Symbol = :mean,
)
    length(segs_x) == length(segs_y) || error("segs_x and segs_y must align")
    n_freqs = nfft ÷ 2 + 1
    acc_xx = zeros(Float64, n_freqs)
    acc_yy = zeros(Float64, n_freqs)
    acc_xy = zeros(ComplexF64, n_freqs)
    total_k = 0
    for (xs, ys) in zip(segs_x, segs_y)
        length(xs) == length(ys) || error("segment x/y lengths differ")
        length(xs) < nfft && continue
        k = _welch_accumulate!(acc_xx, acc_yy, acc_xy, xs, ys, fs,
            nfft, noverlap, window, detrend)
        total_k += k
    end
    total_k == 0 && return (Float64[], Float64[])
    inv_k = 1.0 / total_k
    @inbounds for i in 1:n_freqs
        acc_xx[i] *= inv_k
        acc_yy[i] *= inv_k
        acc_xy[i] *= inv_k
    end
    gamma_sq = Vector{Float64}(undef, n_freqs)
    @inbounds for i in 1:n_freqs
        denom = acc_xx[i] * acc_yy[i]
        gamma_sq[i] = denom > 0 ? min(abs2(acc_xy[i]) / denom, 1.0) : 0.0
    end
    return (_freq_axis(nfft, fs), gamma_sq)
end

function _stft_psd(
    x::AbstractVector{<:Real},
    masked::AbstractVector{Bool},
    fs::Real;
    nfft::Integer,
    noverlap::Integer = nfft ÷ 2,
    window::Symbol = :hann,
    detrend::Symbol = :mean,
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
    w = _make_window(window, nfft)
    sum_w_sq = sum(abs2, w)
    norm_psd = 1.0 / (fs * sum_w_sq)
    spec = Matrix{Float64}(undef, n_freqs, n_times)
    times = Vector{Float64}(undef, n_times)
    buf = Vector{Float64}(undef, nfft)
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
        @inbounds for i in 1:nfft
            buf[i] = Float64(x[start + i - 1])
        end
        _detrend_segment!(buf, detrend)
        @inbounds for i in 1:nfft
            buf[i] *= w[i]
        end
        Xf = rfft(buf)
        @inbounds for i in 1:n_freqs
            p = abs2(Xf[i]) * norm_psd
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
    freqs = Float64[(i - 1) * fs / nfft for i in 1:n_freqs]
    return (freqs, times, spec)
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
