# Precompile.jl - precompilation workload.
# Author: @pankajkmishra
#
# Exercises the paths a user hits first - masking, cleaning, segmenting, Welch
# and STFT estimation, and building an app with each spectral view - on a small
# synthetic series, so the methods involved are compiled into the package image
# rather than on first use. Runs at build time only; nothing here is called at
# run time.

@setup_workload begin
    t0 = DateTime(2020, 1, 1)
    n = 512
    comps = [:bx, :by, :bz, :e1, :e2]
    times = [t0 + Second(i - 1) for i in 1:n]
    vals = Matrix{Float64}(undef, n, length(comps))
    @inbounds for i in 1:n
        phase = 2π * (i - 1) / 64
        vals[i, 1] = sin(phase)
        vals[i, 2] = cos(phase)
        vals[i, 3] = sin(phase / 2)
        vals[i, 4] = 0.25 * cos(phase / 3)
        vals[i, 5] = 0.25 * sin(phase / 3)
    end
    metadata = Dict{Symbol, Any}(
        :site => "precompile",
        :instrument => "synthetic",
        :source_format => :lemi424,
        :sample_rate => 1.0,
        :start_time => first(times),
        :units => Dict(c => component_units(c) for c in comps),
    )
    ta = TimeArray(times, vals, comps, metadata)

    @compile_workload begin
        mask = TimekeeperMask(ta)
        mask_interval!(mask, times[101], times[140])
        cleaned_timearray(ta, mask)
        good_segments(ta, mask; min_samples = 64)
        sample_weights(mask)

        ws = SpectralWorkspace(256, 1.0)
        x = view(vals, :, 1)
        _welch_psd(x, 1.0; nfft = 256, workspace = ws)
        _stft_psd(x, mask.masked, 1.0; nfft = 256, workspace = ws)

        app = TKApp(ta; size = (900, 620))
        app.view_mode[] = :time_spectra
        _build_axes!(app, app.data, app.axes)
        _recompute_spectra!(app)
        app.view_mode[] = :time_spectrogram
        _build_axes!(app, app.data, app.axes)
        _recompute_spectra!(app)
    end
end
