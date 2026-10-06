# Precompile.jl - precompilation workload.
# Author: @pankajkmishra
#
# This file runs the code that a user needs first, on a small synthetic series:
# masks, clean data, segments, Welch estimation, an app with the spectra view,
# the survey dashboard and a transfer function with its window. Thus, the
# package image contains these compiled methods, and the first use is fast.
# This file runs only when the package builds. Nothing here runs when you use
# the package

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

        fbins = collect(range(0.01, 0.5; length = 64))
        pbins = Float64[1 / f for f in fbins]
        _nearest_index(fbins, 0.1)
        _peak_window(fbins, pbins, 0.1, 0.09, 0.11)
        _psd_cursor_text(0.1, 1.0e-3)
        _pin_label_text(0.1)

        app = TKApp(ta; size = (900, 620))
        app.view_mode[] = :time_spectra
        _build_axes!(app, app.data)
        _recompute_spectra!(app)
        _set_spectral_pin!(app, 0.1)
        _clear_spectral_pin!(app)
        # A survey of two sites for the dashboard. It is made from records and
        # not from a scan. Thus, it needs no files
        sruns(t) = [SurveyRun("", 8.0, t, t + Hour(2), 8 * 7200, [:e1, :e2, :bx, :by, :bz])]
        survey = Survey("", [SurveySite("a", "a", :metronix, 48.60, 7.60, 0.0, sruns(t0)),
                             SurveySite("b", "b", :metronix, 48.61, 7.62, 0.0, sruns(t0 + Hour(1)))])
        dash = TKDash(survey; size = (900, 620))
        _focus!(dash, 2)
        _toggle!(dash, 1)
        reference_plan(dash)
        # A transfer function of a synthetic site with a remote site, and the
        # TKProc window that shows it
        nt = 8192
        hx = cumsum(sin.(0.37 .* (1:nt)) .+ cos.(0.11 .* (1:nt)))
        hy = cumsum(cos.(0.29 .* (1:nt)) .+ sin.(0.05 .* (1:nt)))
        mk(site, cols) = TimekeeperRun(site, "synthetic", :synthetic,
            Dict(c => TimekeeperChannel(c, v, 1.0, t0, "", "", Dict{String, Any}()) for (c, v) in cols),
            Dict{Symbol, Any}(:latitude => 48.0, :longitude => 7.0))
        loc = mk("loc", Dict(:e1 => 2 .* hy, :e2 => -hx, :bz => 0.1 .* hx, :bx => hx, :by => hy))
        tf = estimate_tf(loc; remote = mk("rem", Dict(:bx => hx, :by => hy)), window = 128)
        apparent_resistivity(tf)
        proc = TKProc(; size = (900, 620))
        proc.results["loc"] = tf
        plot_tf(tf)
        # The process_interaction methods need a live viewport and a real mouse
        # event. Thus, they compile at the first hover
    end
end
