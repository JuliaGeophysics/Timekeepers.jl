# Spectral Views

The **View** menu in [TKApp](tkapp.md) adds a diagnostic panel next to the
time series. It is there to answer the question you actually have while
masking: *is this stretch of record usable?* Broadband noise, a mains harmonic
that comes and goes, a sensor that stopped responding above some frequency —
all of these are far easier to see in a spectrum than in a trace.

## `Time | Spectra` — Welch PSD

Each channel gets a power spectral density panel beside its trace, estimated by
Welch's method: the visible window is split into overlapping segments, each is
mean-detrended, tapered with a Hann window, transformed, and the resulting
periodograms are averaged.

Crucially, the PSD is computed over the *good segments only*. Masked intervals
are excluded rather than zero-filled, so cutting a spike immediately cleans the
spectrum instead of replacing it with the spectrum of a step edge. Both axes are
logarithmic; the y axis is amplitude²/Hz in the channel's own units.

The header line under the plots reports the configuration in use — transform
length, frequency resolution `df = fs/nfft`, Nyquist frequency, and the
segment duration — along with how many good segments were averaged.

## How the transform length is chosen

`nfft` is derived from the visible window rather than fixed, so the panel stays
informative as you zoom:

- the window length in samples is divided by 8, giving roughly eight segments
  across the view;
- that is rounded down to a power of two;
- the result is clamped to `[256, 8192]`.

With **Window** set to `All`, `nfft` is 8192. Overlap is always `nfft ÷ 2`.

Widening the window therefore buys frequency resolution, and narrowing it buys
time resolution — the usual trade, driven by the same control you already use
to scroll.

If the visible window is shorter than `nfft`, or every good segment in it is,
the panel reports that instead of drawing a misleading spectrum. Widen the
window or unmask something.

## Performance

Spectral estimation runs through a reusable workspace holding the taper, the
FFT plan and its scratch buffers, keyed by
`(nfft, fs, noverlap, window, detrend)`. Workspaces are cached on the app and
reused as you scroll, so repeated estimates at one configuration allocate
nothing beyond the output arrays. Recomputation is also debounced behind a
timer, so dragging the scroll slider does not queue one FFT pass per frame.

!!! note "Internal API"
    The estimators themselves (`Timekeepers._welch_psd`,
    `Timekeepers._welch_psd_segments` and `Timekeepers.SpectralWorkspace`) are
    internal and not covered by semantic versioning. For spectral analysis in
    your own code, take clean segments out of Timekeepers with
    [`good_segments`](@ref) and use a dedicated package such as
    [DSP.jl](https://github.com/JuliaDSP/DSP.jl):

    ```julia
    segments = good_segments(ta, mask; min_samples = 4096)
    ```
