# Spectral Views

The **View** menu in [TKApp](tkapp.md) adds a diagnostic panel to the time
series. It is there to answer the question you actually have while masking: *is
this stretch of record usable?* Broadband noise, a mains harmonic that comes and
goes, a sensor that stopped responding above some frequency — all of these are
far easier to see in a spectrum than in a trace.

## `Time | Spectra` — Welch PSD

Each channel gets a power spectral density panel beside its trace, estimated by
Welch's method: the visible window is split into overlapping segments, each is
mean-detrended, tapered with a Hann window, transformed, and the resulting
periodograms are averaged.

Crucially, the PSD is computed over the *unmasked stretches only*. Masked
intervals are excluded rather than zero-filled, so cutting a spike immediately
cleans the spectrum instead of replacing it with the spectrum of a step edge.
A gap between two runs joined end to end — high-rate Metronix runs, see
[Metronix Sites](metronix.md) — ends a stretch the same way, so no segment
joins two recordings; a window covering several runs averages the segments of
each. Both axes are logarithmic; the y axis is amplitude²/Hz in the channel's
own units.

The header line under the plots reports the configuration in use — transform
length, frequency resolution `df = fs/nfft`, Nyquist frequency, and the
segment duration — along with how many segments were averaged, from how many
unmasked stretches and runs. The **i** badge at its left swaps that line for a
plain-language gloss of each term; the full method is below, under
[How the spectra are computed](#How-the-spectra-are-computed).

## `Spectra` — the same panels, full width

The same estimate with the traces switched off, so the PSD panels take the whole
window. Everything else is unchanged: **Window** and **Scroll** still choose
which samples are estimated, masks still exclude their samples, and the cursor
and pin behave the same. It suits the stretches of work where the spectrum is
what you are reading — comparing channels, or chasing a harmonic down to where
it disappears into the noise floor.

## How the spectra are computed

Each channel's panel is one power spectral density, estimated by Welch's
method from the samples on screen. The steps:

1. **Samples.** Only the samples inside the visible **Window** are used. Masked
   samples are left out — not replaced by zeros, because a block of zeros adds
   step edges whose power spreads across every frequency. What remains is a set
   of *unmasked stretches*: each runs from one masked sample (or edge of the
   window) to the next. A gap between two runs joined end to end — high-rate
   Metronix runs, see [Metronix Sites](metronix.md) — ends a stretch just as a
   mask does, so no segment joins the end of one recording to the start of the
   next, and every run contributes its own segments.
2. **Segments.** Each stretch is cut into segments of `nfft` samples, each
   starting half a segment after the previous one (50% overlap). A stretch
   shorter than `nfft` contributes nothing, and the samples at the end of a
   stretch that do not fill a whole segment are not used. How `nfft` is
   chosen is described in the next section.
3. **Each segment.** Its mean is removed, it is tapered with a Hann window
   `w`, and Fourier transformed. `|X(f)|²` is scaled by `1 / (fs · Σ w²)`, and
   every bin except 0 Hz and Nyquist is doubled to fold in the negative
   frequencies — a one-sided PSD in amplitude²/Hz, in the channel's own units.
4. **Average.** All segments, from every stretch and every run in the window,
   are averaged with equal weight. More segments give a smoother, steadier
   estimate; a narrow window with few segments gives a noisier one.

Reading the result: bins are `df = fs/nfft` apart, from `df` up to the Nyquist
frequency `fs/2`; the 0 Hz bin, which holds only the mean, is not drawn. Both
axes are logarithmic, and the y axis spans at most eight decades below the
peak, so an instrument's anti-alias roll-off in the last bins before Nyquist
does not flatten the rest of the band. Masking or unmasking recomputes the
spectra at once, so the effect of a cut shows immediately.

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

If the visible window is shorter than `nfft`, or every unmasked stretch in it
is, the panel reports that instead of drawing a misleading spectrum. Widen the
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
