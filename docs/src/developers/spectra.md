# Spectral Views

Pick a spectral view in TKApp's **View** menu to spot noise, mains harmonics or a dead sensor that
a trace hides.

| View | What you get |
|:---|:---|
| `Time \| Spectra` | a PSD panel next to each trace |
| `Spectra` | the same PSD panels at full width, no traces |

**Window** and **Scroll** choose the samples, masks remove theirs, and the spectra update as soon as
you mask or unmask. The line below the plots shows `nfft`, `df = fs/nfft`, Nyquist, segment length
and segment count; the **i** badge explains each term.

## How the spectra are computed

Welch's method on the visible window:

1. Masked samples are **removed, not zero-filled**, which leaves good stretches. A gap between
   joined runs (e.g. high-rate Metronix) also ends a stretch, so no segment spans two recordings.
2. Each stretch is cut into `nfft`-sample segments with 50% overlap. Leftovers shorter than `nfft`
   are dropped.
3. Each segment: remove the mean, apply a Hann window `w`, FFT, scale `|X(f)|²` by
   `1 / (fs · Σ w²)`, double all bins except 0 Hz and Nyquist. The result is a one-sided PSD in
   amplitude²/Hz.
4. Average all segments with equal weight.

Both axes are logarithmic. The 0 Hz bin is hidden, and the y axis shows at most eight decades below
the peak.

## How the transform length is selected

`nfft` = window length in samples ÷ 8, rounded down to a power of two, clamped to `[256, 8192]`
(8192 for `All`). Widen the window for finer frequency resolution, narrow it for finer time
resolution. If no stretch is as long as `nfft`, the panel says so: widen the window or remove a mask.

## Performance

Workspaces (taper, FFT plan, buffers) are cached by `(nfft, fs, noverlap, window, detrend)`, and a
short timer delays recomputation while you scroll.

!!! note "Internal API"
    `Timekeepers._welch_psd`, `Timekeepers._welch_psd_segments` and
    `Timekeepers.SpectralWorkspace` are internal and may change. For your own spectral analysis, get
    clean segments and use a dedicated package:

    ```julia
    segments = good_segments(ta, mask; min_samples = 4096)
    ```
