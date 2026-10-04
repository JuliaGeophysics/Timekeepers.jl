# Spectral Views

The **View** menu in [TKApp](tkapp.md) adds a diagnostic panel to the time
series. The panel helps you during the mask work: *can you use this part of
the record?* Some problems are much easier to see in a spectrum than in a
trace:

- broadband noise;
- a mains harmonic that starts and stops;
- a sensor that does not respond above a frequency.

## `Time | Spectra`: Welch PSD

Each channel gets a power spectral density (PSD) panel next to its trace. The
app uses the Welch method:

1. It cuts the visible window into segments that overlap.
2. It removes the mean of each segment.
3. It applies a Hann window to each segment.
4. It does the Fourier transform of each segment.
5. It calculates the average of the periodograms.

The PSD uses *only the stretches without a mask*. The app removes the masked
intervals. It does not fill them with zeros. Thus, when you cut a spike, the
spectrum becomes clean immediately. The spectrum of a step edge does not
replace it.

Some runs are joined end to end, for example high-rate Metronix runs (refer
to [Metronix Sites](metronix.md)). A gap between two such runs also ends a
stretch. Thus, no segment joins two recordings. For a window that covers more
than one run, the app calculates the average of the segments of each run.

The two axes are logarithmic. The y axis is amplitude²/Hz in the units of the
channel.

The header line below the plots gives the configuration:

- the transform length;
- the frequency resolution `df = fs/nfft`;
- the Nyquist frequency;
- the segment duration;
- the number of segments in the average, and the number of stretches and
  runs that they come from.

The **i** badge at the left of the line replaces it with a plain explanation
of each term. For the full method, refer to
[How the spectra are computed](#How-the-spectra-are-computed).

## `Spectra`: the same panels at full width

This view gives the same estimate without the traces. Thus, the PSD panels use
the full window. All the other functions stay the same:

- **Window** and **Scroll** select the samples for the estimate.
- The masks remove their samples.
- The cursor and the pin operate in the same way.

Use this view when you read mainly the spectrum. For example, compare
channels, or follow a harmonic down to the noise floor.

## How the spectra are computed

The panel of each channel shows one power spectral density. The app
calculates it with the Welch method from the samples on the screen. Do these
steps:

1. **Samples.** The app uses only the samples in the visible **Window**. It
   removes the masked samples. It does not replace them with zeros, because a
   block of zeros adds step edges, and their power goes into all frequencies.
   The result is a set of *stretches without a mask*. Each stretch goes from
   one masked sample (or edge of the window) to the next. Some runs are joined
   end to end, for example high-rate Metronix runs (refer to
   [Metronix Sites](metronix.md)). A gap between two such runs ends a stretch,
   as a mask does. Thus, no segment joins the end of one recording to the
   start of the next. Each run gives its own segments.
2. **Segments.** The app cuts each stretch into segments of `nfft` samples.
   Each segment starts half a segment after the previous one (50% overlap). A
   stretch shorter than `nfft` gives no segment. The app does not use the
   samples at the end of a stretch that do not fill a full segment. The next
   section tells how the app selects `nfft`.
3. **Each segment.** The app removes the mean of the segment. It applies a
   Hann window `w`. It then does the Fourier transform. It multiplies `|X(f)|²`
   by `1 / (fs · Σ w²)`. It doubles each bin, but not 0 Hz and Nyquist, to add
   the negative frequencies. The result is a one-sided PSD in amplitude²/Hz,
   in the units of the channel.
4. **Average.** The app calculates the average of all the segments, from each
   stretch and each run in the window. Each segment has the same weight. More
   segments give a smoother and more stable estimate. A narrow window with few
   segments gives more noise.

How to read the result:

- The bins are `df = fs/nfft` apart, from `df` up to the Nyquist frequency
  `fs/2`.
- The app does not show the 0 Hz bin, because it holds only the mean.
- The two axes are logarithmic.
- The y axis shows at most eight decades below the peak. Thus, the anti-alias
  roll-off of an instrument in the last bins before Nyquist does not make the
  rest of the band flat.
- When you mask or unmask, the app calculates the spectra again immediately.
  Thus, you see the effect of a cut immediately.

## How the transform length is selected

The app calculates `nfft` from the visible window. The value is not fixed.
Thus, the panel stays useful when you zoom:

1. The app divides the window length in samples by 8. This gives
   approximately eight segments across the view.
2. It rounds the result down to a power of two.
3. It limits the result to `[256, 8192]`.

When **Window** is `All`, `nfft` is 8192. The overlap is always `nfft ÷ 2`.

Thus, a wider window gives better frequency resolution. A narrower window
gives better time resolution. You control this balance with the same control
that you use to scroll.

Sometimes the visible window, or each stretch without a mask in it, is shorter
than `nfft`. Then the panel tells you this and does not show an incorrect
spectrum. Make the window wider, or remove a mask.

## Performance

The spectral estimate uses a workspace that the app uses again. The workspace
holds the taper, the FFT plan and its scratch buffers. Its key is
`(nfft, fs, noverlap, window, detrend)`. The app keeps the workspaces and uses
them again when you scroll. Thus, more estimates at one configuration use no
new memory, except for the output arrays. A timer also delays each new
calculation. Thus, when you move the scroll slider, the app does not do one
FFT pass for each frame.

!!! note "Internal API"
    The estimators (`Timekeepers._welch_psd`,
    `Timekeepers._welch_psd_segments` and `Timekeepers.SpectralWorkspace`)
    are internal. Semantic versioning does not apply to them. For spectral
    analysis in your own code, get clean segments from Timekeepers with
    [`good_segments`](@ref). Then use a dedicated package such as
    [DSP.jl](https://github.com/JuliaDSP/DSP.jl):

    ```julia
    segments = good_segments(ta, mask; min_samples = 4096)
    ```
