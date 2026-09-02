# TKApp Explorer

TKApp is the interactive half of Timekeepers: a native GLMakie window for
scrolling through a long record, marking bad intervals by eye, and writing the
result back out in the format it came from.

![TKApp time series view](assets/ts.png)

## Opening the window

[`run_tkapp`](@ref) builds the app and blocks until the window closes, then
returns the [`TKApp`](@ref) so the mask you made survives the session:

```julia
using Timekeepers

run_tkapp()                                          # start empty, load from the toolbar
run_tkapp("data/LEMI090.txt")                        # open one file
run_tkapp("data/RK137")                              # open a whole site directory
app = run_tkapp(load_lemi424("data/LEMI090.txt"))    # open an in-memory TimeArray
```

To build the app without blocking — for scripting, or to inspect it before
showing it — construct a [`TKApp`](@ref) and `display` it:

```julia
app = TKApp("data/LEMI090.txt")
display(app)
```

A bundled launcher is available from a repository clone:

```bash
julia --project=. examples/tkapp.jl
```

!!! warning "Needs a real display"
    GLMakie requires a desktop session with OpenGL 3.3 or newer. See
    [Getting Started](getting_started.md#Checking-your-OpenGL-setup) for a
    one-line smoke test. Everything outside the app works headless.

## The toolbar

| Control | What it does |
|:---|:---|
| **Load Run…** | Open one run — a `.txt`, `.dat`, `.lem` or `.xyz` file, or a Metronix `meas_*` directory (see below) |
| **Load Site…** | Open a directory: every run in it is read, ordered by start time and concatenated, with `NaN` filling the gaps between runs |
| **Mask** | Mark the current selection bad |
| **Unmask** | Mark the current selection good again |
| **Clear** | Drop every mask on the record |
| **Write** | Export the cleaned data and the mask (see below) |
| **View** | `Time`, `Spectra`, or `Time \| Spectra` |
| **Window** | Visible span: 30 seconds through 7 days, or `All` |

`Spectra` drops the traces and gives the PSD panels the full width, for the
stretches of work where the spectrum is what you are reading. The **Window**
menu and **Scroll** still choose which samples are estimated, so the spectra
follow the window exactly as they do beside a trace.

## Opening a Metronix run

A Metronix run is a whole `meas_*` directory — several `.ats` channel files plus
the XML header — and file dialogs select files. So **Load Run…** takes you into
the `meas_*` directory and you pick any `.ats` inside it; the run around that
file loads, every channel and the header together. **Load Site…** accepts a
`meas_*` directory too, for anyone who prefers the folder picker.

Point **Load Site…** at the site above those runs and it reads the whole site.
A record carries one sampling rate, so when the site holds more than one the app
asks which to load before it starts, and a **Write** afterwards covers only the
runs at that rate. See [Metronix Sites](metronix.md).

Below the plots, **Scroll** moves the visible window through the record — drag
the slider, or step a window at a time with **&lt;** and **&gt;**.

The status line at the bottom reports what was loaded, what was written, and
any error, so a failed load or write does not disappear into the REPL.

## Selecting and masking

Within a time-series panel:

- **left-drag** — select a time interval; the span highlights across every
  channel at once
- **right-click** — mask the current selection, the same as pressing **Mask**
- **right-drag** — pan
- **scroll** — zoom the y axis

Masked spans render as shaded bands, and the samples inside them are drawn in
grey rather than the channel colour, so you can always see what you cut and
undo it with **Unmask**.

Selecting across the full window and pressing **Mask** is the fast way to
discard a whole bad day; narrowing the **Window** to a minute is the way to cut
a single spike cleanly.

## Writing results

**Write** picks its behaviour from the source format.

For **LEMI-424, GEOMAG and `.xyz`** records it writes two files next to the
source, named after it:

- `<name>_clean.<ext>` — the cleaned record in the original format
- `<name>_mask.csv` — the masked intervals, replayable with [`read_mask`](@ref)

Loading a *site directory* writes `<site>_combined_clean.<ext>` and
`<site>_combined_mask.csv` inside that directory.

For **Metronix** sites, blanking is not an option — the format has no `NaN`, and
downstream tools expect continuous runs. Instead the app calls
[`write_metronix_site_masked`](@ref), which amputates the masked intervals and
writes the surviving stretches as separate `meas_*` directories under
`<site>.W`, with a `README.md` recording each write session. A progress console
opens while this runs. See [Metronix Sites](metronix.md).

## Continuing in code

Every masking accessor takes the app directly, so a session done by hand
continues in the REPL without unpacking anything:

```julia
app = run_tkapp("data/LEMI090.txt")   # mask a few intervals, then close the window

app.data                                    # the loaded TimeArray
app.mask                                    # the TimekeeperMask
masked_samples(app.mask)

cleaned  = cleaned_timearray(app)           # NaN in masked rows
segments = good_segments(app; min_samples = 256)
weights  = sample_weights(app)

write_cleaned("out/clean.csv", app)
write_mask("out/mask.csv", app)
```

See [Masking & Cleaning](masking.md) for what each of those produces.
