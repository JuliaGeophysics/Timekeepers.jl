# TKApp Explorer

TKApp is the interactive part of Timekeepers. It is a native GLMakie window.
Use it to:

- scroll through a long record;
- find bad intervals by eye and mark them;
- write the result in the format of the source.

## Open the window

[`run_tkapp`](@ref) makes the app and waits until you close the window. It
then returns the [`TKApp`](@ref). Thus, you keep the mask that you made:

```julia
using Timekeepers

run_tkapp()                                          # start empty, load from the toolbar
run_tkapp("data/LEMI090.txt")                        # open one file
run_tkapp("data/RK137")                              # open a whole site directory
app = run_tkapp(load_lemi424("data/LEMI090.txt"))    # open an in-memory TimeArray
```

You can make the app and not wait, for example in a script, or to examine the
app before it shows. Make a [`TKApp`](@ref) and `display` it:

```julia
app = TKApp("data/LEMI090.txt")
display(app)
```

A clone of the repository has a launcher:

```bash
julia --project=. examples/tkapp.jl
```

!!! warning "Needs a real display"
    GLMakie needs a desktop session with OpenGL 3.3 or newer. Refer to
    [Getting Started](getting_started.md#Checking-your-OpenGL-setup) for a
    smoke test of one line. All the functions outside the app operate without
    a display.

## The toolbar

| Control | What it does |
|:---|:---|
| **Rate** (left) | The sampling rate on the screen. For a Metronix site with more than one rate, the menu gives each rate |
| **Load Run…** | Opens one run: a `.txt`, `.dat`, `.lem` or `.xyz` file, or a Metronix run by its `.xml` (refer to the subsequent section) |
| **Load Site…** | Opens a directory. The app reads all the runs in it, puts them in order of start time and joins them. It puts `NaN` in the gaps between runs |
| **Mask** | Marks the current selection as bad |
| **Unmask** | Marks the current selection as good again |
| **Clear** | Removes all the masks from the record |
| **Write** | Writes the clean data and the mask (refer to the subsequent section) |
| **View** | `Time`, `Spectra`, or `Time \| Spectra` |
| **Window** | The visible span. Type a whole number and push Enter. Then select its unit: `seconds`, `minutes`, `hours` or `days`. Or select `All` for the full record. A new app shows `All` |

`Spectra` removes the traces and gives the PSD panels the full width. Use it
when you read mainly the spectrum. The **Window** span and **Scroll** still
select the samples for the estimate. Thus, the spectra follow the window in
the same way as next to a trace. When you change the view, the new layout
shows immediately. The app makes new spectrum panels a short time later.

The checkbox next to each channel turns the channel off and on. A channel that
is off keeps its position. Its panels become grey, and its traces and spectrum
are not shown. Thus, a record with only electric channels, or without a
vertical field, is easy to read, and the layout does not move. This switch
changes only the view. The masks and **Write** still apply to all the
channels.

Some channels have no data from the logger. Such a channel is all `NaN`, or
has one constant value, for example the zeros that the logger writes for an
input that is not connected. Such a channel starts in the off state, and the
app does not mask it automatically.

## Open a Metronix run

A Metronix run is a set of `.ats` channel files and an XML header. One `meas_*`
directory can hold more than one run at more than one rate. Thus, **Load
Run…** takes you into the `meas_*` directory. Select the `.xml` of the run, or
one of its `.ats` files. The app loads that run, with all its channels and the
header.

**Load Site…** accepts the site above the `meas_*` directories, or one
`meas_*` directory. It finds all the runs and writes nothing. The app holds
one sampling rate at a time. Thus, for a site with more than one rate, it asks
which rate to import. Later, use the rate menu at the left of the toolbar to
change to a different rate. Each rate keeps its own masked intervals.
**Write** includes each rate that you masked. Refer to
[Metronix Sites](metronix.md).

Below the plots, **Scroll** moves the visible window through the record.
Move the slider, or push **&lt;** and **&gt;** to move one window at a time.

The status line at the bottom tells you what the app loaded and wrote, and
each error. Thus, you see each failed load or write, and the error does not
go only to the REPL.

## Select and mask

In a time-series panel:

- **left-drag**: select a time interval. The span becomes bright across all
  the channels.
- **right-click**: mask the current selection. This is the same as **Mask**.
- **right-drag**: pan.
- **scroll**: zoom the y axis.

Masked spans show as shaded bands. The app draws the samples in them in grey,
not in the color of the channel. Thus, you can always see what you cut, and
undo it with **Unmask**.

To remove a full bad day quickly, select the full window and push **Mask**. To
cut one spike cleanly, make the **Window** one minute long.

## Write the results

**Write** selects its behavior from the source format.

For **LEMI-424, GEOMAG and `.xyz`** records, it writes two files next to the
source. The names come from the source:

- `<name>_clean.<ext>`: the clean record in the original format.
- `<name>_mask.csv`: the masked intervals. Use [`read_mask`](@ref) to apply
  them again.

If you loaded a *site directory*, **Write** writes
`<site>_combined_clean.<ext>` and `<site>_combined_mask.csv` in that
directory.

For **Metronix** sites, blanks are not an option. The format has no `NaN`,
and the next tools need continuous runs. Thus, the app calls
[`write_metronix_site_masked`](@ref). This function cuts the masked intervals
out of the rate directories of the site, `<site>.128`, `<site>.4096`, ..., next
to it. The site itself never changes.

- The runs without a cut stay as they are.
- Each part after a cut becomes a new run number in its `meas_*` directory.
  It has its own XML, edited from the XML of the run, and the `.kml` next to
  it.
- The function ignores parts that are too short to keep as a run, and gives a
  warning.
- A `mask.csv` in each rate directory gives each part that was cut.

The app then shows the record after the cut, with its masks applied. Refer to
[Metronix Sites](metronix.md).

## Continue in code

Each mask function accepts the app directly. Thus, you can make masks by hand
and continue in the REPL. You do not have to unpack the app:

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

Refer to [Masking & Cleaning](masking.md) for the output of each function.
