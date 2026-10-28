# TKApp in Detail

```julia
using Timekeepers

run_tkapp()                                          # start empty, load from the toolbar
run_tkapp("data/LEMI090.txt")                        # one file
run_tkapp("data/RK137")                              # a whole site directory
app = run_tkapp(load_lemi424("data/LEMI090.txt"))    # an in-memory TimeArray

app = TKApp("data/LEMI090.txt"); display(app)        # don't wait for the window to close
```

Or from a clone: `julia --project=. examples/tkapp.jl`.

!!! warning "Needs a display"
    GLMakie needs OpenGL 3.3 or newer ([check it](../getting_started.md)).
    Everything outside the windows works without one.

## Toolbar

| Control | What it does |
|:---|:---|
| **Rate** (left) | Sampling rate on screen; a Metronix site with several rates lists them all |
| **Load Run…** | One `.txt`, `.dat`, `.lem` or `.xyz` file, or a Metronix run (pick its `.xml` or an `.ats` inside `meas_*`) |
| **Load Site…** | A directory: all runs joined in time order, `NaN` in the gaps. A Metronix site loads one rate at a time |
| **Mask** / **Unmask** | Mark the selection bad / good |
| **Clear** | Remove every mask |
| **Write** | Save clean data and the mask (below) |
| **View** | `Time`, `Spectra` or `Time \| Spectra` |
| **Window** | Visible span: a number (press Enter) and a unit, or `All` |

## Mouse

| | |
|:---|:---|
| left-drag | select an interval |
| right-click | mask the selection |
| right-drag | pan |
| scroll | zoom the y axis |

Tips:

- Cut a whole bad day: select the full window, press **Mask**. Cut a spike: set **Window** to one
  minute first.
- The checkbox next to a channel only hides it; masks and **Write** still apply. Empty channels (all
  `NaN` or constant) start hidden.
- Each Metronix rate keeps its own masks, and **Write** cuts every rate you masked.

## What Write saves

| Source | Files |
|:---|:---|
| LEMI-424, GEOMAG, `.xyz` | `<name>_clean.<ext>` (masked rows `NaN`) and `<name>_mask.csv`, next to the source |
| A site directory | `<site>_combined_clean.<ext>` and `<site>_combined_mask.csv` in it |
| Metronix | cuts in `<site>.<rate>` copies via [`write_metronix_site_masked`](@ref); the site is never changed. See [Metronix Sites](metronix.md) |

## Continue in code

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
