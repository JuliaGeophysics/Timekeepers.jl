# TKApp

TKApp shows a recording, lets you mask bad intervals, and writes the clean data in the original
format.

![TKApp showing a five-channel Metronix record](assets/TK.png)

```julia
using Timekeepers

run_tkapp()                         # empty; load from the toolbar
run_tkapp("data/LEMI090.txt")       # one file (LEMI-424 or GEOMAG-02)
run_tkapp("data/RK137")             # a Metronix site directory
```

## Use

| Action | How |
|:---|:---|
| Select an interval | left-drag on a trace |
| Mask the selection | right-click, or **Mask** |
| Undo a mask | select it, then **Unmask** |
| Change the visible span | **Window**, then **Scroll** below the plots |
| Show spectra | **View**: `Time \| Spectra` or `Spectra` |
| Save | **Write** |

## What Write saves

- **LEMI-424 and GEOMAG-02**: `<name>_clean.<ext>` in the original format, with `NaN` in masked
  rows, and the mask as `<name>_mask.csv`, next to the source.
- **Metronix**: the original site is never changed. The masked intervals are cut out of a copy named
  `<site>.<rate>` next to it.
