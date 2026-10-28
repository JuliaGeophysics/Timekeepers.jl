# TKDash in Detail

```julia
using Timekeepers

dash = run_tkdash("data/survey")       # blocks until the window closes
write_reference_plan("plan.txt", dash) # keep the choices made in the window
```

`run_tkdash()` without a path opens a folder dialog. From a clone:
`julia --project=. examples/tkdash.jl`.

!!! tip "Tutorial data"
    The public Metronix tutorial data of BRGM (Alsace, IMAGE FP7 project) has six sites that
    recorded together at 128 Hz. The screenshot uses it.

![TKDash with site002 in focus: two base sites in blue, three remote sites in amber](../assets/TKDash.png)

- **Base** sites recorded with the site within **Base ≤** km (interstation processing).
- **Remote** sites recorded with it at **Remote ≥** km or more (remote reference).
- Only sites with Hx and Hy can be base or remote. Electric-only sites show as `E only`.
- The scan reads headers only, so large surveys open in seconds. It looks four levels deep and
  skips rate directories like `site002.128`.

## Controls

| Control | What it does |
|:---|:---|
| **\|< · < Prev · site · Next > · >\|** | Step through or pick a site |
| **Overview** | Whole survey, no site in focus |
| **Rate** | Only runs at this rate pair up; `All rates` pairs any rates. Opens at 128 Hz if present |
| **Base ≤** / **Remote ≥** | Distances in km; sites in between are neither |
| **Overlap ≥** | Minimum hours recorded together |
| **Restore** | Bring back dropped sites |
| **Open…** | Scan another directory |
| **Export…** | Write the plan file |

| Mouse | |
|:---|:---|
| click a site or chart row | focus it |
| right-click a base/remote site | drop it (hollow on the map); right-click again to restore |
| hover | name, role, distance, hours, channels |
| scroll / left-drag box / right-drag on map | zoom / zoom to box / pan |
| scroll on a chart | zoom time |

On the map: **magenta star** = site, **blue** = base, **amber** = remote (darker = longer
overlap). Bars end with `[together/own]` hours. The teal **common** band is when the site and all
its kept sites recorded at once ([`common_window`](@ref) in code). Charts show the five best of each;
the plan keeps them all.

## The plan file

```
site      base                        overlap (h)          remote                                        overlap (h)
site002   site004, site006            11.77, 9.00          site099, site100                              11.77, 11.77
site004   site009, site002, site006   11.87, 11.77, 9.00   site099, site100                              12.00, 12.00
site099   -                           -                    site100, site004, site009, site002, site006   22.98, 12.00, 11.87, 11.77, 9.00
```

Longest overlap first; `-` is an empty list. [`read_reference_plan`](@ref) returns
`(site, base, base_hours, remote, remote_hours)` per row.

## Without the window

```julia
survey = scan_survey("data/survey")
survey["site004"]                                  # one site, its runs and position
refs = site_references(survey, "site004"; rate = 128, base_km = 5, remote_km = 20)
refs.base                                          # longest overlap first
overlap_intervals(survey["site004"], survey["site100"]; rate = 128)
write_reference_plan("plan.txt", survey; rate = 128, base_km = 5, remote_km = 20)
```

The scan can't see a gap inside one LEMI-424 or GEOMAG file (it reads only the first and last
lines). Split the file at the gap, or look at it in [TKApp](tkapp.md).
