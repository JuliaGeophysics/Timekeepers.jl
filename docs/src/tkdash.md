# TKDash

See which sites recorded together, and pick the base and remote sites of each site.

![TKDash with site002 in focus: two base sites in blue, three remote sites in amber](assets/TKDash.png)

```julia
using Timekeepers
run_tkdash("data/survey")          # or run_tkdash() to pick the folder
```

| To | Do |
|:---|:---|
| Focus a site | click it on the map, or **< Prev** / **Next >** |
| See the whole survey | **Overview** |
| Choose the sampling rate | **Rate** |
| Set the distances (km) | **Base ≤** (interstation) and **Remote ≥** (remote reference) |
| Drop a base or remote site | right-click it; **Restore** brings it back |
| Save the plan for TKProc | **Export…** |

Mixed Metronix, LEMI-424 and GEOMAG surveys work. Next: clean a site in [TKApp](tkapp.md).
