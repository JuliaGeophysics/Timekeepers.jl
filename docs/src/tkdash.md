# TKDash

TKDash scans a survey and shows, for each site, which other sites recorded at the same time and the
same sampling rate:

- **base** sites: within the **Base ≤** distance, for interstation processing
- **remote** sites: beyond the **Remote ≥** distance, for remote-reference processing

![TKDash with site002 in focus: two base sites in blue, three remote sites in amber](assets/TKDash.png)

```julia
using Timekeepers

run_tkdash("data/survey")          # or run_tkdash() to pick the folder
```

## Use

| Action | How |
|:---|:---|
| Focus a site | click it on the map, or **< Prev** / **Next >** |
| Back to the whole survey | **Overview** |
| Choose the sampling rate | **Rate** |
| Change the distances | **Base ≤** and **Remote ≥** (km) |
| Drop a base or remote site | right-click it; **Restore** brings it back |
| Save the plan | **Export…** |

One survey can mix Metronix, LEMI-424 and GEOMAG sites. The plan file lists the base and remote sites of every site. Open a site in [TKApp](tkapp.md) to
clean it.
