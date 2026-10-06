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

<<<<<<< HEAD
```bash
julia --project=. examples/tkdash.jl /path/to/survey
```

One survey can contain Metronix, LEMI-424 and GEOMAG sites together. The scan
reads **only the headers**: the `.ats` headers of a Metronix site, and the
first and the last lines of LEMI-424 and GEOMAG files. Thus, a large survey
opens in one or two seconds.

- A site is a directory that holds Metronix `meas_*` directories or LEMI-424
  or GEOMAG files.
- The scan looks four levels deep in directories that are not sites.
- The scan ignores the rate directories of a site (`site002.128` next to
  `site002`). Thus, it does not count a site two times.

## The window

The window has four parts:

1. a header;
2. the map;
3. a thin button that collapses the map;
4. the charts.

| Control | What it does |
=======
| Action | How |
>>>>>>> origin/main
|:---|:---|
| Focus a site | click it on the map, or **< Prev** / **Next >** |
| Back to the whole survey | **Overview** |
| Choose the sampling rate | **Rate** |
| Change the distances | **Base ≤** and **Remote ≥** (km) |
| Drop a base or remote site | right-click it; **Restore** brings it back |
| Save the plan | **Export…** |

One survey can mix Metronix, LEMI-424 and GEOMAG sites. The plan file lists the base and remote sites of every site. Open a site in [TKApp](tkapp.md) to
clean it.
