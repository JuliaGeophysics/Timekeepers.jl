# TKDash Survey Dashboard

Processing a magnetotelluric site against others needs to know which sites
were recording **at the same time and the same sampling rate**, and how far
away they stood. TKDash scans a survey directory and shows exactly that, for
one site at a time:

- the **site** in focus;
- its **base** sites — sites that recorded with it, within a chosen distance
  (for interstation processing);
- its **remote** sites — sites that recorded with it, beyond a second, larger
  distance (for remote-reference processing).

![TKDash with site002 in focus: two base sites in blue, three remote sites in amber](assets/TKDash.png)

## Opening the dashboard

```julia
using Timekeepers

dash = run_tkdash("data/survey")       # blocks until the window closes
write_reference_plan("plan.txt", dash) # the choices made in the window survive the session
```

[`run_tkdash`](@ref) without a path opens a folder dialog. A bundled launcher
takes the survey directory as its argument:

```bash
julia --project=. examples/tkdash.jl /path/to/survey
```

Metronix, LEMI-424 and GEOMAG sites can sit side by side in one survey. The
scan reads **headers only** — the `.ats` headers of a Metronix site, and
the first and last lines of LEMI-424 and GEOMAG files — so a large survey
opens in a second or two. A site is any directory holding Metronix `meas_*`
directories or LEMI-424 / GEOMAG files; directories that are not sites are
searched four levels deep. Rate directories split off a site (`site002.128`
beside `site002`) are skipped, so nothing is counted twice.

## The window

The window follows the layout of MTGeophysics' data dashboard: a header, then
the map, a thin button that collapses it, and the charts.

| Control | What it does |
|:---|:---|
| **\|< · < Prev · site · Next > · >\|** | Step through the sites, or pick one |
| **Overview** | Back to the whole survey, no site in focus |
| **Restore** | Bring back the sites dropped for the site in focus, or for every site in the overview |
| **Rate** | Which recordings pair: `All rates` — any two sites recording at the same time, whatever their rates (for mixed LEMI, GEOMAG and Metronix surveys); `Shared rate` — only at a rate both sites recorded; or one rate. A mixed survey opens at `All rates`, any other at the rate with the most overlap |
| **Base ≤** | Sites recording with the site within this distance (km) are base sites |
| **Remote ≥** | Sites recording with the site this far (km) or farther are remote sites; sites between the two distances are neither |
| **Overlap ≥** | Sites recording with the site for less than this (hours) are neither |
| **Open…** | Scan another directory |
| **Export…** | Write every site's base and remote sites as a text table (see below) |

**Overview.** One chart shows every run of every site, one row per site, with
each site's recording hours at the end of its bar; runs at other rates are faint. The map
marks every site.

**A site in focus.** Click a site on the map or a row of a chart, or step to
it, and everything turns to it. The map marks the site with a **magenta
star**, its base sites with **blue** circles and its remote sites with
**amber** circles, each darker the longer it recorded with the site; the rest
of the survey stays as small grey dots, and a legend names them. The charts
split in two: **base** — the site over its base sites — and **remote** — the
site over its remote sites — with the time each recorded with the site
painted in its shade. Where each bar ends, `[1.5h/9.5h]` gives the hours a
base or remote site recorded with the site, then its own hours; the site's
own row reads `[22.4h]`. Only the best five base and five remote sites are
shown.

**Common window.** A light teal band, labelled `common [9.0h]`, marks the
time the site and *all* the base and remote sites it keeps were recording
together — the stretch processing can use with every one of them at once.
It counts every kept site, those beyond the five on screen too, so dropping a
site that recorded only briefly widens it. The export writes it for each
site, and [`common_window`](@ref) computes it in code.

**Hovering** over a site on the map or a row of a chart shows its name, its
role and distance from the site, the hours it recorded with the site and on
its own, and, over a bar, the run's file, rate and time span.

**Electric-only sites.** Base and remote sites lend their magnetic field, so
only sites that recorded Hx and Hy are ever base or remote. A telluric site —
electric channels only — can still be the site in focus, paired with the
magnetic fields of its base and remote sites; its chart rows read
`site003 · E only`, and hovering any site lists its channels. For LEMI-424 and
GEOMAG files a channel that is zero or blank at both ends of the file counts
as not recorded.

**Dropping a site.** Right-click a base or remote site, on the map or a chart,
to drop it from the site's lists — a sensor you know was faulty, say. A
dropped site leaves the comparison: it disappears from the charts, the next
best site takes its place, and it no longer narrows the common window or
appears in the export. It stays on the map as a hollow circle; right-click it
there to bring it back, or press **Restore** to bring back every site
dropped for the site in focus. Each site keeps its own drops. The line under the
charts says what to do in grey, confirms each drop, scan and export in cyan,
and reports a failure in red.

**Zooming.** The map keeps the survey's shape at true scale and grows with the
window. Scroll over it to zoom, drag a box with the left button to zoom to the
box, and drag with the right button to pan — useful when a remote site is far
and the base sites crowd together. The zoom stays as you step between sites;
the **Reset Zoom** button under the map brings back the whole survey and
reframes the charts. Scroll over a chart to zoom its
time axis. The `‹` button collapses the map and gives the charts the full
width.

## The plan file

**Export…** and [`write_reference_plan`](@ref) write a plain-text table — one
header line, one row per site, columns aligned:

```
site      base                        overlap (h)          remote                                        overlap (h)
site002   site004, site006            11.77, 9.00          site099, site100                              11.77, 11.77
site004   site009, site002, site006   11.87, 11.77, 9.00   site099, site100                              12.00, 12.00
site099   -                           -                    site100, site004, site009, site002, site006   22.98, 12.00, 11.87, 11.77, 9.00
```

Base and remote sites are listed **longest overlap first**, so the first is the
best, and each `overlap (h)` column gives their hours with the site in the
same order; `-` marks an empty list. The table holds every base and remote
site the site keeps, not only the five the window shows.
[`read_reference_plan`](@ref) reads it back as one
`(site, base, base_hours, remote, remote_hours)` per row, for processing
steps that pair a site with its base or remote sites.

## Without the window

```julia
survey = scan_survey("data/survey")
survey["site004"]                                  # one site, its runs and position
refs = site_references(survey, "site004"; rate = 128, base_km = 5, remote_km = 20)
refs.base                                          # longest overlap first
overlap_intervals(survey["site004"], survey["site100"]; rate = 128)
write_reference_plan("plan.txt", survey; rate = 128, base_km = 5, remote_km = 20)
```

A LEMI-424 or GEOMAG file is taken as recording from its first line to its
last, so a gap inside one file is not seen by the scan; split such records
into files at the gap, or check them in [TKApp](tkapp.md).
