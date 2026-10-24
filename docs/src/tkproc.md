# TKProc

TKProc estimates the impedance and the tipper of one site, with its base and remote sites, and
writes them as an EDI file.

!!! warning "Experimental"
    The processing is new. Its options and its output can change.

![TKProc showing the apparent resistivity, phase and tipper of a site](assets/TKproc.png)

```julia
using Timekeepers

run_tkproc("data/survey")   # opens the first site; start Julia with `julia -t auto`
```

## Use

| Action | How |
|:---|:---|
| Open a survey | **Survey…** |
| Open a site with its base and remote sites | **Site** menu |
| Choose the set-up | **Base**, **Remote**, **Rate** and the options; **Save set-up** writes them to `procset.dat` of the site |
| Estimate | **Process**; it writes the set-up that it used to `procset.dat` of the site |
| Save | **Export…** writes `<name>.edi` and `<name>.png` |

TKProc uses the plan that you exported from [TKDash](tkdash.md), if it finds one in the survey.

The minimum distance of the remote sites in "all ≥ …" (`min_distance`), the period range (`min_period`, `max_period`), the time span (`span`) and the rotation of the plots (`rotate`) are set in `procset.dat` of the site only. TKProc reads them when you select the site and each time you press **Process**.

The first entry of **Base** is the site itself: the site gives its own Hx, Hy.

## The survey directory

```
survey/
├── s/              coil calibration files, e.g. MFS07e160.txt
├── d/dipoles.dat   dipole distances of each site
├── site002/
└── site004/
```

`d/dipoles.dat` has one header line, then one line for each site: the distances in metres from the
centre of the site to the N, S, E and W electrodes. Ex is N + S long and Ey is E + W. `-` is a
distance that is not known. The file has no comments.

```
site       N      S      E      W
site002   24.5   48.5   50.9   51.5
site004   49.0   51.2   51.0   51.0
```

TKProc takes each dipole length from `d/dipoles.dat` first, then from the `.ats` header. If neither
gives it, each electrode is 50 m from the centre. LEMI-424 and GEOMAG files have no electrode
positions, so give their sites a line.

## The exported name

`<name>` holds the values that made the estimate, joined by `-`, so that you can make it again:

```
site-base-remote-rate-window-overlap-bands-prewhiten-top-harmonic-windows-levels-method-huber-leverage-jackknife
site002-site002-site099-all-256-0.5-8-3-auto-4-8-12-ct2004-1.5-1-50.edi
```

The base is the site itself when it gives its own Hx, Hy, and the remote is `none` without a remote
site. The INFO block of the EDI file lists each dipole length with where it came from, each
calibration file, and all the options.
