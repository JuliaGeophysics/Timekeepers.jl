# TKProc in Detail

The window calls functions that you can use without it:

```julia
using Timekeepers

survey = scan_survey("data/survey")
tf = estimate_tf(survey, "site004"; remote = :auto, rate = 128)

flip_check!(survey, tf).message      # which channel, if any, is flipped
export_tf("out", tf)                 # <name>.edi and <name>.png, see tf_filename
```

- No `base` and no `remote` gives a single-site estimate. `base` takes Hx, Hy from another site, and
  `remote` gives the reference.
- The keywords of [`estimate_tf`](@ref) are the options of the window.
- A survey directory holds the coil calibrations in `s/` and the dipole table in
  `d/dipoles.dat`. Coil calibrations come from `s/` or a directory with "cal" in its name
  near the site, or from `calibration`.
- Dipole lengths come from `dipole`, then `d/dipoles.dat` (one line per site with the
  distances in metres from the centre to the N, S, E and W electrodes, Ex = N + S,
  Ey = E + W), then the `.ats` header. Without any of them, or for a `-` distance, each
  electrode is 50 m from the centre. The INFO block of the EDI file lists the source of each
  length and every calibration file and dipole table that the estimate used.

  ```
  site       N      S      E      W
  site002  24.5  48.5  50.9  51.5
  ```
- `Z` is in (mV/km)/nT with `exp(+iωt)`. The errors have no floor.
