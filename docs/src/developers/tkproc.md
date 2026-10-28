# TKProc in Detail

```julia
using Timekeepers

survey = scan_survey("data/survey")
tf = estimate_tf(survey, "site004"; remote = :auto, rate = 128)

export_tf("out", tf)                 # <name>.edi and <name>.png, see tf_filename
```

- No `base` and no `remote`: single-site estimate. `base` takes Hx, Hy from another site;
  `remote` is the reference.
- The keywords of [`estimate_tf`](@ref) are the options of the window.
- Coil calibrations come from `s/`, a nearby directory with "cal" in its name, or `calibration`.
- Dipole lengths come from `dipole`, then `d/dipoles.dat`, then the `.ats` header, else 50 m per
  electrode. The EDI's INFO block records every source used.

  ```
  site       N      S      E      W
  site002  24.5  48.5  50.9  51.5
  ```

- `Z` is in (mV/km)/nT with `exp(+iωt)`. Errors have no floor.
- E = −V/L, where V = V(+) − V(−) and + is the north (east) electrode. A correctly wired site over a
  1D or 2D earth gives φxy in the first quadrant and φyx in the third.
