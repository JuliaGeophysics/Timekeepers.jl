# TKProc in Detail

The window calls functions that you can use without it:

```julia
using Timekeepers

survey = scan_survey("data/survey")
tf = estimate_tf(survey, "site004"; remote = :auto, rate = 128)

flip_check!(survey, tf).message      # which channel, if any, is flipped
export_tf("out", tf)                 # EDI, ModEM, PNG and the record
```

- No `base` and no `remote` gives a single-site estimate. `base` takes Hx, Hy from another site, and
  `remote` gives the reference.
- The keywords of [`estimate_tf`](@ref) are the options of the window.
- Coil calibrations come from a directory with "cal" in its name near the site, or from
  `calibration`.
- `Z` is in (mV/km)/nT with `exp(+iωt)`. The errors have no floor.
