# TKProc

TKProc estimates the impedance and the tipper of one site, with its base and remote sites, and
writes them as EDI and ModEM files.

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
| Choose the set-up | **Base**, **Remote** and **Rate** |
| Estimate | **Process** |
| Find a flipped channel | **FlipCheck** |
| Save | **Export…** writes `<site>.edi`, `<site>.dat`, `<site>.png` and `<site>.md` |

TKProc uses the plan that you exported from [TKDash](tkdash.md), if it finds one in the survey.
