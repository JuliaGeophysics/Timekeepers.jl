# TKProc

TKProc estimates the impedance and the tipper of one site, with its base and remote sites, and
writes them as EDI and ModEM files.

!!! warning "Experimental"
    The processing is new. Its options and its output can change.

```julia
using Timekeepers

run_tkproc("data/survey/site004")   # start Julia with `julia -t auto`
```

## Use

| Action | How |
|:---|:---|
| Open a site | **Load Site…** |
| Choose the set-up | **Base**, **Remote** and **Rate** |
| Estimate | **Process** |
| Find a flipped channel | **FlipCheck** |
| Save | **Export…** writes `<site>.edi`, `<site>.dat`, `<site>.png` and `<site>.md` |

TKProc uses the plan that you exported from [TKDash](tkdash.md), if it finds one in the survey.
