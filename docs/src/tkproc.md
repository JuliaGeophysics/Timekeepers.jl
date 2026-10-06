# TKProc Transfer Functions

TKProc estimates the impedance and the tipper of one magnetotelluric site at a
time. It finds the base and remote sites of the site in the survey around it,
processes the site with the set-up that you select, checks the channels for a
flipped sign and writes the result as EDI, ModEM, a plot and a record of the
processing.

!!! warning "Experimental"
    The processing is new. Its options and its output can change.

## Open the window

```julia
using Timekeepers

proc = run_tkproc("data/survey/site004")   # blocks until the window closes
export_tf("out", proc.results["site004"])  # the estimates survive the session
```

Start Julia with more than one thread (`julia -t auto`). The estimate then
runs on another thread, and the window stays live. A launcher accepts the site
directory as its argument:

```bash
julia --project=. -t auto examples/tkproc.jl /path/to/survey/site004
```

[`run_tkproc`](@ref) scans the directory above the site with
[`scan_survey`](@ref). If you do not give a site, the window opens empty: use
**Load Site…**.

If you exported a plan from [TKDash](tkdash.md) (`reference_plan.txt` in the
survey directory or up to two directories above it), TKProc uses it. The base
and remote menus then hold only the sites of the plan. The **TKDash plan**
switch turns it off.

## The window

| Control | What it does |
|:---|:---|
| **Load Site…** | Open a site directory and scan the survey around it |
| **Base** | The site that gives the inputs Hx, Hy. With none, the site gives its own |
| **Remote** | The far site whose Hx, Hy are the reference. With none, the estimate is single site |
| **Rate** | The sampling rate to process. `All rates` keeps the estimate with the smallest error at each period |
| **Process** | Estimate the site with the menus and the options |
| **FlipCheck** | Compare Hx, Hy of the inputs with up to two base sites to find a flipped channel |
| **Export…** | Write the EDI, the ModEM file, the plot and the record of the site into a directory |
| **Clear** | Empty the screen |

Two rows below the header hold the processing options:

- the window length, the overlap, the bands in each decade, the AR order, the
  top frequency, the lowest harmonic, the fewest windows and the periods;
- the method (`CT2004` or `EB1986`), the Huber threshold, the jackknife blocks
  and the leverage weights.

The ⓘ next to each option explains it. The keywords of
[`estimate_tf`](@ref) give each option and its default.

The plots show the apparent resistivity above the phase, and Tzx above Tzy,
against the period. Zxy is red and Zyx blue. The **Full tensor** switch adds
Zxx (green) and Zyy (lilac), and the **Error bars** switch hides or shows the
errors. The line below the plots gives the check of the channels, and the
status line shows each step of the processing.

## Without the window

```julia
survey = scan_survey("data/survey")
tf = estimate_tf(survey, "site004"; remote = :auto, rate = 128)

rho, rho_err = apparent_resistivity(tf)
phi, phi_err = impedance_phase(tf)
flip_check!(survey, tf).message      # which channel, if any, is flipped

write_edi("site004.edi", tf)
write_modem("survey.dat", [tf]; components = :offdiagonal)
export_tf("out", tf)                 # EDI, ModEM, PNG and the record
```

The set-ups of [`estimate_tf`](@ref):

- no `base`, no `remote`: single site;
- `base`: the inputs come from a base site (for a site without magnetic
  channels, or with poor ones);
- `remote`: one or more remote sites give the reference;
- `base` and `remote` together.

[`default_references`](@ref) gives the set-up that TKProc selects when it
loads a site. For coils, the processing reads the Metronix calibration files
([`read_calibration`](@ref)) from the directory that you give, or from a
directory with "cal" in its name near the site.

`Z` is in (mV/km)/nT with the time dependence `exp(+iωt)`. The errors are the
estimated errors, with no floor. Set the error floors of an inversion in its
own set-up.
