```@raw html
---
layout: home

hero:
  name: "Timekeepers.jl"
  tagline: Timeseries analysis in Electromagnetic Geophysics
  actions:
    - theme: brand
      text: Getting Started
      link: /getting_started
    - theme: alt
      text: GitHub
      link: https://github.com/JuliaGeophysics/Timekeepers.jl
---

<div class="tk-features">
  <a class="tk-feature" href="./tkdash">
    <h2>TKDash</h2>
    <p>Scan a survey. See which sites recorded together, and find the base and remote sites of each site.</p>
    <img src="./assets/TKDash.png" alt="TKDash showing the sites of a survey on a map">
    <span class="tk-link">Open TKDash guide →</span>
  </a>
  <a class="tk-feature" href="./tkapp">
    <h2>TKApp</h2>
    <p>Open a recording. Scroll through it, mask bad intervals, check the spectra and write clean data.</p>
    <img src="./assets/TK.png" alt="TKApp showing a five-channel Metronix record">
    <span class="tk-link">Open TKApp guide →</span>
  </a>
  <a class="tk-feature" href="./tkproc">
    <h2>TKProc</h2>
    <p>Process a site. Estimate its impedance and tipper with base and remote sites, and write EDI and ModEM files.</p>
    <img src="./assets/TKproc.png" alt="TKProc showing the apparent resistivity, phase and tipper of a site">
    <span class="tk-link">Open TKProc guide →</span>
  </a>
</div>
```

Timekeepers.jl reads magnetotelluric and geomagnetic recordings (LEMI-424, GEOMAG-02, Metronix ADU),
lets you mark bad intervals, writes clean data back in the same format, and estimates transfer
functions. It is part of the
[JuliaGeophysics](https://github.com/JuliaGeophysics) ecosystem.

```bash
git clone https://github.com/JuliaGeophysics/Timekeepers.jl
cd Timekeepers.jl
julia --project=. -e "using Pkg; Pkg.instantiate()"
```

New to Julia? [Getting Started](getting_started.md) goes through each step, from installing Julia
to the first window.

To cite Timekeepers.jl, cite the repository: <https://github.com/JuliaGeophysics/Timekeepers.jl>.
