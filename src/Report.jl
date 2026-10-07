# Report.jl - a Markdown record of the processing of a site.
# Author: @pankajkmishra
#
# The EDI file holds the transfer function. The record next to it holds how
# it was made, for a reader: the set-up (the base and remote sites), the data
# that each rate used, the sensors with their calibration files and dipole
# lengths (and where each length came from), the files that gave them, all the
# options, the check of the channels, and a table of the estimates. Export
# writes it as <site>.md beside <site>.edi

_md_cell(x::Real) = isfinite(x) ? @sprintf("%.4g", x) : "–"
_md_pm(x::Real, e::Real) = isfinite(x) ? (isfinite(e) ? @sprintf("%.4g ± %.2g", x, e) : @sprintf("%.4g", x)) : "–"

function _md_value(v)
    v === nothing && return "–"
    v isa AbstractString && return isempty(v) ? "–" : "`$v`"
    v isa Bool && return v ? "yes" : "no"
    v isa Real && return isinf(v) ? "none" : string(v)
    v isa Symbol && return string(v)
    v isa AbstractDict && return isempty(v) ? "–" : join(("$k = $(v[k])" for k in sort!(collect(keys(v)); by = string)), ", ")
    v isa Tuple && length(v) == 2 && return "$(v[1]) to $(v[2])"
    v isa AbstractVector && return isempty(v) ? "–" : join(_md_value.(v), ", ")
    return string(v)
end

const _MD_OPTIONS = [
    (:window, "FFT window (samples)"), (:overlap, "Window overlap"), (:bands_per_decade, "Bands per decade"),
    (:prewhiten, "AR prewhitening order"), (:nyquist_fraction, "Top frequency (× Nyquist)"),
    (:min_harmonic, "Lowest harmonic (last level)"), (:min_windows, "Windows a level needs"),
    (:max_levels, "Most decimation levels"), (:min_period, "Shortest period (s)"),
    (:max_period, "Longest period (s)"), (:method, "Method"), (:huber, "Huber threshold"),
    (:leverage, "Leverage weights"),
    (:jackknife_groups, "Jackknife blocks"),
    (:calibration, "Calibration directories"), (:dipole, "Dipole lengths (m)"),
    (:azimuths, "Sensor azimuths (°)"), (:masks, "Masked intervals"), (:span, "Time span"),
]

"""
    write_tf_report(path, tf::TransferFunction; files = String[]) -> String

Write a Markdown record of how `tf` was made: the site and its set-up (base
and remote sites), the data of each rate (hours, spans, start and end, the
decimation levels), each sensor with its calibration file or its dipole
length and where the length came from (`d/dipoles.dat`, the header or the
default), the calibration files and dipole tables that the estimate used,
all the processing options, the check of the channels
([`check_polarity`](@ref)), the files of the export and a table of the
estimates (apparent resistivity and phase with their errors, the tipper, the
coherences and the windows at each period). A PNG among `files` shows in the
record. The function returns `path`.
"""
function write_tf_report(path::AbstractString, tf::TransferFunction; files = String[])
    md = tf.metadata
    rho, rerr = apparent_resistivity(tf)
    phi, perr = impedance_phase(tf)
    check = check_polarity(tf; magnetic = get(md, :flipcheck, false))
    mkpath(dirname(abspath(path)))
    open(path, "w") do io
        println(io, "# ", tf.site, " – transfer function\n")
        @printf(io, "Processed %s UTC with Timekeepers.jl %s (Julia %s).\n\n",
                Dates.format(now(UTC), "yyyy-mm-dd HH:MM"), pkgversion(Timekeepers), VERSION)

        println(io, "## Site\n")
        println(io, "| | |\n|:---|:---|")
        pos = isfinite(tf.latitude) ? @sprintf("%.5f°, %.5f°, %.0f m", tf.latitude, tf.longitude, tf.elevation) : "unknown"
        println(io, "| Position | ", pos, " |")
        println(io, "| Set-up | ", _setup_text(tf), " |")
        println(io, "| Inputs Hx, Hy from | ", isempty(tf.base) ? tf.site : tf.base, " |")
        println(io, "| Reference | ", isempty(tf.remote) ? "the inputs (no remote site)" : join(tf.remote, ", "), " |")
        w = get(md, :witnesses, String[])
        isempty(w) || println(io, "| Witnesses of H | ", join(w, ", "), " (for the channel check only) |")
        haskey(md, :start) && println(io, "| Data | ", Dates.format(md[:start], "yyyy-mm-dd HH:MM:SS"), " to ",
                                      Dates.format(md[:stop], "yyyy-mm-dd HH:MM:SS"), " |")
        println(io, "| Sign convention | exp(+iωt), Z in (mV/km)/nT, x north, y east |")
        println(io)

        if haskey(md, :rates)
            println(io, "## Data and sensors\n")
            for (fs, info) in sort!(collect(md[:rates]); by = first, rev = true)
                @printf(io, "**%s**: %.2f h in %d span%s, %d decimation level%s, %s to %s.\n\n", _fs_label(fs),
                        info[:hours], info[:spans], info[:spans] == 1 ? "" : "s", info[:levels],
                        info[:levels] == 1 ? "" : "s", Dates.format(info[:start], "yyyy-mm-dd HH:MM:SS"),
                        Dates.format(info[:stop], "yyyy-mm-dd HH:MM:SS"))
                nq = get(info, :nyquist_fraction, NaN)
                edge = get(info, :edge, NaN)
                if isfinite(nq)
                    @printf(io, "Top frequency %.4g Hz (%.2f of Nyquist)", nq * fs / 2, nq)
                    isfinite(edge) ? @printf(io, ": the coherence of E with H ends at %.4g Hz.\n\n", edge) :
                                     println(io, ".\n")
                end
                println(io, "| Channel | Sensor |\n|:---|:---|")
                for (k, v) in sort!(collect(info[:sensors]); by = first)
                    println(io, "| ", k, " | ", v, " |")
                end
                println(io)
            end
            used = sort!(unique(String[f for info in values(md[:rates]) for f in get(info, :files, String[])]))
            if !isempty(used)
                println(io, "Files that gave the calibrations and the dipole lengths:
")
                foreach(f -> println(io, "- `", f, "`"), used)
                println(io)
            end
            fails = get(md, :failures, String[])
            isempty(fails) || (println(io, "Rates without an estimate:\n"); foreach(f -> println(io, "- ", f), fails); println(io))
        end

        if haskey(md, :options)
            o = md[:options]
            println(io, "## Processing options\n")
            println(io, "| Option | Value |\n|:---|:---|")
            for (key, label) in _MD_OPTIONS
                hasproperty(o, key) || continue
                v = getproperty(o, key)
                text = key === :calibration && v === nothing ?
                       "s/ of the survey and directories with \"cal\" in their names, in the site directory and up to three above it (see the sensors)" : _md_value(v)
                println(io, "| ", label, " | ", text, " |")
            end
            println(io)
        end

        println(io, "## Result\n")
        if !isempty(tf.periods)
            rates = sort!(unique(filter(isfinite, tf.sample_rate)); rev = true)
            @printf(io, "%d periods from %.4g s to %.4g s", length(tf.periods), first(tf.periods), last(tf.periods))
            isempty(rates) || print(io, ", from ", join(_fs_label.(rates), " and "))
            println(io, ". The errors are the estimated errors, with no floor.\n")
        end
        println(io, "**Channel check:** ", check.message, "\n")

        if !isempty(files)
            println(io, "## Files\n")
            foreach(f -> println(io, "- `", basename(f), "`"), files)
            println(io)
            png = findfirst(f -> endswith(lowercase(f), ".png"), files)
            png === nothing || println(io, "![", tf.site, "](", basename(files[png]), ")\n")
        end

        println(io, "## Estimates\n")
        println(io, "ρa in Ω·m, φ in degrees, errors one standard deviation. Coherences of Ex, Ey and Hz with the inputs.\n")
        println(io, "| T (s) | rate | windows | ρxy | φxy | ρyx | φyx | Re Tzx | Im Tzx | Re Tzy | Im Tzy | coh Ex | coh Ey | coh Hz |")
        println(io, "|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|")
        for k in eachindex(tf.periods)
            cells = [_md_cell(tf.periods[k]), isfinite(tf.sample_rate[k]) ? _fs_label(tf.sample_rate[k]) : "–",
                     string(tf.n_windows[k]),
                     _md_pm(rho[1, 2, k], rerr[1, 2, k]), _md_pm(phi[1, 2, k], perr[1, 2, k]),
                     _md_pm(rho[2, 1, k], rerr[2, 1, k]), _md_pm(phi[2, 1, k], perr[2, 1, k]),
                     _md_cell(real(tf.T[1, k])), _md_cell(imag(tf.T[1, k])),
                     _md_cell(real(tf.T[2, k])), _md_cell(imag(tf.T[2, k])),
                     _md_cell(tf.coherence[1, k]), _md_cell(tf.coherence[2, k]), _md_cell(tf.coherence[3, k])]
            println(io, "| ", join(cells, " | "), " |")
        end
    end
    return path
end
