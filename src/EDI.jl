# EDI.jl - SEG EDI files of transfer functions.
# Author: @pankajkmishra
#
# The SEG EDI standard (Wight 1987) holds the impedance and the tipper of one
# site as text. The file has:
# - >HEAD: the site, its position and the program
# - >INFO: free text; here the set-up and the options of the estimate
# - >=DEFINEMEAS with one >HMEAS or >EMEAS line for each channel
# - >=MTSECT, then one block for each quantity: the frequencies, the rotation,
#   the real and the imaginary part and the variance of each element of Z and
#   of the tipper
# Z is in (mV/km)/nT with the time dependence exp(+iωt). A variance is the
# variance E|δZ|² of the complex element. A missing value is the EMPTY value

const _EDI_EMPTY = 1.0e32

_edi_dms(x::Real) = begin
    s = x < 0 ? "-" : ""
    a = abs(x)
    d = floor(Int, a)
    m = floor(Int, (a - d) * 60)
    sec = ((a - d) * 60 - m) * 60
    @sprintf("%s%d:%02d:%05.2f", s, d, m, sec)
end

function _edi_degrees(x)
    x === nothing && return NaN
    t = strip(x)
    p = abs.(parse.(Float64, split(t, ':')))
    v = sum(p[k] / 60.0^(k - 1) for k in eachindex(p))
    return startswith(t, "-") ? -v : v
end

function _edi_block(io, name, values; extra = "")
    n = length(values)
    println(io, ">", name, extra, " //", n)
    for (k, v) in enumerate(values)
        x = isfinite(v) ? v : _EDI_EMPTY
        @printf(io, "%16.6E", x)
        (k % 6 == 0 || k == n) && println(io)
    end
    println(io)
end

"""
    write_edi(path, tf::TransferFunction) -> String

Write `tf` to a SEG EDI file. The frequencies are in decreasing order. The
file holds `Z` (ZXX … ZYY: real, imaginary, variance), the tipper (TX, TY:
real, imaginary, variance) if the site has Hz, and the coherences of Ex, Ey
and Hz with the inputs. The `>INFO` block lists the set-up (base and remote
sites), the sensors and the options of the estimate. A missing value is
`1.0E+32`. The function returns `path`.
"""
function write_edi(path::AbstractString, tf::TransferFunction)
    n = length(tf.periods)
    n == 0 && error("$(tf.site) has no periods to write")
    order = n:-1:1                                      # decreasing frequency
    freqs = 1 ./ tf.periods[order]
    md = tf.metadata
    lat = isfinite(tf.latitude) ? tf.latitude : 0.0
    lon = isfinite(tf.longitude) ? tf.longitude : 0.0
    elev = isfinite(tf.elevation) ? tf.elevation : 0.0
    acq = get(md, :start, nothing)
    date(t) = t === nothing ? Dates.format(now(UTC), "yyyy-mm-dd") : Dates.format(t, "yyyy-mm-dd")
    dip = get(md, :dipoles, Dict{Symbol, Float64}())
    mkpath(dirname(abspath(path)))
    open(path, "w") do io
        println(io, ">HEAD")
        println(io, "    DATAID=\"", tf.site, "\"")
        println(io, "    ACQBY=\"Timekeepers.jl\"")
        println(io, "    FILEBY=\"Timekeepers.jl\"")
        println(io, "    ACQDATE=", date(acq))
        get(md, :stop, nothing) === nothing || println(io, "    ENDDATE=", date(md[:stop]))
        println(io, "    FILEDATE=", Dates.format(now(UTC), "yyyy-mm-dd"))
        println(io, "    LAT=", _edi_dms(lat))
        println(io, "    LONG=", _edi_dms(lon))
        @printf(io, "    ELEV=%.1f\n", elev)
        println(io, "    UNITS=M")
        println(io, "    STDVERS=\"SEG 1.0\"")
        println(io, "    PROGVERS=\"Timekeepers.jl $(pkgversion(Timekeepers))\"")
        println(io, "    PROGDATE=", Dates.format(now(UTC), "yyyy-mm-dd"))
        println(io, "    MAXSECT=999")
        println(io, "    EMPTY=1.0E+32")
        println(io)
        println(io, ">INFO")
        println(io, "    MAXINFO=999")
        println(io, "    SITE: ", tf.site)
        println(io, "    SETUP: ", _setup_text(tf))
        m = haskey(md, :options) ? md[:options].method : :ct2004
        println(io, "    METHOD: ", m === :eb1986 ? "robust M-estimate with Huber weights (after Egbert & Booker 1986)" :
                                   "two-stage robust estimate with bounded influence (after Chave & Thomson 2004)",
                ", jackknife variances")
        println(io, "    SIGN CONVENTION: exp(+i omega t)")
        w = get(md, :witnesses, String[])
        isempty(w) || println(io, "    H WITNESSES: ", join(w, ", "), " (channel check only)")
        println(io, "    POLARITY CHECK: ", replace(check_polarity(tf; magnetic = get(md, :flipcheck, false)).message, "°" => " deg"))
        if haskey(md, :rates)
            for (fs, info) in sort!(collect(md[:rates]); by = first)
                @printf(io, "    RATE %s: %.2f h in %d spans, %d levels, top frequency %.4g Hz (%.2f of Nyquist)\n",
                        _fs_label(fs), info[:hours], info[:spans], info[:levels], info[:nyquist_fraction] * fs / 2,
                        info[:nyquist_fraction])
                for (k, v) in sort!(collect(info[:sensors]); by = first)
                    println(io, "    SENSOR ", k, ": ", replace(v, "°" => " deg"))
                end
            end
        end
        if haskey(md, :options)
            o = md[:options]
            println(io, "    WINDOW: ", o.window, " samples, overlap ", o.overlap, ", AR prewhitening order ", o.prewhiten)
            println(io, "    BANDS: ", o.bands_per_decade, " per decade, Huber ", o.huber, ", jackknife errors")
            println(io, "    RANGE: top frequency ", o.nyquist_fraction, " of Nyquist, lowest harmonic ",
                    o.min_harmonic, ", ", o.min_windows, " windows or more")
        end
        println(io)
        println(io, ">=DEFINEMEAS")
        println(io, "    MAXCHAN=", 5 + 2 * length(tf.remote))
        println(io, "    MAXRUN=999")
        println(io, "    MAXMEAS=9999")
        println(io, "    UNITS=M")
        println(io, "    REFTYPE=CART")
        println(io, "    REFLAT=", _edi_dms(lat))
        println(io, "    REFLONG=", _edi_dms(lon))
        @printf(io, "    REFELEV=%.1f\n", elev)
        println(io)
        lx = get(dip, :e1, 0.0) / 2
        ly = get(dip, :e2, 0.0) / 2
        println(io, ">HMEAS ID=1001.001 CHTYPE=HX X=0.0 Y=0.0 Z=0.0 AZM=0.0")
        println(io, ">HMEAS ID=1002.001 CHTYPE=HY X=0.0 Y=0.0 Z=0.0 AZM=90.0")
        println(io, ">HMEAS ID=1003.001 CHTYPE=HZ X=0.0 Y=0.0 Z=0.0 AZM=0.0")
        @printf(io, ">EMEAS ID=1004.001 CHTYPE=EX X=%.1f Y=0.0 Z=0.0 X2=%.1f Y2=0.0\n", -lx, lx)
        @printf(io, ">EMEAS ID=1005.001 CHTYPE=EY X=0.0 Y=%.1f Z=0.0 X2=0.0 Y2=%.1f\n", -ly, ly)
        for k in eachindex(tf.remote)
            @printf(io, ">HMEAS ID=%d.001 CHTYPE=RX X=0.0 Y=0.0 Z=0.0 AZM=0.0\n", 1004 + 2k)
            @printf(io, ">HMEAS ID=%d.001 CHTYPE=RY X=0.0 Y=0.0 Z=0.0 AZM=90.0\n", 1005 + 2k)
        end
        println(io)
        println(io, ">=MTSECT")
        println(io, "    SECTID=\"", tf.site, "\"")
        println(io, "    NFREQ=", n)
        println(io, "    HX=1001.001")
        println(io, "    HY=1002.001")
        println(io, "    HZ=1003.001")
        println(io, "    EX=1004.001")
        println(io, "    EY=1005.001")
        if !isempty(tf.remote)
            println(io, "    RX=1006.001")
            println(io, "    RY=1007.001")
        end
        println(io)
        _edi_block(io, "FREQ", freqs; extra = " ORDER=DEC")
        rot = fill(get(md, :rotation, 0.0), n)
        _edi_block(io, "ZROT", rot)
        for (e, (i, j)) in enumerate(((1, 1), (1, 2), (2, 1), (2, 2)))
            name = _Z_NAMES[e]
            z = tf.Z[i, j, order]
            _edi_block(io, name * "R", real.(z); extra = " ROT=ZROT")
            _edi_block(io, name * "I", imag.(z); extra = " ROT=ZROT")
            _edi_block(io, name * ".VAR", tf.Z_var[i, j, order]; extra = " ROT=ZROT")
        end
        if any(isfinite, tf.T)
            _edi_block(io, "TROT", rot)
            for (j, name) in enumerate(_T_NAMES)
                t = tf.T[j, order]
                _edi_block(io, name * "R.EXP", real.(t); extra = " ROT=TROT")
                _edi_block(io, name * "I.EXP", imag.(t); extra = " ROT=TROT")
                _edi_block(io, name * "VAR.EXP", tf.T_var[j, order]; extra = " ROT=TROT")
            end
        end
        for (k, name) in ((1, "EX"), (2, "EY"), (3, "HZ"))
            any(isfinite, tf.coherence[k, :]) || continue
            _edi_block(io, name * ".COH", tf.coherence[k, order]; extra = " ROT=ZROT")
        end
        println(io, ">END")
    end
    return path
end

"""
    read_edi(path) -> TransferFunction

Read the impedance and the tipper of a SEG EDI file (the blocks that
[`write_edi`](@ref) writes, and the same blocks from other programs). The
values equal to the `EMPTY` value of the file are `NaN`.
"""
function read_edi(path::AbstractString)
    text = read(path, String)
    head(key) = (m = match(Regex("^\\s*$key\\s*=\\s*\"?([^\"\\n]*)\"?", "m"), text); m === nothing ? nothing : strip(m.captures[1]))
    empty = something(tryparse(Float64, something(head("EMPTY"), "1.0E+32")), _EDI_EMPTY)
    blocks = Dict{String, Vector{Float64}}()
    name = nothing
    for line in eachline(IOBuffer(text))
        s = strip(line)
        if startswith(s, ">")
            m = match(r"^>\s*([A-Za-z0-9.]+)", s)
            name = m === nothing ? nothing : uppercase(m.captures[1])
            name !== nothing && occursin("//", s) && (blocks[name] = Float64[])
            name !== nothing && !occursin("//", s) && (name = nothing)
            continue
        end
        name === nothing && continue
        for tok in split(s)
            v = tryparse(Float64, tok)
            v === nothing || push!(blocks[name], v)
        end
    end
    haskey(blocks, "FREQ") || error("No >FREQ block in $path")
    f = blocks["FREQ"]
    clean(v) = [abs(x) >= 0.99 * empty ? NaN : x for x in v]
    get_block(k) = haskey(blocks, k) ? clean(blocks[k]) : fill(NaN, length(f))
    order = sortperm(f; rev = true)
    n = length(f)
    Z = fill(complex(NaN, NaN), 2, 2, n)
    ZV = fill(NaN, 2, 2, n)
    T = fill(complex(NaN, NaN), 2, n)
    TV = fill(NaN, 2, n)
    for (e, (i, j)) in enumerate(((1, 1), (1, 2), (2, 1), (2, 2)))
        nm = _Z_NAMES[e]
        Z[i, j, :] = complex.(get_block(nm * "R"), get_block(nm * "I"))[order]
        ZV[i, j, :] = get_block(nm * ".VAR")[order]
    end
    for (j, nm) in enumerate(_T_NAMES)
        T[j, :] = complex.(get_block(nm * "R.EXP"), get_block(nm * "I.EXP"))[order]
        TV[j, :] = get_block(nm * "VAR.EXP")[order]
    end
    coh = fill(NaN, 3, n)
    for (k, nm) in ((1, "EX"), (2, "EY"), (3, "HZ"))
        coh[k, :] = get_block(nm * ".COH")[order]
    end
    site = something(head("DATAID"), splitext(basename(path))[1])
    elev = something(tryparse(Float64, something(head("ELEV"), "NaN")), NaN)
    return TransferFunction(String(site), :single, "", String[], _edi_degrees(head("LAT")), _edi_degrees(head("LONG")), elev,
                            1 ./ f[order], Z, ZV, T, TV, coh, fill(NaN, 2, n), zeros(Int, n), zeros(Int, n),
                            fill(NaN, n), Dict{Symbol, Any}(:source => abspath(path)))
end
