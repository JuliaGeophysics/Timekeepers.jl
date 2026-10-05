# ModEM.jl - ModEM data files.
# Author: @pankajkmishra
#
# The ModEM inversion reads its data from a text file with one block for each
# data type. This file writes the impedance (Full_Impedance) and the tipper
# (Full_Vertical_Components) of a set of sites. Each block has six header
# lines: the type, the time convention, the units, the orientation of the
# axes, the origin of the grid and the counts of periods and sites. Then each
# line holds one element at one period:
#
#   period code lat lon x y z component real imag error
#
# x is north and y is east, in metres from the origin. The error is the
# standard deviation of the real part (and of the imaginary part). It also
# reads such a file again, back into TransferFunction values

const _MODEM_COLUMNS = "# Period(s) Code GG_Lat GG_Lon X(m) Y(m) Z(m) Component Real Imag Error"
const _Z_NAMES = ("ZXX", "ZXY", "ZYX", "ZYY")
const _T_NAMES = ("TX", "TY")

# WGS84
const _WGS84_A = 6378137.0
const _WGS84_E2 = 6.69437999014e-3

# Metres north and east of (lat0, lon0) on the WGS84 ellipsoid, with the radii
# of curvature at the origin. Over a survey of tens of kilometres, the error
# is below a metre for each kilometre
function _local_xy(lat, lon, lat0, lon0)
    φ = deg2rad(lat0)
    s2 = sin(φ)^2
    N = _WGS84_A / sqrt(1 - _WGS84_E2 * s2)
    M = _WGS84_A * (1 - _WGS84_E2) / (1 - _WGS84_E2 * s2)^1.5
    return (M * deg2rad(lat - lat0), N * cos(φ) * deg2rad(lon - lon0))
end

function _modem_origin(tfs)
    located = [tf for tf in tfs if isfinite(tf.latitude) && isfinite(tf.longitude)]
    isempty(located) && return (0.0, 0.0)
    lat = [tf.latitude for tf in located]
    lon = [tf.longitude for tf in located]
    return ((minimum(lat) + maximum(lat)) / 2, (minimum(lon) + maximum(lon)) / 2)
end

"""
    write_modem(path, tfs; components = :full, units = :field, sign = :plus, origin = nothing,
                elevation = false, description = "") -> String

Write [`TransferFunction`](@ref)s (one or a vector) to a ModEM data file.

- `components` -- `:full` (the four elements of `Z` and the tipper),
  `:offdiagonal` (only `Zxy`, `Zyx` and the tipper) or `:impedance` (no
  tipper).
- `units` -- `:field` writes `[mV/km]/[nT]`. `:ohm` writes `Z` in Ω (`[V/m]/[A/m]`).
- `sign` -- `:plus` writes `exp(+i\\omega t)`, the convention of the estimate.
  `:minus` writes the complex conjugate with `exp(-i\\omega t)`.
- `origin` -- `(lat, lon)` of the grid origin. The default is the centre of the
  sites.
- `elevation` -- `true` writes `z = -elevation`. With `false` (the default),
  `z = 0`, and the mesh puts the sites on its surface.

The error of an element is its estimated standard deviation along one axis,
`√(var/2)`, with no floor: set the error floors of an inversion in its own
set-up. The file leaves out the elements without an estimate or with an error
that is not finite. The function returns `path`.
"""
function write_modem(path::AbstractString, tfs; components::Symbol = :full, units::Symbol = :field,
                     sign::Symbol = :plus, origin = nothing, elevation::Bool = false,
                     description::AbstractString = "")
    tfs = tfs isa TransferFunction ? [tfs] : collect(TransferFunction, tfs)
    isempty(tfs) && error("No transfer functions to write")
    components in (:full, :offdiagonal, :impedance) || error("components must be :full, :offdiagonal or :impedance")
    units in (:field, :ohm) || error("units must be :field or :ohm")
    sign in (:plus, :minus) || error("sign must be :plus or :minus")
    lat0, lon0 = origin === nothing ? _modem_origin(tfs) : origin
    # (mV/km)/nT = 1e-6 V/m / (1e-9 T) = 1e3 V/m/T; Ω = (V/m)/(A/m) = μ0 (V/m)/T
    zscale = units === :ohm ? 1.0e3 * _MU0 : 1.0
    zunits = units === :ohm ? "[V/m]/[A/m]" : "[mV/km]/[nT]"
    conv(z) = sign === :minus ? conj(z) : z
    zlines, tlines = String[], String[]
    zper, tper = Set{Float64}(), Set{Float64}()
    zsites, tsites = Set{String}(), Set{String}()
    for tf in tfs
        x, y = isfinite(tf.latitude) ? _local_xy(tf.latitude, tf.longitude, lat0, lon0) : (0.0, 0.0)
        z = elevation && isfinite(tf.elevation) ? -tf.elevation : 0.0
        lat = isfinite(tf.latitude) ? tf.latitude : 0.0
        lon = isfinite(tf.longitude) ? tf.longitude : 0.0
        head(T) = @sprintf("%13.6E %-12s %10.5f %10.5f %13.3f %13.3f %11.3f", T, tf.site, lat, lon, x, y, z)
        for k in eachindex(tf.periods)
            T = tf.periods[k]
            for (e, (i, j)) in enumerate(((1, 1), (1, 2), (2, 1), (2, 2)))
                components === :offdiagonal && i == j && continue
                v = tf.Z[i, j, k]
                err = sqrt(tf.Z_var[i, j, k] / 2)
                isfinite(v) && isfinite(err) || continue
                err *= zscale
                v = conv(v) * zscale
                push!(zlines, head(T) * @sprintf(" %-4s %14.6E %14.6E %14.6E", _Z_NAMES[e], real(v), imag(v), err))
                push!(zper, T)
                push!(zsites, tf.site)
            end
            components === :impedance && continue
            for j in 1:2
                v = tf.T[j, k]
                err = sqrt(tf.T_var[j, k] / 2)
                isfinite(v) && isfinite(err) || continue
                v = conv(v)
                push!(tlines, head(T) * @sprintf(" %-4s %14.6E %14.6E %14.6E", _T_NAMES[j], real(v), imag(v), err))
                push!(tper, T)
                push!(tsites, tf.site)
            end
        end
    end
    isempty(zlines) && isempty(tlines) && error("The transfer functions have no finite estimates to write")
    expo = sign === :minus ? "> exp(-i\\omega t)" : "> exp(+i\\omega t)"
    desc = isempty(description) ? "Timekeepers.jl transfer functions: " * join(unique(tf.site for tf in tfs), ", ") : description
    mkpath(dirname(abspath(path)))
    open(path, "w") do io
        if !isempty(zlines)
            println(io, "# ", desc)
            println(io, _MODEM_COLUMNS)
            println(io, "> Full_Impedance")
            println(io, expo)
            println(io, "> ", zunits)
            println(io, "> 0.00")
            @printf(io, "> %.6f %.6f\n", lat0, lon0)
            @printf(io, "> %d %d\n", length(zper), length(zsites))
            foreach(l -> println(io, l), zlines)
        end
        if !isempty(tlines)
            println(io, "# ", desc)
            println(io, _MODEM_COLUMNS)
            println(io, "> Full_Vertical_Components")
            println(io, expo)
            println(io, "> []")
            println(io, "> 0.00")
            @printf(io, "> %.6f %.6f\n", lat0, lon0)
            @printf(io, "> %d %d\n", length(tper), length(tsites))
            foreach(l -> println(io, l), tlines)
        end
    end
    return path
end

"""
    read_modem(path) -> Vector{TransferFunction}

Read the `Full_Impedance`, `Off_Diagonal_Impedance` and
`Full_Vertical_Components` blocks of a ModEM data file. The result has one
[`TransferFunction`](@ref) for each site, in the order of the file, in
(mV/km)/nT with the convention `exp(+iωt)`. The variance of each element is
`2·error²`. The fields that a data file does not hold (coherences, counts)
are `NaN` or zero.
"""
function read_modem(path::AbstractString)
    order = String[]
    pos = Dict{String, NTuple{3, Float64}}()
    vals = Dict{String, Dict{Tuple{Float64, String}, Tuple{ComplexF64, Float64}}}()
    kind = :none
    header = 0
    scale = 1.0
    minus = false
    for line in eachline(path)
        s = strip(line)
        isempty(s) && continue
        startswith(s, "#") && continue
        if startswith(s, ">")
            body = strip(s[2:end])
            if body in ("Full_Impedance", "Off_Diagonal_Impedance", "Full_Vertical_Components")
                kind = body == "Full_Vertical_Components" ? :tipper : :impedance
                header = 1
            elseif header == 1
                minus = occursin("-", body)
                header = 2
            elseif header == 2
                scale = occursin("[V/m]/[A/m]", body) || lowercase(body) == "ohm" ? 1 / (1.0e3 * _MU0) :
                        occursin("[V/m]/[T]", body) ? 1.0e-3 : 1.0
                header = 3
            else
                header += 1
            end
            continue
        end
        parts = split(s)
        length(parts) >= 11 || continue
        T = parse(Float64, parts[1])
        code = String(parts[2])
        lat, lon = parse(Float64, parts[3]), parse(Float64, parts[4])
        z = parse(Float64, parts[7])
        comp = uppercase(parts[8])
        v = complex(parse(Float64, parts[9]), parse(Float64, parts[10]))
        err = parse(Float64, parts[11])
        minus && (v = conj(v))
        if kind === :impedance
            v *= scale
            err *= scale
        end
        haskey(vals, code) || (push!(order, code); vals[code] = Dict{Tuple{Float64, String}, Tuple{ComplexF64, Float64}}())
        pos[code] = (lat, lon, -z)
        vals[code][(T, comp)] = (v, err)
    end
    out = TransferFunction[]
    for code in order
        d = vals[code]
        periods = sort!(unique(first.(collect(keys(d)))))
        n = length(periods)
        Z = fill(complex(NaN, NaN), 2, 2, n)
        ZV = fill(NaN, 2, 2, n)
        Tz = fill(complex(NaN, NaN), 2, n)
        TV = fill(NaN, 2, n)
        for (k, T) in enumerate(periods)
            for (e, (i, j)) in enumerate(((1, 1), (1, 2), (2, 1), (2, 2)))
                x = get(d, (T, _Z_NAMES[e]), nothing)
                x === nothing || (Z[i, j, k] = x[1]; ZV[i, j, k] = 2 * x[2]^2)
            end
            for j in 1:2
                x = get(d, (T, _T_NAMES[j]), nothing)
                x === nothing || (Tz[j, k] = x[1]; TV[j, k] = 2 * x[2]^2)
            end
        end
        lat, lon, elev = pos[code]
        push!(out, TransferFunction(code, :single, "", String[], lat, lon, elev, periods, Z, ZV, Tz, TV,
                                    fill(NaN, 3, n), fill(NaN, 2, n), zeros(Int, n), zeros(Int, n),
                                    fill(NaN, n), Dict{Symbol, Any}(:source => abspath(path))))
    end
    return out
end
