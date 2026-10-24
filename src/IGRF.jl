# IGRF.jl - the main geomagnetic field from IGRF-14, for the declination of a survey.
# Author: @pankajkmishra
#
# A fluxgate that was turned to magnetic north by hand needs the declination
# to give its axes in geographic coordinates (the x axis from the static field).
# This file computes the declination from the 14th generation International
# Geomagnetic Reference Field (Alken et al., 2025), with the coefficients of
# igrf14coeffs.txt beside this file (NOAA NCEI / IAGA V-MOD), degree 13. The model holds
# from 1900 to 2030: linear between the 5-year models, and with the secular
# variation after 2025. The positions are geodetic (WGS84).

const _IGRF_FILE = joinpath(@__DIR__, "igrf14coeffs.txt")
const _IGRF_NMAX = 13
const _IGRF_A = 6371.2

struct _IGRFModel
    epochs::Vector{Float64}
    g::Array{Float64, 3}
    h::Array{Float64, 3}
    gsv::Matrix{Float64}
    hsv::Matrix{Float64}
end

const _IGRF = Ref{Union{Nothing, _IGRFModel}}(nothing)

function _igrf_model()
    _IGRF[] === nothing || return _IGRF[]
    lines = readlines(_IGRF_FILE)
    head = findfirst(l -> startswith(l, "g/h"), lines)
    head === nothing && error("$(_IGRF_FILE) has no 'g/h n m' header line")
    tok = split(lines[head])
    epochs = parse.(Float64, tok[4:(end - 1)])
    ne = length(epochs)
    g = zeros(_IGRF_NMAX, _IGRF_NMAX + 1, ne); h = zeros(size(g))
    gsv = zeros(_IGRF_NMAX, _IGRF_NMAX + 1); hsv = zeros(size(gsv))
    for l in lines[(head + 1):end]
        t = split(l)
        length(t) == ne + 4 || continue
        n, m = parse(Int, t[2]), parse(Int, t[3])
        v = parse.(Float64, t[4:end])
        if t[1] == "g"
            g[n, m + 1, :] .= v[1:ne]; gsv[n, m + 1] = v[end]
        elseif t[1] == "h"
            h[n, m + 1, :] .= v[1:ne]; hsv[n, m + 1] = v[end]
        end
    end
    _IGRF[] = _IGRFModel(epochs, g, h, gsv, hsv)
    return _IGRF[]
end

_decimal_year(t::DateTime) = (y = year(t); y + Dates.value(t - DateTime(y)) / Dates.value(DateTime(y + 1) - DateTime(y)))
_decimal_year(t::Date) = _decimal_year(DateTime(t))
_decimal_year(t::Real) = Float64(t)

function _igrf_coefficients(y::Float64)
    M = _igrf_model()
    e = M.epochs
    (e[1] <= y <= e[end] + 5) || error(@sprintf("IGRF-14 holds from %.0f to %.0f; got %.2f", e[1], e[end] + 5, y))
    if y >= e[end]
        dt = y - e[end]
        return M.g[:, :, end] .+ dt .* M.gsv, M.h[:, :, end] .+ dt .* M.hsv
    end
    k = searchsortedlast(e, y)
    w = (y - e[k]) / (e[k + 1] - e[k])
    return (1 - w) .* M.g[:, :, k] .+ w .* M.g[:, :, k + 1], (1 - w) .* M.h[:, :, k] .+ w .* M.h[:, :, k + 1]
end

"""
    igrf_field(latitude, longitude, t; altitude_km = 0.0) -> NamedTuple

The main geomagnetic field of IGRF-14 at a geodetic (WGS84) position and time
`t` (a `DateTime`, a `Date` or a decimal year, 1900-2030). The result has the
north, east and down components (nT), `declination` (degrees east of true
north), `inclination` (degrees, down positive) and `intensity` (nT).

See also [`igrf_declination`](@ref), [`survey_declination`](@ref).
"""
function igrf_field(latitude::Real, longitude::Real, t; altitude_km::Real = 0.0)
    (-90 <= latitude <= 90) || error("latitude $latitude is not in -90…90")
    g, h = _igrf_coefficients(_decimal_year(t))
    a2, b2 = 6378.137^2, 6356.752314245^2
    colat = 90.0 - latitude
    ct, st = cosd(colat), sind(colat)
    one, two = a2 * st^2, b2 * ct^2
    three = one + two
    rho = sqrt(three)
    r = sqrt(altitude_km * (altitude_km + 2rho) + (a2 * one + b2 * two) / three)
    cd = (altitude_km + rho) / r
    sd = (a2 - b2) / rho * ct * st / r
    cθ = ct * cd - st * sd
    sθ = st * cd + ct * sd
    sθ = max(sθ, 1.0e-10)
    φ = deg2rad(longitude)
    N = _IGRF_NMAX
    P = zeros(N + 1, N + 1); dP = zeros(N + 1, N + 1)
    P[1, 1] = 1.0
    for n in 1:N, m in 0:n
        if m == n
            k = n == 1 ? 1.0 : sqrt((2n - 1) / 2n)
            P[n + 1, n + 1] = k * sθ * P[n, n]
            dP[n + 1, n + 1] = k * (sθ * dP[n, n] + cθ * P[n, n])
        else
            c1 = sqrt((n - 1)^2 - m^2)
            p2 = n >= 2 ? P[n - 1, m + 1] : 0.0
            dp2 = n >= 2 ? dP[n - 1, m + 1] : 0.0
            d = sqrt(n^2 - m^2)
            P[n + 1, m + 1] = ((2n - 1) * cθ * P[n, m + 1] - c1 * p2) / d
            dP[n + 1, m + 1] = ((2n - 1) * (cθ * dP[n, m + 1] - sθ * P[n, m + 1]) - c1 * dp2) / d
        end
    end
    Br = Bt = Bp = 0.0
    for n in 1:N
        f = (_IGRF_A / r)^(n + 2)
        for m in 0:n
            cm, sm = cos(m * φ), sin(m * φ)
            gnm, hnm = g[n, m + 1], h[n, m + 1]
            Br += (n + 1) * f * (gnm * cm + hnm * sm) * P[n + 1, m + 1]
            Bt -= f * (gnm * cm + hnm * sm) * dP[n + 1, m + 1]
            Bp += f * m * (gnm * sm - hnm * cm) * P[n + 1, m + 1] / sθ
        end
    end
    X, Y, Z = -Bt, Bp, -Br
    X, Z = X * cd + Z * sd, Z * cd - X * sd
    H = hypot(X, Y)
    return (north = X, east = Y, down = Z, declination = atand(Y, X), inclination = atand(Z, H),
            intensity = hypot(H, Z))
end

"""
    igrf_declination(latitude, longitude, t; altitude_km = 0.0) -> Float64

The declination of IGRF-14 (degrees east of true north) at a geodetic
position and time. Refer to [`igrf_field`](@ref).
"""
igrf_declination(latitude::Real, longitude::Real, t; altitude_km::Real = 0.0) =
    igrf_field(latitude, longitude, t; altitude_km).declination

"""
    survey_declination(survey) -> NamedTuple

The declination that IGRF-14 gives for a survey: at the median position of its
sites with a position, at the middle of its recordings. The result has
`declination` (degrees east), `latitude`, `longitude` and `time`, or
`declination = NaN` if no site has a position.
"""
function survey_declination(s::Survey)
    pos = [(x.latitude, x.longitude) for x in s.sites if isfinite(x.latitude) && isfinite(x.longitude)]
    times = [t for x in s.sites for r in x.runs for t in (r.start, r.stop) if year(r.start) >= 1900]
    isempty(pos) && return (declination = NaN, latitude = NaN, longitude = NaN, time = nothing)
    lat, lon = median(first.(pos)), median(last.(pos))
    t = isempty(times) ? now(UTC) : sort(times)[(length(times) + 1) ÷ 2]
    d = try
        igrf_declination(lat, lon, t)
    catch
        NaN
    end
    return (declination = d, latitude = lat, longitude = lon, time = t)
end
