# Polarity.jl - a check of the channels of a site: sign, layout and labels.
# Author: @pankajkmishra
#
# A sensor that is connected the other way (or an electrode pair laid out the
# other way round) turns the sign of its channel. Two sensors that point the
# same way record the same signal. Two channels with swapped labels exchange
# their roles. The check finds these from the estimate, with three kinds of
# evidence:
#
# 1. The impedance. Zxy = Ex/Hy and Zyx = Ey/Hx. With the time dependence
#    exp(+iωt), over a 1D or 2D earth the phase of Zxy is near +45° and the
#    phase of Zyx near -135°. A Zxy with the opposite sign tells that Ex or Hy
#    is flipped, and the impedance alone cannot tell which: the two give the
#    same Z. A flipped channel moves the phase by 180° at every period. 3D
#    structure moves it too, but mostly at the long periods, which sense the
#    deep ground. Thus, the short periods decide: the first quarter of the
#    clear periods. Three quarters of them must agree, so a few points that
#    noise moves out of the quadrant do not change the verdict. The longer
#    periods only add a note:
#    out of their quadrant, they show 3D structure, not a flipped channel. A
#    sign that is wrong at the short periods and right at the long ones is
#    not a flip either: it points to the calibration or the filters at the
#    high frequencies, or to shallow 3D structure.
# 2. The magnetic field at other sites: the remote sites of the estimate and
#    the witness sites (often base sites, which take no part in the
#    estimate). The field is nearly the same over tens of kilometres, so Hx
#    of a site follows Hx of another site with the same sign. This is the
#    second witness that the impedance needs: if H agrees with the other sites
#    and Z has the wrong sign, E is flipped. With two or more other sites, the
#    site that disagrees with the others is the flipped one. The electric field has no such witness, because it depends on the
#    ground near each site.
# 3. The two channels of each pair. Natural fields change their direction, so
#    over a long record Hx and Hy (and Ex and Ey) are far from fully coherent.
#    A coherence near 1 at almost every period, with an impedance whose rows
#    are in proportion, tells that the two sensors are parallel. An impedance
#    with Zxx and Zyy larger than Zxy and Zyx tells that the labels x and y
#    are swapped (or that the layout is turned about 90°)

const POL_MIN_PERIODS = 3
const POL_SHORT_SHARE = 0.25                # the short periods: this share of the clear periods
const POL_SHARE = 0.75                      # the share of the periods for a verdict
const POL_ZCOH = 0.5
const POL_ZREL = 0.3
const POL_HCOH = 0.6
const POL_PARALLEL_COH = 0.95
const POL_PARALLEL_DET = 0.15
const POL_SWAP_RATIO = 2.0

# The verdict of a test that is true or false at each period: :yes, :no or
# :unclear, with the count of true periods and of all periods
function _pol_share(flags)
    n = length(flags)
    n < POL_MIN_PERIODS && return (:unclear, 0, n)
    k = count(flags)
    return (k >= POL_SHARE * n ? :yes : k <= (1 - POL_SHARE) * n ? :no : :unclear, k, n)
end

const QUAD_XY = "the first quadrant (0° to 90°)"
const QUAD_YX = "the third quadrant (-180° to -90°)"

_pol_deg(a) = isfinite(a) ? @sprintf("%+.0f°", a) : "?"

_pol_off(a, expected) = abs(rem(a - expected, 360, RoundNearest)) > 90

# |det A| against the size of its columns: near 0 when the two columns point
# the same way
# "a", "a and b", "a, b and c"
_pol_list(v) = length(v) <= 1 ? join(v) : join(v[1:(end - 1)], ", ") * " and " * v[end]

_pol_det(A) = (a = abs(A[1, 1] * A[2, 2] - A[1, 2] * A[2, 1]);
               b = norm(A[:, 1]) * norm(A[:, 2]); b > 0 ? a / b : NaN)

"""
    check_polarity(tf::TransferFunction) -> NamedTuple

Check the channels of the site of `tf` for a flipped sign, for parallel
sensors and for swapped labels. The result has:
- `zxy`, `zyx` -- `:normal`, `:flipped` (the opposite sign) or `:unclear`,
  from the first quarter of the clear periods (the shortest), three quarters
  of which must agree. Longer periods out of their quadrant are 3D structure:
  the message notes them;
With `magnetic = false`, the check uses only the impedance and the pairs of
channels: it tells that something is flipped and which channels can be at
fault. With `magnetic = true` (the default), it also compares H of the inputs
with the other sites that the estimate holds (the remote sites, and the
witnesses of [`estimate_tf`](@ref)) to tell which.

- `h` -- for each other site (the remote sites and the witnesses of
  [`estimate_tf`](@ref)): `(site, hx, hy)`, each `:agree`, `:flipped` or
  `:unclear`;
- `parallel_e`, `parallel_h` -- `true` when Ex and Ey (or Hx and Hy of the
  inputs) record the same signal;
- `swapped` -- `true` when Zxx and Zyy are larger than Zxy and Zyx;
- `flipped` -- the channels that the check finds sign-flipped, e.g.
  `["site004 Ex", "site004 Ey"]`;
- `ok` -- `true` when the check finds no problem;
- `message` -- the result in plain words. When the evidence cannot tell which
  channel is at fault, it names the candidates and what would tell them
  apart.

The sign test expects the time dependence `exp(+iωt)` of
[`estimate_tf`](@ref) and shallow ground close to 1D or 2D. Thus, read a
verdict on one site together with its neighbours. The check reports; it does
not change the estimate.
"""
function check_polarity(tf::TransferFunction; magnetic::Bool = true)
    site = tf.site
    T = tf.periods
    n = length(T)
    mag = get(tf.metadata, :magnetic, nothing)
    input = mag === nothing ? (isempty(tf.base) ? site : tf.base) : mag.input
    pcoh = get(tf.metadata, :pair_coherence, fill(NaN, 2, n))
    clear(i, j, c, k) = isfinite(tf.Z[i, j, k]) && isfinite(tf.Z_var[i, j, k]) &&
                        sqrt(tf.Z_var[i, j, k] / 2) / abs(tf.Z[i, j, k]) < POL_ZREL &&
                        (!isfinite(tf.coherence[c, k]) || tf.coherence[c, k] >= POL_ZCOH)
    # the periods for the layout tests: all four elements, with the outputs
    # coherent with the inputs. Not the error of Zxy and Zyx: with swapped
    # labels they are the small elements
    full = [k for k in 1:n if all(isfinite, tf.Z[:, :, k]) &&
            all(c -> !isfinite(tf.coherence[c, k]) || tf.coherence[c, k] >= POL_ZCOH, 1:2)]
    notes = String[]
    flipped = String[]
    ratio(k, n) = "$k of $n periods"

    #---- the layout: parallel sensors and swapped labels
    pe = _pol_share([pcoh[1, k] > POL_PARALLEL_COH && _pol_det(transpose(tf.Z[:, :, k])) < POL_PARALLEL_DET
                     for k in 1:n if isfinite(pcoh[1, k]) && all(isfinite, tf.Z[:, :, k])])
    ph = _pol_share([pcoh[2, k] > POL_PARALLEL_COH for k in 1:n if isfinite(pcoh[2, k])])
    sw = _pol_share([abs(tf.Z[1, 1, k]) + abs(tf.Z[2, 2, k]) > POL_SWAP_RATIO * (abs(tf.Z[1, 2, k]) + abs(tf.Z[2, 1, k]))
                     for k in full])
    parallel_e, parallel_h, swapped = pe[1] === :yes, ph[1] === :yes, sw[1] === :yes
    parallel_e && push!(notes, "Ex and Ey of $site record the same signal ($(ratio(pe[2], pe[3]))): the two dipoles " *
                               "are parallel, or one channel is a copy of the other. Only one direction of the " *
                               "electric field was measured, so only one row of Z is real.")
    parallel_h && push!(notes, "Hx and Hy of $input record the same signal ($(ratio(ph[2], ph[3]))): the coils are " *
                               "parallel, or one channel is a copy of the other. Z from these inputs is not " *
                               "reliable: use a base site with good Hx, Hy for the inputs.")

    #---- the magnetic field against each other site: sign and labels
    hrefs = NamedTuple[]
    hswap = Dict{String, Bool}()
    if magnetic && mag !== nothing && !isempty(mag.references)
        for (m, ref) in enumerate(mag.references)
            ok(j, k) = isfinite(mag.coherence[j, m, k]) && mag.coherence[j, m, k] >= POL_HCOH &&
                       all(isfinite, mag.B[:, :, m, k])
            v = map(1:2) do j
                verdict, _, _ = _pol_share([_pol_off(rad2deg(angle(mag.B[j, j, m, k])), 0.0) for k in 1:n if ok(j, k)])
                verdict === :yes ? :flipped : verdict === :no ? :agree : :unclear
            end
            s, _, _ = _pol_share([abs(mag.B[1, 2, m, k]) + abs(mag.B[2, 1, m, k]) >
                                  POL_SWAP_RATIO * (abs(mag.B[1, 1, m, k]) + abs(mag.B[2, 2, m, k]))
                                  for k in 1:n if ok(1, k) && ok(2, k)])
            hswap[ref] = s === :yes
            push!(hrefs, (site = ref, hx = hswap[ref] ? :unclear : v[1], hy = hswap[ref] ? :unclear : v[2]))
        end
    end

    # Hx, Hy labels swapped between the inputs and the other sites
    swapped_refs = [r for (r, s) in hswap if s]
    if !isempty(swapped_refs)
        if length(swapped_refs) == length(hrefs) && length(hrefs) >= 2
            push!(notes, "Hx and Hy of $input are swapped: each other site ($(join(swapped_refs, ", "))) sees " *
                         "its Hx as their Hy.")
        elseif length(swapped_refs) == length(hrefs)
            push!(notes, "Hx and Hy are swapped between $input and $(only(swapped_refs)); a third site " *
                         "would tell which site has the swapped labels.")
        else
            push!(notes, "Hx and Hy of $(join(swapped_refs, ", ")) are swapped: they disagree with $input and " *
                         "the other sites.")
        end
    end

    # the sign of each H of the inputs: :ok, :flipped or :unknown. The notes
    # name Hx and Hy together when the two agree
    hstate = Dict{Symbol, Symbol}()
    findings = Dict{Symbol, Any}()
    for key in (:hx, :hy)
        agree = [h.site for h in hrefs if getfield(h, key) === :agree]
        disagree = [h.site for h in hrefs if getfield(h, key) === :flipped]
        if isempty(agree) && isempty(disagree)
            hstate[key], findings[key] = :unknown, nothing
        elseif isempty(disagree)
            hstate[key], findings[key] = :ok, nothing
        elseif isempty(agree) && length(disagree) >= 2
            hstate[key], findings[key] = :flipped, (:input, disagree, agree)
        elseif isempty(agree)
            hstate[key], findings[key] = :unknown, (:pair, disagree, agree)
        else
            hstate[key], findings[key] = :ok, (:refs, disagree, agree)
        end
    end
    comps = findings[:hx] == findings[:hy] && findings[:hx] !== nothing ? [("Hx and Hy", findings[:hx], ["Hx", "Hy"])] :
            [(c, findings[k], [c]) for (c, k) in (("Hx", :hx), ("Hy", :hy)) if findings[k] !== nothing]
    for (label, (kind, disagree, agree), cs) in comps
        they = length(cs) == 2 ? "they go" : "it goes"
        if kind === :input
            append!(flipped, ["$input $c" for c in cs])
            push!(notes, "$input has its $label sign-flipped: $they opposite to every other site " *
                         "($(_pol_list(disagree))).")
        elseif kind === :pair
            push!(notes, "$input and $(only(disagree)) have $label of opposite sign: one of the two sites is " *
                         "flipped. A third site (remote or base) would tell which.")
        else
            for r in disagree
                append!(flipped, ["$r $c" for c in cs])
                push!(notes, "$r has its $label sign-flipped: $they opposite to $(_pol_list(vcat(input, agree))).")
            end
        end
    end

    #---- the sign of the impedance: Zxy = Ex/Hy, Zyx = Ey/Hx. The short
    # periods decide, the long ones add a note
    function sign_test(i, j, c, expected)
        ks = [k for k in 1:n if clear(i, j, c, k)]
        off = Dict(k => _pol_off(rad2deg(angle(tf.Z[i, j, k])), expected) for k in ks)
        isempty(ks) && return (verdict = (:unclear, 0, 0), state = :unclear, upto = NaN, long = (0, 0), from = NaN,
                               phase = NaN)
        short = ks[1:min(length(ks), max(POL_MIN_PERIODS, ceil(Int, POL_SHORT_SHARE * length(ks))))]
        long = setdiff(ks, short)
        v = _pol_share([off[k] for k in short])
        nlong_off = count(k -> off[k], long)
        from = nlong_off == 0 ? NaN : T[first(k for k in long if off[k])]
        st = v[1] === :yes ? :flipped : v[1] === :no ? :normal : :unclear
        # a flip turns every period. Wrong at the short periods but right at
        # most of the long ones is something else
        if st === :flipped && length(long) >= POL_MIN_PERIODS && nlong_off <= (1 - POL_SHARE) * length(long)
            st = :conflict
        end
        # the mean direction of the phase at the short periods
        phase = rad2deg(angle(sum(tf.Z[i, j, k] / abs(tf.Z[i, j, k]) for k in short)))
        return (verdict = v, state = st, upto = T[last(short)], long = (nlong_off, length(long)), from = from,
                phase = phase)
    end
    zxy = sign_test(1, 2, 1, 45.0)
    zyx = sign_test(2, 1, 2, -135.0)
    zs = (xy = zxy.state === :conflict ? :unclear : zxy.state, yx = zyx.state === :conflict ? :unclear : zyx.state)
    period(t) = t >= 1 ? @sprintf("%.3g s", t) : @sprintf("%.2g s", t)
    info = String[]
    zused = false
    if parallel_e || parallel_h || swapped || !isempty(swapped_refs)
        if swapped
            who = any(values(hswap)) ? "" :
                  isempty(hrefs) ? " (Ex and Ey, or Hx and Hy; a remote or base site would tell which)" :
                  " (Ex and Ey: Hx, Hy keep their labels against the other sites)"
            push!(notes, "Zxx and Zyy are larger than Zxy and Zyx ($(ratio(sw[2], sw[3]))): the x and y labels " *
                         "of a pair are swapped$who, or the layout is turned about 90°.")
        end
        push!(notes, "The sign test of Z waits until the layout is right.")
    else
        e_flip, e_cancel, both_amb = String[], Tuple{String, String, String}[], String[]
        for (z, e, h, key, label, expect) in ((zxy, "Ex", "Hy", :hy, "Zxy", "+45°"), (zyx, "Ey", "Hx", :hx, "Zyx", "-135°"))
            st = z.state
            if st === :conflict
                zused = true
                push!(notes, "$label has the wrong sign at the short periods (up to $(period(z.upto))) but the " *
                             "right sign at most long periods. A flipped channel would turn every period, so " *
                             "this is more likely the calibration or the filters at the high frequencies, or " *
                             "shallow 3D structure.")
            elseif st === :flipped && hstate[key] === :ok
                push!(e_flip, e)
            elseif st === :flipped && hstate[key] === :flipped
                zused = true
                push!(notes, "The flipped $h of $input explains the sign of $label: $e of $site is fine.")
            elseif st === :flipped
                push!(both_amb, "$e of $site or $h of $input")
            elseif st === :normal && hstate[key] === :flipped
                # E and H flipped together: their effects on Z cancel
                push!(e_cancel, (e, h, label))
            elseif st === :unclear && z.verdict[3] >= POL_MIN_PERIODS
                zused = true
                push!(notes, "$label is near $expect at only $(z.verdict[3] - z.verdict[2]) of the " *
                             "$(z.verdict[3]) short periods (up to $(period(z.upto))): no verdict on its sign " *
                             "(noise or shallow 3D structure can do this).")
            end
            # out of the quadrant at the long periods only: the deep ground
            if st === :normal && z.long[1] >= max(1, (1 - POL_SHARE) * z.long[2])
                push!(info, "$label is in its quadrant at the short periods and leaves it at $(z.long[1]) of " *
                            "$(z.long[2]) longer periods (from $(period(z.from))): 3D structure, not a flipped channel.")
            end
        end
        if !isempty(e_flip)
            zused = true
            append!(flipped, ["$site $e" for e in e_flip])
            upto = period(max(zxy.upto, zyx.upto))
            seen = join((e == "Ex" ? "Zxy is near $(_pol_deg(zxy.phase)), not in $QUAD_XY" :
                         "Zyx is near $(_pol_deg(zyx.phase)), not in $QUAD_YX" for e in e_flip), " and ")
            push!(notes, "$site has its $(join(e_flip, " and ")) sign-flipped: at the short periods (up to " *
                         "$upto) $seen, while H agrees with the other sites.")
        end
        for (e, h, label) in e_cancel
            zused = true
            push!(flipped, "$site $e")
            push!(notes, "$site has its $e sign-flipped too: its $h is flipped, yet $label has the right sign, " *
                         "so $e must be flipped as well.")
        end
        if !isempty(both_amb)
            zused = true
            how = !magnetic ? "Comparing H with other sites (FlipCheck) tells which." :
                  isempty(hrefs) ? "No other site recorded H at the same time to tell which." :
                  "A third site (remote or base) would tell which."
            upto = period(max(zxy.upto, zyx.upto))
            if length(both_amb) == 2
                push!(notes, "Something is flipped. Zxy should lie in $QUAD_XY and Zyx in $QUAD_YX, but at the " *
                             "short periods (up to $upto) Zxy is near $(_pol_deg(zxy.phase)) and Zyx near " *
                             "$(_pol_deg(zyx.phase)): both have the opposite sign. Either Ex and Ey of $site are " *
                             "flipped, or Hx and Hy of $input are; Z alone cannot tell E from H. $how")
            else
                xy = startswith(only(both_amb), "Ex")
                push!(notes, "Something is flipped. $(xy ? "Zxy" : "Zyx") should lie in $(xy ? QUAD_XY : QUAD_YX), " *
                             "but at the short periods (up to $upto) it is near " *
                             "$(_pol_deg(xy ? zxy.phase : zyx.phase)): the opposite sign. Either $(only(both_amb)) is " *
                             "flipped; Z alone cannot tell E from H. $how")
            end
        end
    end
    zused && push!(notes, "(The sign test of Z assumes the shallow ground is close to 1D or 2D.)")

    nz = zxy.verdict[3] + zyx.verdict[3]
    ok = isempty(notes) && zs.xy === :normal && zs.yx === :normal
    message = ok ? "No sign of a flipped channel: at the short periods Zxy is near $(_pol_deg(zxy.phase)) and Zyx " *
                   "near $(_pol_deg(zyx.phase)), in their quadrants; the channels of each pair record different " *
                   "signals" * (isempty(hrefs) ? "." : "; and Hx, Hy agree with the other sites.") :
              isempty(notes) ? (nz < 2POL_MIN_PERIODS ? "Too few clear periods to check the channels." :
                                "No clear verdict on the channels.") :
              join(notes, " ")
    isempty(info) || (message *= " " * join(info, " "))
    return (zxy = zs.xy, zyx = zs.yx, h = hrefs, parallel_e, parallel_h, swapped, flipped = unique(flipped),
            ok, message)
end

"""
    flip_check!(survey, tf::TransferFunction; witnesses = 2, progress = nothing) -> NamedTuple

Find which channel is flipped: compare Hx, Hy of the inputs of `tf` with up
to `witnesses` base sites of the site in `survey` (the sites that recorded
with it, near it, longest overlap first; not the base or remote sites of the
estimate). The comparison takes no part in the estimate. It adds the
witnesses to the evidence of `tf`, next to its remote sites, and returns
[`check_polarity`](@ref)`(tf)`. The EDI and the record of `tf` then report
the full check.
"""
function flip_check!(s::Survey, tf::TransferFunction; witnesses::Integer = 2, progress = nothing)
    say = progress === nothing ? (_ -> nothing) : progress
    input = isempty(tf.base) ? tf.site : tf.base
    used = Set(vcat(tf.site, tf.base, tf.remote, get(tf.metadata, :witnesses, String[])))
    refs = site_references(s, s[tf.site])
    names = first([c.site for c in refs.base if !(c.site in used) && has_magnetic(s[c.site])], witnesses)
    if !isempty(names)
        rates = collect(keys(get(tf.metadata, :rates, Dict(NaN => nothing))))
        _add_witnesses!(tf, s[input], Any[s[n] for n in names], rates, tf.metadata[:options]; progress = say)
    end
    tf.metadata[:flipcheck] = true
    return check_polarity(tf)
end
