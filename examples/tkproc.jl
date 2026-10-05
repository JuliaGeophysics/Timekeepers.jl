# tkproc.jl - starts the transfer function window.
# Author: @pankajkmishra
#
# This script opens TKProc on one site. TKProc scans the survey around the site
# (the directory above it) and loads its base sites (they recorded with it,
# near it) and its remote sites (they recorded with it, far from it) into its
# menus, with the best set-up selected. If you selected the base and remote
# sites of each site in TKDash and exported them (reference_plan.txt in the
# survey directory), TKProc uses that plan: the menus hold the sites of the
# plan, one base and one remote site at a time. The TKDash plan switch turns
# it off.
#
# Set the processing options in the second row and press Process. The status
# line shows each step, and the line above it tells if a channel looks
# flipped. Load Site opens the next site. Export writes the site as an EDI
# file, a ModEM data file, a PNG of the impedance (off-diagonal, or the full
# tensor with the switch) and the tipper, and <site>.md, the record of the
# processing, into a directory that you select. Give the site directory as the first argument. If you give no
# argument, the script uses site004 of the BRGM Alsace survey next to this
# repository. If that survey is not there, the window opens empty
#
#   julia --project=. -t auto examples/tkproc.jl /path/to/survey/site004
#
# When you close the window, the script prints the estimates

using Timekeepers

const DEFAULT_SITE = normpath(joinpath(@__DIR__, "..", "..", "MetronixATS", "site004"))

site = !isempty(ARGS) ? ARGS[1] : isdir(DEFAULT_SITE) ? DEFAULT_SITE : nothing
proc = site === nothing ? run_tkproc() : run_tkproc(site)

for (name, tf) in sort!(collect(proc.results); by = first)
    println(rpad(name, 16), tf)
end
