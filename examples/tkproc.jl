# tkproc.jl - starts the transfer function window.
# Author: @pankajkmishra
#
# This script opens TKProc on a survey. TKProc scans the survey directory and
# lists its sites with electric channels in the Site menu. When you select a
# site, TKProc loads its base sites (they recorded with it, near it) and its
# remote sites (they recorded with it, far from it) into their menus, with the
# best set-up selected. If you selected the base and remote sites of each site
# in TKDash and exported them (reference_plan.txt in the survey directory),
# TKProc uses that plan: the menus hold the sites of the plan, one base and
# one remote site at a time. Without a plan, TKProc makes the lists itself
#
# Set the processing options in the second row and press Process. The status
# line shows each step, and the line above it tells if a channel looks
# flipped. Select the next site in the Site menu, or open another survey with
# Survey…. Export writes the site as an EDI file and a PNG of the impedance
# (off-diagonal, or the full tensor with the switch) and the tipper, into a
# directory that you select. Their name holds the site, its base and remote
# sites and the values of the options (tf_filename). The window opens on the first site of the survey in
# alphabetical order. Set SURVEY to the full path of your survey directory, or
# give it as the argument. If that survey is not there, the window opens empty
#
#   julia --project=. -t auto examples/tkproc.jl D:/path/to/survey
#
# When you close the window, the script prints the estimates

using Timekeepers

# the full path of the survey directory:
const SURVEY = ""

survey = !isempty(ARGS) ? ARGS[1] : isdir(SURVEY) ? SURVEY : nothing
proc = survey === nothing ? run_tkproc() : run_tkproc(survey)

for (name, tf) in sort!(collect(proc.results); by = first)
    println(rpad(name, 16), tf)
end
