# tkproc.jl - starts the transfer function window.
# Author: @pankajkmishra
#
# This script opens TKProc without a survey. Use its Survey… button to select
# a survey directory. TKProc then lists its sites with electric channels in
# the Site menu. When you select a site, TKProc loads its base sites (they
# recorded with it, near it) and its remote sites (they recorded with it, far
# from it) into their menus, with the best set-up selected. If you selected
# the base and remote sites of each site in TKDash and exported them
# (reference_plan.txt in the survey directory), TKProc uses that plan: the
# menus hold the sites of the plan, one base and one remote site at a time.
# Without a plan, TKProc makes the lists itself
#
# Start Julia with more than one thread so the window stays live while a site
# is processed:
#
#   julia --project=. -t auto examples/tkproc.jl

using Timekeepers

run_tkproc()
