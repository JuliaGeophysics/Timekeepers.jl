# tkdash.jl - starts the survey dashboard.
# Author: @pankajkmishra
#
# This script opens TKDash. A folder dialog asks for the survey directory,
# then TKDash scans it and shows when each site recorded. For the site that
# you select, it shows the base sites (they recorded with it, near it) and
# the remote sites (they recorded with it, far from it)
#
#   julia --project=. examples/tkdash.jl
#
# Export in the window writes the base and remote sites of each site as a
# table

using Timekeepers

run_tkdash()
