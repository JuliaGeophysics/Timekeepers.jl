# tkdash.jl - starts the survey dashboard.
# Author: @pankajkmishra
#
# This script scans a survey directory and opens TKDash. TKDash shows when each
# site recorded. For the site that you select, it shows the base sites (they
# recorded with it, near it) and the remote sites (they recorded with it, far
# from it). Give the survey directory as the first argument. If you give no
# argument, the script uses the BRGM Alsace survey next to this repository. If
# that survey is not there, a folder dialog asks for a directory
#
#   julia --project=. examples/tkdash.jl /path/to/survey
#
# Export in the window writes the base and remote sites of each site as a
# table. When you close the window, the script prints them

using Timekeepers

const DEFAULT_SURVEY = normpath(joinpath(@__DIR__, "..", "..", "MetronixATS"))

root = !isempty(ARGS) ? ARGS[1] : isdir(DEFAULT_SURVEY) ? DEFAULT_SURVEY : nothing
dash = root === nothing ? run_tkdash() : run_tkdash(root)

for r in reference_plan(dash)
    println(rpad(r.site, 16), "base: ", rpad(isempty(r.base) ? "–" : join(r.base, ", "), 36),
            "remote: ", isempty(r.remote) ? "–" : join(r.remote, ", "))
end
