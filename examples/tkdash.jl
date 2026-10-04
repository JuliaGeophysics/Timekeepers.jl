# tkdash.jl - launch the survey dashboard.
# Author: @pankajkmishra
#
# Scans a survey directory and opens TKDash: when every site recorded and, for
# the site picked, its base sites (recorded with it, close by) and
# remote sites (recorded with it, far away). Pass the survey directory as the
# first argument; with none, the BRGM Alsace survey beside this repository is
# used, or a folder dialog asks
#
#   julia --project=. examples/tkdash.jl /path/to/survey
#
# Export in the window writes every site's base and remote sites as a table;
# closing the window prints them

using Timekeepers

const DEFAULT_SURVEY = normpath(joinpath(@__DIR__, "..", "..", "MetronixATS"))

root = !isempty(ARGS) ? ARGS[1] : isdir(DEFAULT_SURVEY) ? DEFAULT_SURVEY : nothing
dash = root === nothing ? run_tkdash() : run_tkdash(root)

for r in reference_plan(dash)
    println(rpad(r.site, 16), "base: ", rpad(isempty(r.base) ? "–" : join(r.base, ", "), 36),
            "remote: ", isempty(r.remote) ? "–" : join(r.remote, ", "))
end
