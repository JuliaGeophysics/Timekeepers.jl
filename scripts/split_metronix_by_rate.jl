#!/usr/bin/env julia
#
# Split a raw Metronix ADU site (a folder of meas_* directories that may mix
# sampling rates) into one single-rate site directory per sampling rate:
#
#     RK137/  ->  RK137.TK128/   (all 128 Hz meas_ dirs)
#                 RK137.TK4096/  (all 4096 Hz meas_ dirs)
#                 RK137.TK131072/
#
# meas_ directories are copied verbatim (lossless), grouped by the sampling
# rate detected from each run's ATS header. The Timekeepers GUI / loaders
# expect one of these single-rate .TK<rate> directories.
#
# Usage:
#     julia --project=. scripts/split_metronix_by_rate.jl <site_dir> [<site_dir> ...]
#     julia --project=. scripts/split_metronix_by_rate.jl --move <site_dir>   # move instead of copy

using Timekeepers

"""
    split_metronix_site_by_rate(site_dir; move=false, dest_parent=dirname(site_dir)) -> Vector{String}

Group the `meas_*` directories of `site_dir` by sampling rate and copy (or
move) each into a sibling `<sitename>.TK<rate>` directory. Returns the created
directories.
"""
function split_metronix_site_by_rate(site_dir::AbstractString;
                                     move::Bool = false,
                                     dest_parent::Union{Nothing, AbstractString} = nothing)
    site_dir = rstrip(abspath(site_dir), ['/', '\\'])
    isdir(site_dir) || error("Not a directory: $site_dir")
    site = basename(site_dir)
    runs = metronix_site_runs(site_dir)
    isempty(runs) && (@warn "No Metronix meas_ directories found" site_dir; return String[])
    parent = dest_parent === nothing ? dirname(site_dir) : String(dest_parent)

    dests = String[]
    for (rate, meas_dirs) in sort(collect(runs); by = first)
        dest = joinpath(parent, "$(site).TK$(string(round(Int, rate)))")
        mkpath(dest)
        for meas in sort(meas_dirs)
            target = joinpath(dest, basename(meas))
            isdir(target) && rm(target; recursive = true)
            move ? mv(meas, target) : cp(meas, target)
        end
        println("  $(round(Int, rate)) Hz -> $(basename(dest))  ($(length(meas_dirs)) run(s))")
        push!(dests, dest)
    end
    return dests
end

function _main(args)
    move = false
    sites = String[]
    for a in args
        if a == "--move"
            move = true
        elseif a in ("-h", "--help")
            println("usage: julia --project=. scripts/split_metronix_by_rate.jl [--move] <site_dir> ...")
            return
        else
            push!(sites, a)
        end
    end
    if isempty(sites)
        println("usage: julia --project=. scripts/split_metronix_by_rate.jl [--move] <site_dir> ...")
        return
    end
    for s in sites
        println("Splitting $s by sampling rate:")
        split_metronix_site_by_rate(s; move = move)
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    _main(ARGS)
end
