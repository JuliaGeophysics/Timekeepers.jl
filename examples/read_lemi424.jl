# read_lemi424.jl - minimal reading example.
# Author: @pankajkmishra
#
# Loads a bundled LEMI-424 file as a TimeArray and prints it, to show the
# shortest path from a file on disk to data you can work with.

using Timekeepers

tv = load_lemi424(joinpath(default_data_dir(), "LEMI090.txt"); site = "LEMI090")

println(tv)
