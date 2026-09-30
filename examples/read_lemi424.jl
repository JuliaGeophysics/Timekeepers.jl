# read_lemi424.jl - minimal reading example.
# Author: @pankajkmishra
#
# Loads data/LEMI090.txt as a TimeArray and prints it, to show the shortest
# path from a file on disk to data you can work with. No recordings ship with
# the package; the docs point to a public LEMI-424 dataset to put there.

using Timekeepers

tv = load_lemi424(joinpath(default_data_dir(), "LEMI090.txt"); site = "LEMI090")

println(tv)
