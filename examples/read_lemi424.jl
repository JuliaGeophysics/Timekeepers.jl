using Timekeepers

tv = load_lemi424(joinpath(default_data_dir(), "LEMI090.txt"); site = "LEMI090")

println(tv)
