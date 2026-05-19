using Timekeepers

tv = load_metronix(joinpath(default_data_dir(), "DF090"); frequency = 128, site = "DF090")

println(tv)
