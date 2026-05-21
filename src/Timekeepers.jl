module Timekeepers

using Dates
using GLMakie
using NativeFileDialog
using Printf
using Statistics
using TimeSeries

include("Types.jl")
include("Utilities.jl")
include("TimeArrayIO.jl")
include("Masking.jl")
include("Spectra.jl")
include("LEMI424.jl")
include("TimekeeperIO.jl")
include("Explorer.jl")

export TimekeeperChannel
export TimekeeperRun
export default_data_dir
export components
export default_components
export sampling_rate
export start_time
export end_time
export duration_seconds
export to_timearray
export from_timearray
export TimekeeperMask
export mask_interval!
export unmask_interval!
export clear_mask!
export masked_samples
export sample_weights
export cleaned_timearray
export good_segments
export combine_masks
export write_cleaned
export write_mask
export read_mask
export read_timekeeper
export write_timekeeper

export LEMI424_COLUMNS
export read_lemi424
export load_lemi424
export write_lemi424
export TKApp
export run_tkapp

end
