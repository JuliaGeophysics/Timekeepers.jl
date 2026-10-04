# Timekeepers.jl - package entry point.
# Author: @pankajkmishra
#
# Timekeepers reads, inspects, edits and rewrites magnetotelluric and
# geomagnetic time series. It supports three instrument formats (LEMI-424,
# GEOMAG and Metronix ATS), a mask/amputation model for marking bad intervals,
# Welch spectral estimation, an interactive GLMakie explorer app, and a survey
# dashboard for choosing the base and remote sites of each local site
#
# This file defines the module, pulls in every component in dependency order
# and declares the public API

module Timekeepers

using Dates
using EzXML
using GLMakie
using NativeFileDialog
using PrecompileTools
using Printf
using Statistics
using TimeSeries

include("Types.jl")
include("Utilities.jl")
include("TimeArrayIO.jl")
include("Masking.jl")
include("Spectra.jl")
include("LEMI424.jl")
include("GEOMAG.jl")
include("MetronixATS.jl")
include("TimekeeperIO.jl")
include("Explorer.jl")
include("Survey.jl")
include("Dashboard.jl")
include("Precompile.jl")

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
export read_geomag
export load_geomag
export write_geomag
export METRONIX_CHANNEL_MAP
export read_metronix
export load_metronix
export write_metronix
export write_metronix_site
export write_metronix_site_masked
export split_metronix_site
export metronix_site_is_split
export metronix_site_rates
export load_metronix_site
export metronix_site_runs
export is_metronix_site
export TKApp
export run_tkapp

export SurveyRun
export SurveySite
export Survey
export scan_survey
export site_rates
export site_components
export has_magnetic
export survey_rates
export recording_intervals
export recording_seconds
export overlap_intervals
export overlap_seconds
export overlap_matrix
export site_distance
export site_references
export common_window
export reference_plan
export write_reference_plan
export read_reference_plan
export TKDash
export run_tkdash

end
