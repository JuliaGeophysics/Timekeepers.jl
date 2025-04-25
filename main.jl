#!/usr/bin/env julia

using Dates
using Plots
using Statistics
using Printf
using UnicodePlots

include("src/EMTimeSeries.jl")
using .EMTimeSeries

# Read the time series data
frequency = 128        
data, site_info, ats_files_map, xml_files = read_ts(frequency)

# Display site information
print_site_info_table(site_info)
plot_all_components_terminal(data, site_info)

# Create an output directory name (you can customize this)
output_dir = "ProcessedDF090_$(Dates.format(now(), "yyyymmdd_HHMMSS"))"

# Call the function to write all files to the output directory
processed_files = write_all_ats_files(output_dir, data, site_info, ats_files_map, xml_files)

println("All files processed successfully and saved to: $output_dir")