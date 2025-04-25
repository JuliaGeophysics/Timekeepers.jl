module EMTimeSeries


using Plots
using Statistics
using Dates
using Printf
using UnicodePlots


include("ATSReader.jl")
include("Helpers.jl")
include("TSReader.jl") 
include("Visualization.jl")
include("ATSWriter.jl")


export read_ats_header
export read_ats_sample
export write_ats_file
export extract_adu_number
export format_duration
export read_ts
export site_info_text
export print_site_info_table
export plot_time_series_terminal
export plot_all_components_terminal
export create_pdf_plot

end 