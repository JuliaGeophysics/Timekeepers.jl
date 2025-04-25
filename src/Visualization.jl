function site_info_text(site_info)
    fields = [
        ("Site folder", get(site_info, "folder_name", "")),
        ("ADU", get(site_info, "adu_number", "")),
        ("Survey ID", get(site_info, "SurveyID", "")),
        ("Client", get(site_info, "Client", "")),
        ("Contractor", get(site_info, "Contractor", "")),
        ("Operator", get(site_info, "Operator", "")),
        ("Area", get(site_info, "Area", "")),
        ("Sampling Freq (Hz)", get(site_info, "sampling_frequency", "")),
        ("Duration", get(site_info, "duration_formatted", "")),
        ("Start Time", haskey(site_info, "start_datetime") && site_info["start_datetime"] !== nothing ? Dates.format(site_info["start_datetime"], "yyyy-mm-dd HH:MM:SS") : ""),
        ("End Time", haskey(site_info, "end_datetime") && site_info["end_datetime"] !== nothing ? Dates.format(site_info["end_datetime"], "yyyy-mm-dd HH:MM:SS") : ""),
        ("Latitude", @sprintf("%.6f", get(site_info, "ADU_Lat", NaN))),
        ("Longitude", @sprintf("%.6f", get(site_info, "ADU_Long", NaN))),
        ("Elevation (m)", @sprintf("%.2f", get(site_info, "ADU_Elev", NaN))),
        ("Weather", get(site_info, "Weather", "")),
        ("Comments", get(site_info, "Comments", ""))
    ]
    info_lines = [@sprintf("%-16s: %s", k, v) for (k, v) in fields if !(v == "" || occursin("NaN", string(v)))]
    return join(info_lines, "\n")
end

function print_site_info_table(site_info; site_name="DF090")
    # Prepare keys and values for display (choose your most relevant fields)
    fields = [
        ("Site name", site_name),
        ("Site folder", get(site_info, "folder_name", "")),
        ("ADU", get(site_info, "adu_number", "")),
        ("Survey ID", get(site_info, "SurveyID", "")),
        ("Client", get(site_info, "Client", "")),
        ("Contractor", get(site_info, "Contractor", "")),
        ("Operator", get(site_info, "Operator", "")),
        ("Area", get(site_info, "Area", "")),
        ("Sampling Freq (Hz)", get(site_info, "sampling_frequency", "")),
        ("Duration", get(site_info, "duration_formatted", "")),
        ("Start Time", haskey(site_info, "start_datetime") && site_info["start_datetime"] !== nothing ? Dates.format(site_info["start_datetime"], "yyyy-mm-dd HH:MM:SS") : ""),
        ("End Time", haskey(site_info, "end_datetime") && site_info["end_datetime"] !== nothing ? Dates.format(site_info["end_datetime"], "yyyy-mm-dd HH:MM:SS") : ""),
        ("Latitude", @sprintf("%.6f", get(site_info, "ADU_Lat", NaN))),
        ("Longitude", @sprintf("%.6f", get(site_info, "ADU_Long", NaN))),
        ("Elevation (m)", @sprintf("%.2f", get(site_info, "ADU_Elev", NaN))),
        ("Comments", get(site_info, "Comments", ""))
    ]

    # Filter out empty or NaN fields
    display_fields = [f for f in fields if !(f[2] == "" || occursin("NaN", string(f[2])))]

    # Compute column widths
    key_width = maximum(length(f[1]) for f in display_fields)
    val_width = maximum(length(f[2]) for f in display_fields)
    border_width = key_width + val_width + 7

    println("==================================================================================================") 
    for (key, val) in display_fields
        println("| " * rpad(key, key_width) * "  " * rpad(val, val_width) * " | ")
    end

    println("-" * repeat("-", border_width - 1))
end

function plot_time_series_terminal(time_axis, y_data, label, max_hours)
    # Create a plot with the data
    plt = lineplot(time_axis, y_data, 
                  width=100, height=6, 
                  xlabel="", 
                  ylabel=label,
                  # Force a specific range
                  xlim=(0, max_hours),
                  grid=false)
    
    display(plt)
end

function plot_all_components_terminal(data, site_info)
    comp_units = Dict("Ex"=>"mV/km", "Ey"=>"mV/km", "Hx"=>"nT", "Hy"=>"nT", "Hz"=>"nT")
    components = ["Ex", "Ey", "Hx", "Hy", "Hz"]
    valid_components = filter(c -> haskey(data, c), components)
    
    # Get the exact duration in hours
    max_hours = get(site_info, "duration_seconds", Inf) / 3600
    
    for comp in valid_components
        y_data = data[comp]["data"]
        fs = data[comp]["fs"]
        
        # Calculate the number of samples that fit in our exact duration
        max_samples = min(length(y_data), Int(round(max_hours * 3600 * fs)))
        
        # Create a fixed number of evenly spaced points
        num_display_points = min(12000, max_samples)
        
        # Use range to create EXACTLY evenly spaced points from 0 to max_hours
        t = range(0, max_hours, length=num_display_points)
        
        # Sample the data at these exact points
        indices = [min(max_samples, Int(round(ti * 3600 * fs)) + 1) for ti in t]
        y = y_data[indices]
        
        unit = get(comp_units, comp, "")
        # Pass max_hours to the plotting function
        plot_time_series_terminal(collect(t), y, "$comp ($unit)", max_hours)
    end
end


function create_pdf_plot(data, site_info, output_file)
    gr()
    comp_colors = Dict(
        "Ex" => :darkblue,
        "Ey" => :darkblue,
        "Hx" => :royalblue2,
        "Hy" => :royalblue2,
        "Hz" => :royalblue2
    )
    components = ["Ex", "Ey", "Hx", "Hy", "Hz"]
    valid_components = filter(c -> haskey(data, c), components)
    n_plots = length(valid_components)

    first_comp = valid_components[1]
    fs = data[first_comp]["fs"]
    total_samples = length(data[first_comp]["data"])
    duration_seconds = total_samples / fs

    max_display_points = 15000
    decimation_factor = max(1, ceil(Int, total_samples / max_display_points))
    time_axis = (0:decimation_factor:total_samples-1) ./ fs / 3600

    # Prepare start and end datetime strings
    if haskey(site_info, "start_datetime") && site_info["start_datetime"] !== nothing
        start_dt = site_info["start_datetime"]
        end_dt = start_dt + Dates.Second(round(Int, duration_seconds))
        xticks = ([time_axis[1], time_axis[end]],
                  [Dates.format(start_dt, "yyyy-mm-dd HH:MM:SS"),
                   Dates.format(end_dt, "yyyy-mm-dd HH:MM:SS")])
        xlabel = "Time (UTC)"
        startstr = Dates.format(start_dt, "yyyy-mm-dd HH:MM:SS")
        endstr = Dates.format(end_dt, "yyyy-mm-dd HH:MM:SS")
    else
        xticks = ([time_axis[1], time_axis[end]], ["0 h", @sprintf("%.2f h", duration_seconds/3600)])
        xlabel = "Elapsed Time (hours)"
        startstr = "Unknown"
        endstr = "Unknown"
    end

    comp_units = Dict("Ex"=>"mV/km", "Ey"=>"mV/km", "Hx"=>"nT", "Hy"=>"nT", "Hz"=>"nT")

    plots = []
    for (i, comp) in enumerate(valid_components)
        y_data = data[comp]["data"][1:decimation_factor:end]
        color = get(comp_colors, comp, :black)
        unit = get(comp_units, comp, "")
        ylims = (minimum(y_data), maximum(y_data))
        p = plot(
            time_axis, y_data,
            color=color, linewidth=0.8, legend=false,
            ylabel=comp*" ("*unit*")",
            xlabel=(i == n_plots ? xlabel : ""),
            xticks=(i == n_plots ? xticks : false),
            yticks=:auto,
            framestyle=:box,
            grid=false,
            minorgrid=false,
            left_margin=20Plots.mm,
            bottom_margin=(i == n_plots ? 15Plots.mm : 2Plots.mm),
            top_margin=(i == 1 ? 8Plots.mm : 2Plots.mm),
            xlims=(time_axis[1], time_axis[end]),
            ylims=ylims,
            tickfontsize=10, labelfontsize=13,
            dpi=300,
            background_color=:white
        )
        push!(plots, p)
    end

    # Only save the time series panels as PDF, no summary page
    final_plot = plot(plots...,
        layout=(n_plots, 1), link=:x,
        size=(1200, 240+150*n_plots),
        margin=8Plots.mm,
    )
    savefig(final_plot, output_file)
    println("PDF plot saved as $output_file")
end