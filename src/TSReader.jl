"""
    read_ts(frequency)
    
Reads all 5 components of magnetotelluric data from DF090 directory
at the specified sampling frequency.

Parameters:
- frequency: The sampling frequency to read (e.g., 128 for 128Hz data)

Returns:
- The data dictionary containing all component data and site info
"""
function read_ts(frequency)
    # Path to the DF090 directory
    base_dir = "DF090"
    
    println("Looking for measurement directories in: $(abspath(base_dir))")
    
    # Check if base directory exists
    if !isdir(base_dir)
        error("Base directory DF090 not found")
    end
    
    # Find measurement folders (starting with "meas_")
    meas_dirs = String[]
    for entry in readdir(base_dir)
        full_path = joinpath(base_dir, entry)
        if isdir(full_path) && startswith(entry, "meas_")
            push!(meas_dirs, full_path)
            println("Found measurement directory: $full_path")
        end
    end
    
    if isempty(meas_dirs)
        error("No measurement directories found in DF090 folder")
    end
    
    # Frequency tag to look for
    freq_tag = "$(frequency)H"
    println("Looking for files with frequency tag: $freq_tag")
    
    # Component mappings
    component_map = Dict(
        "TEx" => "Ex",
        "TEy" => "Ey", 
        "THx" => "Hx",
        "THy" => "Hy",
        "THz" => "Hz"
    )
    
    # Process each measurement directory
    data = Dict{String, Dict{String, Any}}()
    target_dir = nothing
    site_info = Dict{String, Any}()
    
    for dir in meas_dirs
        println("Examining directory: $dir")
        dir_files = readdir(dir)
        
        # Find all .ats files with our target frequency
        ats_files = filter(f -> endswith(lowercase(f), "$(lowercase(freq_tag)).ats"), dir_files)
        
        println("Found $(length(ats_files)) ATS files with frequency $freq_tag")
        
        if !isempty(ats_files)
            target_dir = dir
            println("Found target directory with $(frequency)Hz data: $dir")
            
            # Store folder name for display
            site_info["folder_name"] = basename(dir)
            
            # Try to extract start time from directory name
            start_time_str = ""
            if startswith(basename(dir), "meas_")
                start_time_str = replace(basename(dir), "meas_" => "")
                try
                    site_info["start_datetime"] = Dates.DateTime(start_time_str, "yyyy-mm-dd_HH-MM-SS")
                catch
                    site_info["start_datetime"] = nothing
                end
            end
            
            # Process each component
            first_file = ""
            for (metronix_comp, mt_comp) in component_map
                # Find corresponding file
                comp_files = filter(f -> occursin(metronix_comp, f) && 
                                      occursin(freq_tag, f) && 
                                      endswith(f, ".ats"), 
                                    dir_files)
                
                if !isempty(comp_files)
                    data_file = joinpath(dir, comp_files[1])
                    if first_file == ""
                        first_file = data_file
                    end
                    println("Found $mt_comp file: $data_file")
                    
                    # Read time series data
                    println("Reading time series data from: $data_file")
                    ts_data, fs, start_time_unix, header = read_ats_sample(data_file)
                    println("Successfully read $(length(ts_data)) samples for $mt_comp")
                    
                    # Store in our data dictionary
                    data[mt_comp] = Dict(
                        "header" => header,
                        "data" => ts_data,
                        "fs" => fs,
                        "start_time" => start_time_unix
                    )
                else
                    println("WARNING - No file found for $mt_comp component")
                end
            end
            
            # Extract ADU number from filename
            if first_file != ""
                site_info["adu_number"] = extract_adu_number(first_file)
            end
            
            # Break after finding first directory with target frequency data
            break
        end
    end
    
    if isempty(data)
        error("No $(frequency)Hz data could be found in any meas_ directory")
    end
    
    # Calculate additional site information
    first_comp = first(keys(data))
    site_info["sampling_frequency"] = data[first_comp]["fs"]
    
    # Calculate duration and end time
    total_samples = length(data[first_comp]["data"])
    duration_seconds = total_samples / site_info["sampling_frequency"]
    site_info["duration_seconds"] = duration_seconds
    site_info["duration_formatted"] = format_duration(duration_seconds)
    
    if haskey(site_info, "start_datetime") && site_info["start_datetime"] !== nothing
        site_info["end_datetime"] = site_info["start_datetime"] + Dates.Second(round(Int, duration_seconds))
    else
        site_info["end_datetime"] = nothing
    end
    
    # Copy any other useful info from the header
    if haskey(data, first_comp) && haskey(data[first_comp], "header")
        header = data[first_comp]["header"]
        for field in ["ADU_Lat", "ADU_Long", "ADU_Elev", "Client", "Contractor", 
                     "Area", "SurveyID", "Operator", "Weather", "Comments"]
            if haskey(header, field) && header[field] != ""
                site_info[field] = header[field]
            end
        end
    end
    
    return data, site_info
end