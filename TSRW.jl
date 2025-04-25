using Plots
using Statistics
using Dates
using Printf
using DelimitedFiles
using FileIO
using FFTW
using StatsBase  # For faster statistical operations

"""
    robust_despike_timeseries(data::Vector, threshold_factor::Float64=3.5, window_size::Int=50)
    
A robust implementation of the despiking algorithm that avoids threading issues
and is optimized for speed.
"""
function robust_despike_timeseries(data::Vector, threshold_factor::Float64=3.5, window_size::Int=50)
    n = length(data)
    despiked_data = copy(data)
    spike_indices = Int[]
    
    println("Starting despiking with window size $window_size and threshold $threshold_factor...")
    
    # First pass: detect and remove major spikes
    # Process in chunks to improve memory locality and performance
    chunk_size = 10000
    num_chunks = ceil(Int, n / chunk_size)
    
    for chunk_idx in 1:num_chunks
        start_idx = (chunk_idx - 1) * chunk_size + 1
        end_idx = min(chunk_idx * chunk_size, n)
        
        if chunk_idx % 10 == 0 || chunk_idx == 1 || chunk_idx == num_chunks
            println("Processing chunk $chunk_idx/$num_chunks (indices $start_idx:$end_idx)...")
        end
        
        # For each point in the chunk, check if it's a spike
        for i in start_idx:end_idx
            # Define window boundaries
            half_window = div(window_size, 2)
            w_start = max(1, i - half_window)
            w_end = min(n, i + half_window)
            
            # Get window data excluding the current point
            window_indices = vcat(w_start:(i-1), (i+1):w_end)
            window_data = @view data[window_indices]
            
            # Calculate median and MAD
            window_median = median(window_data)
            window_mad = median(abs.(window_data .- window_median)) * 1.4826  # Convert to equivalent std dev
            
            # Check if the point is an outlier
            if window_mad > 1e-10 && abs(data[i] - window_median) > threshold_factor * window_mad
                despiked_data[i] = window_median
                push!(spike_indices, i)
            end
        end
    end
    
    # Second pass for more subtle spikes (sample-based approach to save time)
    if !isempty(spike_indices)
        second_threshold = threshold_factor * 0.85
        second_window = div(window_size, 2)
        max_samples = min(5000, n ÷ 10)  # Sample at most 5000 points or 10% of data
        
        # Determine points to check (exclude already identified spikes)
        remaining_points = setdiff(1:n, spike_indices)
        
        # Sample a subset of points to check for second-pass spikes
        if length(remaining_points) > max_samples
            sample_points = sample(remaining_points, max_samples, replace=false)
        else
            sample_points = remaining_points
        end
        
        println("Second pass: checking $(length(sample_points)) sample points for remaining spikes...")
        
        new_spikes = Int[]
        for i in sample_points
            # Define window boundaries
            half_window = div(second_window, 2)
            w_start = max(1, i - half_window)
            w_end = min(n, i + half_window)
            
            # Get window indices excluding known spikes and current point
            window_indices = filter(j -> j != i && !(j in spike_indices), w_start:w_end)
            
            if length(window_indices) >= 3  # Need enough data for statistics
                window_data = despiked_data[window_indices]
                window_median = median(window_data)
                window_mad = median(abs.(window_data .- window_median)) * 1.4826
                
                # Check if point is an outlier in second pass
                if window_mad > 1e-10 && abs(despiked_data[i] - window_median) > second_threshold * window_mad
                    despiked_data[i] = window_median
                    push!(new_spikes, i)
                end
            end
        end
        
        append!(spike_indices, new_spikes)
        println("Second pass identified additional $(length(new_spikes)) spikes")
    end
    
    println("Total: Identified and removed $(length(spike_indices)) spikes ($(round(100.0 * length(spike_indices) / n, digits=2))%)")
    
    return despiked_data, spike_indices
end

"""
    read_ats_header(filename::String)
    
Read the header of a Metronix ATS file according to the format specification.
"""
function read_ats_header(filename::String)
    println("Reading ATS header from: $filename")
    
    header = Dict{String, Any}()
    
    open(filename, "r") do f
        # Read header length first (to verify)
        header_length = read(f, UInt16)
        seekstart(f)
        
        # Read the fields sequentially
        header["header_length"] = read(f, UInt16)
        header["header_version"] = read(f, Int16)
        header["sample_length"] = read(f, Int32)
        header["sampling_rate"] = read(f, Float32)
        header["start"] = read(f, Int32)
        header["lsbval"] = read(f, Float64)
        header["GMToffset"] = read(f, Int32)
        header["Res1"] = read(f, Int32)
        header["serial_number_ADU06"] = read(f, UInt16)
        header["serial_number_ADC_board"] = read(f, Int16)
        header["channel_number"] = read(f, Int8)
        header["Res2"] = read(f, Int8)
        
        # Read string fields as byte arrays and convert
        channel_type_bytes = Vector{UInt8}(undef, 2)
        read!(f, channel_type_bytes)
        header["channel_type"] = String(filter(x -> x != 0x00, channel_type_bytes))
        
        sensor_type_bytes = Vector{UInt8}(undef, 6)
        read!(f, sensor_type_bytes)
        header["sensor_type"] = String(filter(x -> x != 0x00, sensor_type_bytes))
        
        header["sensor_serial_number"] = read(f, Int16)
        header["x1"] = read(f, Float32)
        header["y1"] = read(f, Float32)
        header["z1"] = read(f, Float32)
        header["x2"] = read(f, Float32)
        header["y2"] = read(f, Float32)
        header["z2"] = read(f, Float32)
        header["E_field_dipole_length"] = read(f, Float32)
        header["angle"] = read(f, Float32)
        header["rho_probe"] = read(f, Float32)
        header["DC_offset_voltage"] = read(f, Float32)
        header["internal_gain_ampli"] = read(f, Float32)
        header["Res3"] = read(f, Int32)
        header["ADU_Lat"] = read(f, Int32)
        header["ADU_Long"] = read(f, Int32)
        header["ADU_Elev"] = read(f, Int32)
        
        lat_long_type_bytes = Vector{UInt8}(undef, 1)
        read!(f, lat_long_type_bytes)
        header["Lat_Long_TYPE"] = String(filter(x -> x != 0x00, lat_long_type_bytes))
        
        add_coordinates_bytes = Vector{UInt8}(undef, 1)
        read!(f, add_coordinates_bytes)
        header["add_coordinates"] = String(filter(x -> x != 0x00, add_coordinates_bytes))
        
        header["ref_meridian"] = read(f, Int16)
        header["xcoor"] = read(f, Float64)
        header["ycoor"] = read(f, Float64)
        
        gps_clock_status_bytes = Vector{UInt8}(undef, 1)
        read!(f, gps_clock_status_bytes)
        header["gps_clock_status"] = String(filter(x -> x != 0x00, gps_clock_status_bytes))
        
        header["accuracy_GPS"] = read(f, Int8)
        header["offset_UTC"] = read(f, Int16)
        header["Res4a"] = read(f, Int32)
        header["Res4b"] = read(f, Int32)
        header["Res4c"] = read(f, Int32)
        
        survey_header_filename_bytes = Vector{UInt8}(undef, 12)
        read!(f, survey_header_filename_bytes)
        header["survey_header_filename"] = String(filter(x -> x != 0x00, survey_header_filename_bytes))
        
        type_of_meas_bytes = Vector{UInt8}(undef, 4)
        read!(f, type_of_meas_bytes)
        header["type_of_meas"] = String(filter(x -> x != 0x00, type_of_meas_bytes))
        
        logfile_bytes = Vector{UInt8}(undef, 12)
        read!(f, logfile_bytes)
        header["logfile"] = String(filter(x -> x != 0x00, logfile_bytes))
        
        result_selftest_bytes = Vector{UInt8}(undef, 2)
        read!(f, result_selftest_bytes)
        header["result_selftest"] = String(filter(x -> x != 0x00, result_selftest_bytes))
        
        header["Res5"] = read(f, Int16)
        header["number_of_calib_freq"] = read(f, Int16)
        header["length_of_freq_entry"] = read(f, Int16)
        header["version_calib"] = read(f, Int16)
        header["start_addres"] = read(f, Int16)
        header["Res6"] = read(f, Int64)
        
        cal_filename_ADU06_bytes = Vector{UInt8}(undef, 12)
        read!(f, cal_filename_ADU06_bytes)
        header["cal_filename_ADU06"] = String(filter(x -> x != 0x00, cal_filename_ADU06_bytes))
        
        header["datetime_calib"] = read(f, Int32)
        
        cal_sensor_filename_bytes = Vector{UInt8}(undef, 12)
        read!(f, cal_sensor_filename_bytes)
        header["cal_sensor_filename"] = String(filter(x -> x != 0x00, cal_sensor_filename_bytes))
        
        header["datetime_calib_sens"] = read(f, Int32)
        header["powerline1"] = read(f, Float32)
        header["powerline2"] = read(f, Float32)
        header["Res7"] = read(f, Int64)
        header["CSAMT_Tx_freq"] = read(f, Float32)
        header["CSAMT_TS_blocks"] = read(f, Int16)
        header["CSAMT_stacks"] = read(f, Int16)
        header["CSAMT_blk_length"] = read(f, Int32)
        header["Res8"] = read(f, Int32)
        
        client_bytes = Vector{UInt8}(undef, 16)
        read!(f, client_bytes)
        header["Client"] = String(filter(x -> x != 0x00, client_bytes))
        
        contractor_bytes = Vector{UInt8}(undef, 16)
        read!(f, contractor_bytes)
        header["Contractor"] = String(filter(x -> x != 0x00, contractor_bytes))
        
        area_bytes = Vector{UInt8}(undef, 16)
        read!(f, area_bytes)
        header["Area"] = String(filter(x -> x != 0x00, area_bytes))
        
        surveyid_bytes = Vector{UInt8}(undef, 16)
        read!(f, surveyid_bytes)
        header["SurveyID"] = String(filter(x -> x != 0x00, surveyid_bytes))
        
        operator_bytes = Vector{UInt8}(undef, 16)
        read!(f, operator_bytes)
        header["Operator"] = String(filter(x -> x != 0x00, operator_bytes))
        
        # Skip Res9
        res9_bytes = Vector{UInt8}(undef, 112)
        read!(f, res9_bytes)
        
        weather_bytes = Vector{UInt8}(undef, 64)
        read!(f, weather_bytes)
        header["Weather"] = String(filter(x -> x != 0x00, weather_bytes))
        
        comments_bytes = Vector{UInt8}(undef, 512)
        read!(f, comments_bytes)
        header["Comments"] = String(filter(x -> x != 0x00, comments_bytes))
    end

    # Apply corrections as in Python code
    header["ADU_Lat"] *= 1e-3 / 3600
    header["ADU_Long"] *= 1e-3 / 3600
    header["ADU_Elev"] *= 1e-2

    return header
end

"""
    read_ats_sample(filename::String)
    
Read the time series data from a Metronix ATS file.
Returns (data, sampling_rate, start_time, header)
"""
function read_ats_sample(filename::String)
    println("Reading ATS sample data from: $filename")
    
    header = read_ats_header(filename)
    sampling_rate = header["sampling_rate"]
    start_time = header["start"]
    header_length = header["header_length"]
    sample_length = header["sample_length"]
    lsbval = header["lsbval"]

    println("Reading $(sample_length) samples...")
    
    # Read the raw int32 data
    data = Vector{Int32}(undef, sample_length)
    open(filename, "r") do f
        seek(f, header_length)
        read!(f, data)
    end
    
    # Apply LSB scaling factor to convert to physical units
    data_scaled = Float64.(data) .* lsbval
    
    # Print data statistics
    println("Data statistics (unscaled) - min: $(minimum(data)), max: $(maximum(data))")
    println("Data statistics (scaled) - min: $(minimum(data_scaled)), max: $(maximum(data_scaled)), mean: $(mean(data_scaled))")
    
    return data_scaled, sampling_rate, start_time, header
end

"""
    write_ats_file(filename::String, data, header)
    
Write time series data to a Metronix ATS file with the given header.
"""
function write_ats_file(filename::String, data, header)
    println("Writing ATS file: $filename")
    
    # Convert float data back to Int32 using the LSB value
    lsbval = header["lsbval"]
    raw_data = round.(Int32, data ./ lsbval)
    
    open(filename, "w") do f
        # Write header
        write(f, UInt16(header["header_length"]))
        write(f, Int16(header["header_version"]))
        write(f, Int32(header["sample_length"]))
        write(f, Float32(header["sampling_rate"]))
        write(f, Int32(header["start"]))
        write(f, Float64(header["lsbval"]))
        write(f, Int32(header["GMToffset"]))
        write(f, Int32(header["Res1"]))
        write(f, UInt16(header["serial_number_ADU06"]))
        write(f, Int16(header["serial_number_ADC_board"]))
        write(f, Int8(header["channel_number"]))
        write(f, Int8(header["Res2"]))
        
        # Write string fields with proper zero-padding
        channel_type = Vector{UInt8}(undef, 2)
        channel_type .= 0x00
        ct_bytes = Vector{UInt8}(header["channel_type"])
        channel_type[1:min(length(ct_bytes), 2)] = ct_bytes[1:min(length(ct_bytes), 2)]
        write(f, channel_type)
        
        sensor_type = Vector{UInt8}(undef, 6)
        sensor_type .= 0x00
        st_bytes = Vector{UInt8}(header["sensor_type"])
        sensor_type[1:min(length(st_bytes), 6)] = st_bytes[1:min(length(st_bytes), 6)]
        write(f, sensor_type)
        
        write(f, Int16(header["sensor_serial_number"]))
        write(f, Float32(header["x1"]))
        write(f, Float32(header["y1"]))
        write(f, Float32(header["z1"]))
        write(f, Float32(header["x2"]))
        write(f, Float32(header["y2"]))
        write(f, Float32(header["z2"]))
        write(f, Float32(header["E_field_dipole_length"]))
        write(f, Float32(header["angle"]))
        write(f, Float32(header["rho_probe"]))
        write(f, Float32(header["DC_offset_voltage"]))
        write(f, Float32(header["internal_gain_ampli"]))
        write(f, Int32(header["Res3"]))
        
        # Convert lat/long back to original format
        write(f, Int32(round(header["ADU_Lat"] * 3600 * 1e3)))
        write(f, Int32(round(header["ADU_Long"] * 3600 * 1e3)))
        write(f, Int32(round(header["ADU_Elev"] * 1e2)))
        
        lat_long_type = Vector{UInt8}(undef, 1)
        lat_long_type .= 0x00
        if haskey(header, "Lat_Long_TYPE") && !isempty(header["Lat_Long_TYPE"])
            lat_long_type[1] = Vector{UInt8}(header["Lat_Long_TYPE"])[1]
        end
        write(f, lat_long_type)
        
        add_coordinates = Vector{UInt8}(undef, 1)
        add_coordinates .= 0x00
        if haskey(header, "add_coordinates") && !isempty(header["add_coordinates"])
            add_coordinates[1] = Vector{UInt8}(header["add_coordinates"])[1]
        end
        write(f, add_coordinates)
        
        write(f, Int16(header["ref_meridian"]))
        write(f, Float64(header["xcoor"]))
        write(f, Float64(header["ycoor"]))
        
        gps_clock_status = Vector{UInt8}(undef, 1)
        gps_clock_status .= 0x00
        if haskey(header, "gps_clock_status") && !isempty(header["gps_clock_status"])
            gps_clock_status[1] = Vector{UInt8}(header["gps_clock_status"])[1]
        end
        write(f, gps_clock_status)
        
        write(f, Int8(header["accuracy_GPS"]))
        write(f, Int16(header["offset_UTC"]))
        write(f, Int32(header["Res4a"]))
        write(f, Int32(header["Res4b"]))
        write(f, Int32(header["Res4c"]))
        
        # Write string fields with proper zero-padding
        function write_padded_string(fieldname, buf_length)
            buffer = Vector{UInt8}(undef, buf_length)
            buffer .= 0x00
            if haskey(header, fieldname) && !isempty(header[fieldname])
                field_bytes = Vector{UInt8}(header[fieldname])
                buffer[1:min(length(field_bytes), buf_length)] = field_bytes[1:min(length(field_bytes), buf_length)]
            end
            write(f, buffer)
        end
        
        write_padded_string("survey_header_filename", 12)
        write_padded_string("type_of_meas", 4)
        write_padded_string("logfile", 12)
        write_padded_string("result_selftest", 2)
        
        write(f, Int16(header["Res5"]))
        write(f, Int16(header["number_of_calib_freq"]))
        write(f, Int16(header["length_of_freq_entry"]))
        write(f, Int16(header["version_calib"]))
        write(f, Int16(header["start_addres"]))
        write(f, Int64(header["Res6"]))
        
        write_padded_string("cal_filename_ADU06", 12)
        write(f, Int32(header["datetime_calib"]))
        write_padded_string("cal_sensor_filename", 12)
        write(f, Int32(header["datetime_calib_sens"]))
        write(f, Float32(header["powerline1"]))
        write(f, Float32(header["powerline2"]))
        write(f, Int64(header["Res7"]))
        write(f, Float32(header["CSAMT_Tx_freq"]))
        write(f, Int16(header["CSAMT_TS_blocks"]))
        write(f, Int16(header["CSAMT_stacks"]))
        write(f, Int32(header["CSAMT_blk_length"]))
        write(f, Int32(header["Res8"]))
        
        write_padded_string("Client", 16)
        write_padded_string("Contractor", 16)
        write_padded_string("Area", 16)
        write_padded_string("SurveyID", 16)
        write_padded_string("Operator", 16)
        
        # Write Res9 as zeros
        res9 = Vector{UInt8}(undef, 112)
        res9 .= 0x00
        write(f, res9)
        
        write_padded_string("Weather", 64)
        write_padded_string("Comments", 512)
        
        # Write the data
        write(f, raw_data)
    end
    
    println("Successfully wrote $filename with $(length(data)) samples")
end

"""
    extract_adu_number(filename::String)

Extract the ADU number from the filename of an ATS file.
"""
function extract_adu_number(filename::String)
    # Get the base filename
    base = basename(filename)
    # Extract the ADU number (typically the first digits in the filename)
    match_result = match(r"^(\d+)_", base)
    if match_result !== nothing
        return match_result[1]
    else
        return "Unknown"
    end
end

"""
    format_duration(seconds)
    
Format duration in seconds to HH:MM:SS format
"""
function format_duration(seconds)
    hours = Int(floor(seconds / 3600))
    minutes = Int(floor((seconds - hours * 3600) / 60))
    secs = Int(round(seconds - hours * 3600 - minutes * 60))
    return @sprintf("%02d:%02d:%02d", hours, minutes, secs)
end

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
    ats_files_map = Dict{String, String}()  # Store original file paths
    xml_files = String[]  # Store XML file paths
    target_dirs = String[]  # Store directories with data
    site_info = Dict{String, Any}()
    
    for dir in meas_dirs
        println("Examining directory: $dir")
        dir_files = readdir(dir)
        
        # Find all .ats files with our target frequency
        ats_files = filter(f -> endswith(lowercase(f), "$(lowercase(freq_tag)).ats"), dir_files)
        
        println("Found $(length(ats_files)) ATS files with frequency $freq_tag")
        
        if !isempty(ats_files)
            push!(target_dirs, dir)
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
                    
                    # Store the original file path
                    ats_files_map[mt_comp] = data_file
                else
                    println("WARNING - No file found for $mt_comp component")
                end
            end
            
            # Find XML files in the same directory
            xml_files = filter(f -> endswith(lowercase(f), ".xml"), dir_files)
            xml_files = [joinpath(dir, f) for f in xml_files]
            
            # Extract ADU number from filename
            if first_file != ""
                site_info["adu_number"] = extract_adu_number(first_file)
            end
            
            # Don't break after first directory - collect all data
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
    
    return data, site_info, ats_files_map, xml_files, target_dirs
end

"""
    despike_and_save_data(data, site_info, ats_files_map, xml_files, target_dirs, output_dir)
    
Applies despiking to all components in the data and saves to output_dir,
preserving the original directory structure.
"""
function despike_and_save_data(data, site_info, ats_files_map, xml_files, target_dirs, output_dir)
    # Ensure output directory exists
    if !isdir(output_dir)
        println("Creating output directory: $output_dir")
        mkpath(output_dir)
    end
    
    # Create despiked data dictionary
    despiked_data = Dict{String, Dict{String, Any}}()
    despike_stats = Dict{String, Any}()
    
    # Map original paths to target paths
    dir_mapping = Dict{String, String}()
    for dir in target_dirs
        rel_path = relpath(dir, "DF090")
        target_path = joinpath(output_dir, rel_path)
        dir_mapping[dir] = target_path
        
        # Create the target directory if it doesn't exist
        if !isdir(target_path)
            println("Creating directory: $target_path")
            mkpath(target_path)
        end
    end
    
    # Process each component
    for comp in ["Ex", "Ey", "Hx", "Hy", "Hz"]
        if haskey(data, comp)
            # Apply despiking with the robust implementation
            println("Despiking $comp component...")
            original_data = data[comp]["data"]
            despiked, spike_indices = robust_despike_timeseries(original_data)
            
            # Copy original data structure and update with despiked data
            despiked_data[comp] = deepcopy(data[comp])
            despiked_data[comp]["data"] = despiked
            
            # Save statistics about spikes
            spike_count = length(spike_indices)
            spike_percentage = 100.0 * spike_count / length(original_data)
            println("Removed $spike_count spikes from $comp ($(round(spike_percentage, digits=4))% of data)")
            
            # Save spike statistics for reporting
            despike_stats[comp] = Dict(
                "spike_count" => spike_count,
                "spike_percentage" => spike_percentage
            )
            
            # Write the despiked ATS file with "Proc_" prefix
            if haskey(ats_files_map, comp)
                orig_file = ats_files_map[comp]
                orig_dir = dirname(orig_file)
                target_dir = dir_mapping[orig_dir]
                
                orig_filename = basename(orig_file)
                # Create new filename with Proc_ prefix
                new_filename = "Proc_" * orig_filename
                output_file = joinpath(target_dir, new_filename)
                
                # Write the ATS file
                write_ats_file(output_file, despiked, data[comp]["header"])
                println("Saved despiked $comp data to $output_file")
            end
        end
    end
    
    # Copy XML files preserving directory structure
    for xml_file in xml_files
        orig_dir = dirname(xml_file)
        if haskey(dir_mapping, orig_dir)
            target_dir = dir_mapping[orig_dir]
            dest_file = joinpath(target_dir, basename(xml_file))
            cp(xml_file, dest_file, force=true)
            println("Copied XML file to $dest_file")
        else
            # If the XML is in a directory we didn't process, put it in the root output dir
            dest_file = joinpath(output_dir, basename(xml_file))
            cp(xml_file, dest_file, force=true)
            println("Copied XML file to $dest_file")
        end
    end
    
    # Create a summary file with despiking statistics in the root output directory
    summary_file = joinpath(output_dir, "despiking_summary.txt")
    open(summary_file, "w") do file
        println(file, "# Despiking Summary for $(get(site_info, "SurveyID", "Unknown"))")
        println(file, "# Site: $(get(site_info, "SurveyID", "Unknown"))")
        println(file, "# ADU: $(get(site_info, "adu_number", "Unknown"))")
        println(file, "# Processing date: $(Dates.format(now(), "yyyy-mm-dd HH:MM:SS"))")
        println(file, "# ----------------------------------------")
        
        for comp in ["Ex", "Ey", "Hx", "Hy", "Hz"]
            if haskey(despike_stats, comp)
                stats = despike_stats[comp]
                println(file, "# $comp: $(stats["spike_count"]) spikes removed ($(round(stats["spike_percentage"], digits=4))%)")
            else
                println(file, "# $comp: Not processed")
            end
        end
    end
    
    return despiked_data, despike_stats
end

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

"""
    create_comparison_plot(original_data, despiked_data, site_info, despike_stats, output_file)
    
Creates a PDF with plots showing the original and despiked time series for comparison
"""
function create_comparison_plot(original_data, despiked_data, site_info, despike_stats, output_file)
    gr()
    comp_colors = Dict(
        "Ex" => :darkblue,
        "Ey" => :darkblue,
        "Hx" => :royalblue2,
        "Hy" => :royalblue2,
        "Hz" => :royalblue2
    )
    spike_color = :red
    
    components = ["Ex", "Ey", "Hx", "Hy", "Hz"]
    valid_components = filter(c -> haskey(original_data, c) && haskey(despiked_data, c), components)
    n_plots = length(valid_components)

    first_comp = valid_components[1]
    fs = original_data[first_comp]["fs"]
    total_samples = length(original_data[first_comp]["data"])
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
    else
        xticks = ([time_axis[1], time_axis[end]], ["0 h", @sprintf("%.2f h", duration_seconds/3600)])
        xlabel = "Elapsed Time (hours)"
    end

    comp_units = Dict("Ex"=>"mV/km", "Ey"=>"mV/km", "Hx"=>"nT", "Hy"=>"nT", "Hz"=>"nT")

    # Create two plots for each component - original and despiked
    plots = []
    
    for comp in valid_components
        # Original data (with spikes highlighted)
        org_data = original_data[comp]["data"][1:decimation_factor:end]
        desp_data = despiked_data[comp]["data"][1:decimation_factor:end]
        color = get(comp_colors, comp, :black)
        unit = get(comp_units, comp, "")
        
        # Get global min/max for consistent y-axis scaling
        all_data = vcat(org_data, desp_data)
        global_ylims = (minimum(all_data), maximum(all_data))
        
        # Create original plot
        p_org = plot(
            time_axis, org_data,
            color=color, linewidth=0.8, 
            title="Original $(comp)",
            ylabel=comp*" ("*unit*")",
            xlabel="",
            xticks=false,
            yticks=:auto,
            framestyle=:box,
            grid=false,
            minorgrid=false,
            left_margin=15Plots.mm,
            bottom_margin=5Plots.mm,
            xlims=(time_axis[1], time_axis[end]),
            ylims=global_ylims,
            tickfontsize=9, labelfontsize=11,
            legendfontsize=9,
            dpi=300,
            background_color=:white,
            legend=:topright
        )
        
        # Create despiked plot
        p_desp = plot(
            time_axis, desp_data,
            color=color, linewidth=0.8,
            title="Despiked $(comp)",
            ylabel=comp*" ("*unit*")",
            xlabel=(comp == valid_components[end] ? xlabel : ""),
            xticks=(comp == valid_components[end] ? xticks : false),
            yticks=:auto,
            framestyle=:box,
            grid=false,
            minorgrid=false,
            left_margin=15Plots.mm,
            bottom_margin=(comp == valid_components[end] ? 15Plots.mm : 5Plots.mm),
            xlims=(time_axis[1], time_axis[end]),
            ylims=global_ylims,
            tickfontsize=9, labelfontsize=11,
            legendfontsize=9,
            dpi=300,
            background_color=:white,
            legend=:topright
        )
        
        # Add spike percentage annotation if statistics are available
        if haskey(despike_stats, comp)
            spike_percentage = despike_stats[comp]["spike_percentage"]
            annotate!(p_org, [(0.98, 0.95, text("Spikes: $(round(spike_percentage, digits=2))%", 
                           9, :right, :top, :red))])
        end
        
        push!(plots, p_org)
        push!(plots, p_desp)
    end

    # Arrange plots in a grid with two columns (original and despiked)
    layout = @layout [grid(n_plots, 2)]
    
    final_plot = plot(plots...,
        layout=layout,
        size=(1000, 300*n_plots),
        margin=8Plots.mm,
        title="Time Series Comparison - Original vs Despiked"
    )
    
    savefig(final_plot, output_file)
    println("Comparison plot saved as $output_file")
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

    println("╔" * "═"^(key_width+2) * "╦" * "═"^(val_width+2) * "╗")

    for (key, val) in display_fields
        println("║ " * rpad(key, key_width) * " ║ " * rpad(val, val_width) * " ║")
    end

    println("╚" * "═"^(key_width+2) * "╩" * "═"^(val_width+2) * "╝")
end

# Main execution
println("Script starting at $(now())")
println("Current directory: $(pwd())")
println("Calling read_ts(128) to find and plot 128Hz data")

# Read the original data
data, site_info, ats_files_map, xml_files, target_dirs = read_ts(128)
print_site_info_table(site_info)

# Create output directory
output_dir = "DF090-W"
println("Output directory set to: $output_dir")

# Apply despiking and save the despiked data
despiked_data, despike_stats = despike_and_save_data(data, site_info, ats_files_map, xml_files, target_dirs, output_dir)

# Create comparison plot in the root output directory
comparison_plot_file = joinpath(output_dir, "timeseries_comparison.pdf")
create_comparison_plot(data, despiked_data, site_info, despike_stats, comparison_plot_file)

println("Despiking completed successfully.")
println("Despiked data written to $output_dir with preserved directory structure")
println("Comparison plot saved as $comparison_plot_file")