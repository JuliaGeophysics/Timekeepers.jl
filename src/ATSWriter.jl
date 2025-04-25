"""
    write_all_ats_files(output_dir::String, data, site_info, ats_files_map, xml_files)
    
Write all time series components to ATS files in the specified output directory,
preserving original filenames. Also copies XML files and creates a processing summary.

Parameters:
- output_dir: Directory where processed files will be saved
- data: Dictionary of time series data components
- site_info: Dictionary of site information
- ats_files_map: Mapping of component names to original file paths
- xml_files: List of XML files to copy
"""
function write_all_ats_files(output_dir::String, data, site_info, ats_files_map, xml_files)
    # Get current date/time in UTC
    current_datetime = Dates.format(now(Dates.UTC), "yyyy-mm-dd HH:MM:SS")
    current_user = get(ENV, "USER", get(ENV, "USERNAME", "unknown"))
    
    println("Current Date and Time (UTC): $current_datetime")
    println("Current User's Login: $current_user")
    
    # Create output directory if it doesn't exist
    if !isdir(output_dir)
        println("Creating output directory: $output_dir")
        mkpath(output_dir)
    end
    
    # Process each component
    processed_files = String[]
    for component in ["Ex", "Ey", "Hx", "Hy", "Hz"]
        if haskey(data, component) && haskey(ats_files_map, component)
            # Get the original file path
            original_file = ats_files_map[component]
            
            # Extract the base filename
            original_filename = basename(original_file)
            
            # Create a new filename (with same name)
            output_filename = joinpath(output_dir, original_filename)
            
            # Get the time series data and header
            ts_data = data[component]["data"]
            header = data[component]["header"]
            
            # Update header with current processing time if needed
            if haskey(header, "processing_datetime")
                header["processing_datetime"] = current_datetime
            end
            
            # Write the ATS file
            write_ats_file(output_filename, ts_data, header)
            println("Successfully wrote $component data to $output_filename")
            push!(processed_files, output_filename)
        else
            if !haskey(data, component)
                println("Component $component not found in data")
            elseif !haskey(ats_files_map, component)
                println("Original filename for $component not found")
            end
        end
    end
    
    # Copy XML files to the output directory
    for xml_file in xml_files
        xml_filename = basename(xml_file)
        output_xml = joinpath(output_dir, xml_filename)
        cp(xml_file, output_xml, force=true)
        println("Copied XML file: $xml_filename to $output_dir")
        push!(processed_files, output_xml)
    end
    
    
    return processed_files
end