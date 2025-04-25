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