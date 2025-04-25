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