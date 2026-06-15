using Dates
using Printf
using Test
using TimeSeries
using Timekeepers

function _small_timearray(; n = 12)
    t0 = DateTime(2020, 1, 1)
    times = [t0 + Second(i - 1) for i in 1:n]
    vals = Matrix{Float64}(undef, n, 2)
    @inbounds for i in 1:n
        vals[i, 1] = i
        vals[i, 2] = 2i
    end
    return TimeArray(times, vals, [:bx, :by], Dict{Symbol, Any}(:sample_rate => 1.0))
end

function _write_sample_lemi424(path::AbstractString; n = 4, extra_columns = 0, drop_trailing = 0, start::DateTime = DateTime(2020, 1, 1))
    open(path, "w") do io
        for i in 0:(n - 1)
            t = start + Second(i)
            fields = Any[
                Dates.year(t),
                lpad(Dates.month(t), 2, '0'),
                lpad(Dates.day(t), 2, '0'),
                lpad(Dates.hour(t), 2, '0'),
                lpad(Dates.minute(t), 2, '0'),
                lpad(Dates.second(t), 2, '0'),
                1.0 + i,
                2.0 + i,
                3.0 + i,
                10.0,
                11.0,
                4.0 + i,
                5.0 + i,
                0.0,
                0.0,
                12.5,
                100.0,
                6022.00000,
                "N",
                02456.00000,
                "E",
                8,
                1,
                0,
            ]
            for extra in 1:extra_columns
                push!(fields, 1000 + extra)
            end
            if drop_trailing > 0
                resize!(fields, length(fields) - drop_trailing)
            end
            println(io, join(fields, ' '))
        end
    end
    return path
end

function _write_sample_geomag(path::AbstractString; n = 6)
    open(path, "w") do io
        println(io, "; MS:GEOMAG-02  #26-2011")
        println(io, "; Date: 2025/05/23; Time: 00:00:00")
        println(io, "; Sampling: 0.10 sec")
        println(io, "; Latitude: 60 35'14.4\"N;  Longitude: 027 35'09.0\"E;  Altitude: ------")
        println(io, ";    Date        Time    X [nT]   Y [nT]   Z [nT]    Ex[mV]   Ey[mV]   Ts[C] Te[C]")
        println(io, ";")
        for i in 0:(n - 1)
            @printf(
                io,
                "2025 05 23  00 00 %05.2f %+09.3f %+09.3f %+09.3f %+08.3f %+08.3f %+04.1f %+04.1f\n",
                i / 10,
                44.0 + i,
                33.0 + i,
                10.0 + i,
                15.0 + i,
                -7.0 - i,
                9.4,
                18.5,
            )
        end
    end
    return path
end

const _METRONIX_CHANNELS = [
    ("C00", "TEx", "Ex"),
    ("C01", "TEy", "Ey"),
    ("C02", "THx", "Hx"),
    ("C03", "THy", "Hy"),
    ("C04", "THz", "Hz"),
]

function _write_sample_ats(path::AbstractString, channel_type::AbstractString,
                           data::AbstractVector{<:Real}, lsb::Float64, fs::Real, start_unix::Integer)
    hdr = zeros(UInt8, 1024)
    hdr[1:2] = reinterpret(UInt8, [UInt16(1024)])
    hdr[5:8] = reinterpret(UInt8, [Int32(length(data))])
    hdr[9:12] = reinterpret(UInt8, [Float32(fs)])
    hdr[13:16] = reinterpret(UInt8, [Int32(start_unix)])
    hdr[17:24] = reinterpret(UInt8, [Float64(lsb)])
    ctb = Vector{UInt8}(codeunits(channel_type))
    hdr[39:(38 + length(ctb))] = ctb
    raw = round.(Int32, data ./ lsb)
    open(path, "w") do io
        write(io, hdr)
        write(io, raw)
    end
    return path
end

function _write_sample_metronix(meas_dir::AbstractString; n = 80, fs = 8,
                                start_dt::DateTime = DateTime(2025, 4, 1, 7, 0, 6))
    mkpath(meas_dir)
    token = "$(round(Int, fs))H"
    start_unix = round(Int, Dates.datetime2unix(start_dt))
    stop_dt = start_dt + Second((n - 1) ÷ fs)
    data = Dict{String, Vector{Float64}}()
    chan_xml = IOBuffer()
    for (k, (cc, tag, ct)) in enumerate(_METRONIX_CHANNELS)
        lsb = 1.0e-6 * k
        d = [sin(i / 5) + k for i in 1:n]
        fname = "999_V01_$(cc)_R000_$(tag)_BL_$(token).ats"
        _write_sample_ats(joinpath(meas_dir, fname), ct, d, lsb, fs, start_unix)
        data[ct] = d
        print(chan_xml, """
              <channel id="$(k - 1)">
                <start_time>$(Dates.format(start_dt, "HH:MM:SS"))</start_time>
                <start_date>$(Dates.format(start_dt, "yyyy-mm-dd"))</start_date>
                <num_samples>$n</num_samples>
                <ats_data_file>$fname</ats_data_file>
              </channel>
""")
    end
    fmt(dt) = Dates.format(dt, "yyyy-mm-dd_HH-MM-SS")
    xml_name = "999_$(fmt(start_dt))_$(fmt(stop_dt))_R000_$(token).xml"
    xml = """
<?xml version='1.0' encoding='UTF-8' standalone='no'?>
<measurement>
  <recording>
    <start_time>$(Dates.format(start_dt, "HH:MM:SS"))</start_time>
    <stop_time>$(Dates.format(stop_dt, "HH:MM:SS"))</stop_time>
    <start_date>$(Dates.format(start_dt, "yyyy-mm-dd"))</start_date>
    <stop_date>$(Dates.format(stop_dt, "yyyy-mm-dd"))</stop_date>
    <output>
      <ProcessingTree>
        <output>
          <DigitalFilter>
            <output>
              <ATSWriter>
                <configuration>
$(String(take!(chan_xml)))                </configuration>
                <output_file>
                  <ats_file_size>$(1024 + n * 4)</ats_file_size>
                </output_file>
              </ATSWriter>
            </output>
          </DigitalFilter>
        </output>
      </ProcessingTree>
    </output>
  </recording>
</measurement>
"""
    open(joinpath(meas_dir, xml_name), "w") do io
        print(io, xml)
    end
    return data
end

@testset "Metronix ATS read and round-trip" begin
    mktempdir() do root
        meas = joinpath(root, "RK999", "meas_2025-04-01_07-00-05")
        truth = _write_sample_metronix(meas)
        run = read_metronix(meas)

        @test components(run) == [:bx, :by, :bz, :e1, :e2]
        @test isapprox(sampling_rate(run), 8.0)
        @test start_time(run) == DateTime(2025, 4, 1, 7, 0, 6)
        @test run.metadata[:n_samples] == 80
        @test maximum(abs.(run.channels[:e1].data .- truth["Ex"])) < 1e-5
        @test maximum(abs.(run.channels[:bz].data .- truth["Hz"])) < 1e-5

        out = joinpath(root, "single", "meas_out")
        write_metronix(out, run)
        @test count(f -> endswith(f, ".ats"), readdir(out)) == 5
        @test count(f -> endswith(f, ".xml"), readdir(out)) == 1
        rt = read_metronix(out)
        @test rt.channels[:e1].data == run.channels[:e1].data
        @test rt.channels[:by].data == run.channels[:by].data
    end
end

@testset "Metronix amputation writes split meas_ dirs" begin
    mktempdir() do root
        site = joinpath(root, "RK999")
        meas = joinpath(site, "meas_2025-04-01_07-00-05")
        _write_sample_metronix(meas)
        run = read_metronix(meas)

        ta = to_timearray(run)
        mask = TimekeeperMask(ta)
        mask.masked[33:48] .= true

        dest, dirs = write_metronix_site(run; mask = mask)
        @test dest == site * ".W"
        @test isdir(dest)
        @test length(dirs) == 2

        for d in dirs
            @test count(f -> endswith(f, ".ats"), readdir(d)) == 5
            @test count(f -> endswith(f, ".xml"), readdir(d)) == 1
        end

        seg1 = read_metronix(dirs[1])
        seg2 = read_metronix(dirs[2])
        @test length(seg1.channels[:e1].data) == 32
        @test length(seg2.channels[:e1].data) == 32
        @test seg1.channels[:e1].data == run.channels[:e1].data[1:32]
        @test seg2.channels[:e1].data == run.channels[:e1].data[49:80]
        @test start_time(seg1) == DateTime(2025, 4, 1, 7, 0, 6)
        @test start_time(seg2) == DateTime(2025, 4, 1, 7, 0, 12)

        xml2 = only(filter(f -> endswith(f, ".xml"), readdir(dirs[2])))
        @test occursin("07-00-12", xml2)
        xml_text = read(joinpath(dirs[2], xml2), String)
        @test occursin("<num_samples>32</num_samples>", xml_text)
        @test basename(dirs[2]) == "meas_2025-04-01_07-00-12"
    end
end

@testset "Metronix split with two gaps yields three dirs" begin
    mktempdir() do root
        meas = joinpath(root, "RK998", "meas_2025-04-01_07-00-05")
        _write_sample_metronix(meas)
        run = read_metronix(meas)
        ta = to_timearray(run)
        mask = TimekeeperMask(ta)
        mask.masked[17:24] .= true
        mask.masked[49:56] .= true
        _, dirs = write_metronix_site(run; mask = mask)
        @test length(dirs) == 3
    end
end

@testset "Metronix split-by-rate helper script" begin
    include(joinpath(@__DIR__, "..", "scripts", "split_metronix_by_rate.jl"))
    mktempdir() do root
        site = joinpath(root, "RKMULTI")
        _write_sample_metronix(joinpath(site, "meas_2025-04-01_07-00-05"); fs = 8, n = 80)
        _write_sample_metronix(joinpath(site, "meas_2025-04-01_08-00-05"); fs = 8, n = 80,
                               start_dt = DateTime(2025, 4, 1, 8, 0, 6))
        _write_sample_metronix(joinpath(site, "meas_2025-04-01_09-00-05"); fs = 64, n = 640,
                               start_dt = DateTime(2025, 4, 1, 9, 0, 6))

        @test metronix_site_rates(site) == [8.0, 64.0]
        runs = metronix_site_runs(site)
        @test length(runs[8.0]) == 2
        @test length(runs[64.0]) == 1

        dests = split_metronix_site_by_rate(site)
        @test Set(basename.(dests)) == Set(["RKMULTI.TK8", "RKMULTI.TK64"])
        d8 = joinpath(root, "RKMULTI.TK8")
        d64 = joinpath(root, "RKMULTI.TK64")
        @test count(x -> startswith(x, "meas_"), readdir(d8)) == 2
        @test count(x -> startswith(x, "meas_"), readdir(d64)) == 1
        @test metronix_site_rates(d8) == [8.0]
        @test metronix_site_rates(d64) == [64.0]

        src = joinpath(site, "meas_2025-04-01_09-00-05")
        dst = joinpath(d64, "meas_2025-04-01_09-00-05")
        for f in readdir(src)
            @test read(joinpath(src, f)) == read(joinpath(dst, f))
        end
    end
end

@testset "Metronix masked site write round-trips losslessly" begin
    mktempdir() do root
        site = joinpath(root, "RK128")
        meas = joinpath(site, "meas_2025-04-01_07-00-05")
        _write_sample_metronix(meas; fs = 8, n = 80)

        dest = write_metronix_site_masked(site)
        @test dest == site * ".W"
        @test isfile(joinpath(dest, "README.md"))
        out_meas = joinpath(dest, only(filter(x -> startswith(x, "meas_"), readdir(dest))))
        for f in filter(x -> endswith(x, ".ats"), readdir(meas))
            @test read(joinpath(meas, f)) == read(joinpath(out_meas, f))
        end
        a = read_metronix(meas)
        b = read_metronix(out_meas)
        @test a.channels[:e1].data == b.channels[:e1].data
        @test a.channels[:bz].data == b.channels[:bz].data
        @test start_time(a) == start_time(b)
    end
end

@testset "Metronix write preserves rate tag and logs mask history" begin
    mktempdir() do root
        site = joinpath(root, "RK137.TK128")
        meas = joinpath(site, "meas_2025-04-01_07-00-05")
        _write_sample_metronix(meas; fs = 8, n = 80)

        iv = (DateTime(2025, 4, 1, 7, 0, 8), DateTime(2025, 4, 1, 7, 0, 10))
        dest = write_metronix_site_masked(site; intervals = [iv])
        @test dest == site * ".W"
        @test basename(dest) == "RK137.TK128.W"

        readme = joinpath(dest, "README.md")
        @test isfile(readme)
        text = read(readme, String)
        @test occursin("Write session", text)
        @test occursin(string(iv[1]), text)
        @test occursin(string(iv[2]), text)
        @test occursin("RK137.TK128", text)

        write_metronix_site_masked(site; intervals = [iv])
        text2 = read(readme, String)
        @test length(collect(eachmatch(r"## Write session", text2))) == 2
    end
end

@testset "masking and stream writes" begin
    ta = _small_timearray()
    mask = TimekeeperMask(ta)
    times = Timekeepers._ta_timestamps(ta)
    mask_interval!(mask, times[3], times[5])

    cleaned = cleaned_timearray(ta, mask)
    @test all(isnan, Timekeepers._ta_values(cleaned)[3:5, :])
    dropped = cleaned_timearray(ta, mask; mode = :drop)
    @test length(Timekeepers._ta_timestamps(dropped)) == 9
    @test sample_weights(mask; good = 1, bad = 0)[3:5] == [0, 0, 0]

    mktempdir() do dir
        mask_path = joinpath(dir, "mask.csv")
        data_path = joinpath(dir, "cleaned.csv")
        write_mask(mask_path, mask)
        write_cleaned(data_path, ta, mask)
        round_trip = read_mask(mask_path, ta)
        @test round_trip.masked == mask.masked
        @test occursin("timestamp,bx,by", read(data_path, String))
    end
end

@testset "LEMI-424 flexible columns" begin
    mktempdir() do dir
        extra_path = joinpath(dir, "sample_extra.txt")
        short_path = joinpath(dir, "sample_short.txt")
        _write_sample_lemi424(extra_path; extra_columns = 2)
        _write_sample_lemi424(short_path; drop_trailing = 2)

        extra = load_lemi424(extra_path)
        short = load_lemi424(short_path)
        @test size(Timekeepers._ta_values(extra)) == (4, 5)
        @test size(Timekeepers._ta_values(short)) == (4, 5)
        @test Timekeepers._ta_values(extra)[:, 1:3] == Timekeepers._ta_values(short)[:, 1:3]

        run = read_lemi424(short_path)
        @test :bx in components(run)
        @test isnan(run.metadata[:elevation]) == false
        @test run.metadata[:battery_start] == 12.5
        @test all(isnan, run.channels[:time_diff].data)
    end
end

@testset "site load (directory of LEMI-424 files)" begin
    mktempdir() do dir
        # Three runs: 4 s each, with a 6 s gap between files 1 and 2, contiguous between 2 and 3.
        _write_sample_lemi424(joinpath(dir, "run_a.txt"); n = 4, start = DateTime(2020, 1, 1, 0, 0, 0))
        _write_sample_lemi424(joinpath(dir, "run_b.txt"); n = 4, start = DateTime(2020, 1, 1, 0, 0, 10))
        _write_sample_lemi424(joinpath(dir, "run_c.txt"); n = 4, start = DateTime(2020, 1, 1, 0, 0, 14))
        # Non-data file should be ignored.
        write(joinpath(dir, "notes.md"), "ignore me")

        ta, fmt = Timekeepers._load_site_directory(dir)
        @test fmt == :lemi424
        times = Timekeepers._ta_timestamps(ta)
        vals = Timekeepers._ta_values(ta)
        @test first(times) == DateTime(2020, 1, 1, 0, 0, 0)
        @test last(times) == DateTime(2020, 1, 1, 0, 0, 17)
        @test size(vals, 1) == 18
        @test Timekeepers._ta_colnames(ta) == [:bx, :by, :bz, :e1, :e2]
        # Samples 5..10 sit in the gap (indices 5..10 -> seconds 4..9 inclusive).
        @test all(isnan, vals[5:10, :])
        # First and last data points should be finite.
        @test all(isfinite, vals[1, :])
        @test all(isfinite, vals[end, :])
        meta = Timekeepers._ta_meta(ta)
        @test meta[:n_files] == 3
        @test meta[:sample_rate] == 1.0
        @test meta[:site] == basename(dir)

        # Auto-write combined SITENAME.txt and verify re-scan ignores it.
        out_path = Timekeepers._write_combined_site!(ta, dir, fmt)
        @test isfile(out_path)
        @test basename(out_path) == basename(dir) * ".txt"
        files_seen = Timekeepers._list_data_files(dir)
        @test out_path ∉ files_seen
        @test length(files_seen) == 3

        ta2, fmt2 = Timekeepers._load_site_directory(dir)
        @test fmt2 == :lemi424
        @test Timekeepers._ta_meta(ta2)[:n_files] == 3
    end
end

@testset "LEMI-424 IO" begin
    mktempdir() do dir
        in_path = joinpath(dir, "sample.txt")
        out_path = joinpath(dir, "sample_out.txt")
        _write_sample_lemi424(in_path)

        ta = load_lemi424(in_path)
        @test size(Timekeepers._ta_values(ta)) == (4, 5)
        @test Timekeepers._ta_colnames(ta) == [:bx, :by, :bz, :e1, :e2]
        run = read_lemi424(in_path; include_aux = false)
        @test sampling_rate(run) == 1.0
        @test :bx in components(run)

        write_lemi424(out_path, ta)
        ta2 = load_lemi424(out_path)
        @test size(Timekeepers._ta_values(ta2)) == (4, 5)
    end
end

@testset "GEOMAG IO" begin
    mktempdir() do dir
        path = joinpath(dir, "GEOMAG.TXT")
        _write_sample_geomag(path)

        @test Timekeepers._detect_format(path) == :geomag
        ta = load_geomag(path)
        @test size(Timekeepers._ta_values(ta)) == (6, 5)
        @test Timekeepers._ta_colnames(ta) == [:bx, :by, :bz, :e1, :e2]
        @test Timekeepers._sample_rate_from_timearray(ta) == 10.0
        @test Timekeepers._ta_meta(ta)[:instrument_model] == "GEOMAG-02"

        run = read_timekeeper(path)
        @test run.source_format == :geomag
        @test sampling_rate(run) == 10.0
        @test :temperature_e in components(run)
        @test :temperature_h in components(run)
    end
end

@testset "LEMI-424 aux preservation" begin
    mktempdir() do dir
        in_path = joinpath(dir, "sample.txt")
        out_path = joinpath(dir, "sample_out.txt")
        _write_sample_lemi424(in_path)

        ta = load_lemi424(in_path)
        aux = Timekeepers._ta_meta(ta)[:aux_columns]
        @test aux[:temperature_e] == fill(10.0, 4)
        @test aux[:temperature_h] == fill(11.0, 4)
        @test aux[:battery] == fill(12.5, 4)
        @test aux[:elevation] == fill(100.0, 4)
        @test aux[:lat_hemisphere] == fill("N", 4)
        @test aux[:lon_hemisphere] == fill("E", 4)
        @test aux[:n_satellites] == fill(8.0, 4)
        @test aux[:gps_fix] == fill(1.0, 4)

        write_lemi424(out_path, ta)
        round_trip = load_lemi424(out_path)
        ra = Timekeepers._ta_meta(round_trip)[:aux_columns]
        @test ra[:temperature_e] == fill(10.0, 4)
        @test ra[:temperature_h] == fill(11.0, 4)
        @test ra[:battery] == fill(12.5, 4)
        @test ra[:elevation] == fill(100.0, 4)
        @test ra[:lat_hemisphere] == fill("N", 4)
        @test ra[:lon_hemisphere] == fill("E", 4)
        @test ra[:n_satellites] == fill(8.0, 4)
        @test ra[:gps_fix] == fill(1.0, 4)
    end
end

@testset "GEOMAG aux preservation" begin
    mktempdir() do dir
        in_path = joinpath(dir, "GEOMAG.TXT")
        out_path = joinpath(dir, "GEOMAG_clean.TXT")
        _write_sample_geomag(in_path)
        ta = load_geomag(in_path)
        aux = Timekeepers._ta_meta(ta)[:aux_columns]
        @test all(aux[:temperature_h] .== 9.4)
        @test all(aux[:temperature_e] .== 18.5)

        write_geomag(out_path, ta)
        round_trip = load_geomag(out_path)
        ra = Timekeepers._ta_meta(round_trip)[:aux_columns]
        @test all(ra[:temperature_h] .== 9.4)
        @test all(ra[:temperature_e] .== 18.5)
    end
end

@testset "GEOMAG writer" begin
    mktempdir() do dir
        in_path = joinpath(dir, "GEOMAG.TXT")
        out_path = joinpath(dir, "GEOMAG_clean.TXT")
        _write_sample_geomag(in_path)
        ta = load_geomag(in_path)
        write_geomag(out_path, ta)
        @test Timekeepers._detect_format(out_path) == :geomag
        round_trip = load_geomag(out_path)
        @test size(Timekeepers._ta_values(round_trip)) == (6, 5)
    end
end

@testset "spectral workspace" begin
    n = 512
    x = [sin(2π * (i - 1) / 32) for i in 1:n]
    y = [cos(2π * (i - 1) / 32) for i in 1:n]
    masked = falses(n)
    ws = Timekeepers.SpectralWorkspace(128, 1.0)

    freqs, psd = Timekeepers._welch_psd(x, 1.0; nfft = 128, workspace = ws)
    @test length(freqs) == 65
    @test length(psd) == 65
    @test all(isfinite, psd)

    segs = [view(x, 1:256), view(x, 257:512)]
    freqs2, psd2, n_used = Timekeepers._welch_psd_segments(segs, 1.0; nfft = 128, workspace = ws)
    @test freqs2 == freqs
    @test length(psd2) == 65
    @test n_used == 2

    sfreqs, stimes, spec = Timekeepers._stft_psd(x, masked, 1.0; nfft = 128, workspace = ws)
    @test length(sfreqs) == 65
    @test size(spec, 1) == 65
    @test size(spec, 2) == length(stimes)
end

@testset "app icon" begin
    icons = Timekeepers._timekeepers_icons()
    @test size.(icons, 1) == [16, 32, 64, 128]
    @test size.(icons, 2) == [16, 32, 64, 128]
    @test isfile(Timekeepers.TIMEKEEPERS_LOGO_PATH)

    mktempdir() do dir
        path = joinpath(dir, "timekeepers-icon.png")
        Timekeepers._write_timekeepers_icon_png(path; icon_size = 32)
        @test isfile(path)
        @test read(path, 8) == UInt8[0x89, 0x50, 0x4e, 0x47, 0x0d, 0x0a, 0x1a, 0x0a]
    end
end

@testset "startup status" begin
    ta = _small_timearray()
    app = TKApp(ta; size = (700, 420))
    @test occursin("Timekeepers ready", app.status_label.text[])
end
