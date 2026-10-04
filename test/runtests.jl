# runtests.jl - the package test suite.
# Author: @pankajkmishra
#
# The tests examine each format from end to end. Each test writes a synthetic
# file, reads it again and compares the result. The formats are:
# - LEMI-424, also files with more or fewer columns at the end
# - GEOMAG
# - Metronix ATS, with the cuts in a site and meas_ directories that hold more
#   than one run at more than one rate
# The tests also examine:
# - masks, clean output and segments
# - the load of a site with many files, with gaps filled
# - spectral estimation
# - the app: it builds and reports a ready status
# - the survey scan and the TKDash dashboard

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

function _write_sample_metronix(meas_dir::AbstractString; n = 80, fs = 8, run = "R000",
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
        fname = "999_V01_$(cc)_$(run)_$(tag)_BL_$(token).ats"
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
    xml_name = "999_$(fmt(start_dt))_$(fmt(stop_dt))_$(run)_$(token).xml"
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

# each file in root, by its path relative to root
_tree(root) = Dict(replace(relpath(joinpath(d, f), root), '\\' => '/') => joinpath(d, f)
                   for (d, _, files) in walkdir(root) for f in files)

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

@testset "Metronix run cut at a mask, in its rate directory" begin
    mktempdir() do root
        site = joinpath(root, "RK999")
        meas = joinpath(site, "meas_2025-04-01_07-00-05")
        _write_sample_metronix(meas)
        before = Dict(k => read(p) for (k, p) in _tree(site))
        run = read_metronix(meas)

        ta = to_timearray(run)
        mask = TimekeeperMask(ta)
        mask.masked[33:48] .= true

        dest, ids = write_metronix_site(run; mask = mask)
        @test dest == Timekeepers._norm_path(joinpath(site * ".8", basename(meas)))
        @test Dict(k => read(p) for (k, p) in _tree(site)) == before      # the site never changes
        @test length(ids) == 2
        @test count(f -> endswith(f, ".ats"), readdir(dest)) == 10
        @test count(f -> endswith(f, ".xml"), readdir(dest)) == 2

        seg1 = read_metronix(ids[1])
        seg2 = read_metronix(ids[2])
        @test seg1.metadata[:metronix_run_token] == "R000"
        @test seg2.metadata[:metronix_run_token] == "R001"
        @test seg1.channels[:e1].data == run.channels[:e1].data[1:32]
        @test seg2.channels[:e1].data == run.channels[:e1].data[49:80]
        @test start_time(seg1) == DateTime(2025, 4, 1, 7, 0, 6)
        @test start_time(seg2) == DateTime(2025, 4, 1, 7, 0, 12)

        # 32 samples at 8 Hz from 07:00:12 record 4 s. The stop of the ADU is 07:00:16
        @test basename(ids[2]) == "999_2025-04-01_07-00-12_2025-04-01_07-00-16_R001_8H.xml"
        xml_text = read(ids[2], String)
        @test occursin("<num_samples>32</num_samples>", xml_text)

        # the XML of the segment is the XML of the run. Only the fields of the
        # segment change
        template = read(run.metadata[:metronix_xml_path], String)
        a, b = split(template, '\n'), split(xml_text, '\n')
        @test length(a) == length(b)
        changed = [strip(b[i]) for i in eachindex(a) if a[i] != b[i]]
        @test all(l -> occursin(r"^<(start_time|stop_time|start_date|stop_date|num_samples|ats_file_size|ats_data_file)>", l), changed)
        @test "<start_time>07:00:12</start_time>" in changed
        @test "<stop_time>07:00:16</stop_time>" in changed
        @test "<ats_file_size>$(1024 + 32 * 4)</ats_file_size>" in changed
        @test "<ats_data_file>999_V01_C00_R001_TEx_BL_8H.ats</ats_data_file>" in changed

        # the cut of a run from its rate directory occurs in that location
        mask1 = TimekeeperMask(to_timearray(seg1))
        mask1.masked[9:16] .= true
        dest2, ids2 = write_metronix_site(seg1; mask = mask1)
        @test dest2 == dest
        @test [read_metronix(id).metadata[:metronix_run_token] for id in ids2] == ["R000", "R002"]
        @test length(metronix_site_runs(dest)[8.0]) == 3
        @test !isfile(ids[1])                             # the cut run is removed
        @test !any(startswith(".tk_staging_"), readdir(dirname(dest)))
    end
end

@testset "Metronix run cut at two gaps yields three runs" begin
    mktempdir() do root
        meas = joinpath(root, "RK998", "meas_2025-04-01_07-00-05")
        _write_sample_metronix(meas)
        run = read_metronix(meas)
        ta = to_timearray(run)
        mask = TimekeeperMask(ta)
        mask.masked[17:24] .= true
        mask.masked[49:56] .= true
        _, ids = write_metronix_site(run; mask = mask)
        @test [read_metronix(id).metadata[:metronix_run_token] for id in ids] == ["R000", "R001", "R002"]
    end
end

@testset "Metronix skips segments too short to store as a run" begin
    mktempdir() do root
        meas = joinpath(root, "RK997", "meas_2025-04-01_07-00-05")
        _write_sample_metronix(meas)
        run = read_metronix(meas)
        mask = TimekeeperMask(to_timearray(run))
        # 8 Hz: samples 17:21 are 5 good samples, less than one whole second
        mask.masked[9:16] .= true
        mask.masked[22:40] .= true
        _, ids = @test_logs (:warn, r"too short to store") match_mode = :any write_metronix_site(run; mask = mask)
        @test start_time.(read_metronix.(ids)) == [DateTime(2025, 4, 1, 7, 0, 6), DateTime(2025, 4, 1, 7, 0, 11)]

        # min_samples sets a higher limit: the first 8 samples are also removed
        other = joinpath(root, "other", basename(meas))
        _, ids = @test_logs (:warn, r"too short to store") match_mode = :any write_metronix_site(run; mask = mask, min_samples = 16, dest = other)
        @test start_time(read_metronix(only(ids))) == DateTime(2025, 4, 1, 7, 0, 11)
    end
end

@testset "Metronix run without its XML reads and writes" begin
    mktempdir() do root
        meas = joinpath(root, "RK996", "meas_2025-04-01_07-00-05")
        _write_sample_metronix(meas)
        rm(only(filter(f -> endswith(f, ".xml"), readdir(meas; join = true))))
        run = @test_logs (:warn, r"No \.xml") match_mode = :any read_metronix(meas)
        @test run.metadata[:metronix_xml_path] === nothing
        @test run.metadata[:metronix_prefix] == "999"
        @test run.metadata[:metronix_run_token] == "R000"
        @test run.metadata[:metronix_freq_token] == "8H"
        @test run.metadata[:n_samples] == 80

        mask = TimekeeperMask(to_timearray(run))
        mask.masked[33:48] .= true
        dest, ids = write_metronix_site(run; mask = mask)
        @test length(ids) == 2
        @test count(f -> endswith(f, ".ats"), readdir(dest)) == 10
        @test !any(f -> endswith(f, ".xml"), readdir(dest))
        @test read_metronix(ids[2]).channels[:e1].data == run.channels[:e1].data[49:80]
    end
end

@testset "Metronix rate tokens" begin
    @test Timekeepers._freq_token(4096.0) == "4096H"
    @test Timekeepers._freq_token(128) == "128H"
    @test Timekeepers._freq_token(0.125) == "8S"
end

# A site with the shape of a real ADU campaign:
# - one meas_ directory with one high-rate run
# - one meas_ directory with a long low-rate run, two bursts at a high rate
#   (R000, R001), the XML of a burst that did not record (R002) and a .kml
const _BURST_N = 40960                                     # 20 s at 2048 Hz

function _write_sample_mixed_site(site::AbstractString)
    a = joinpath(site, "meas_2025-04-01_06-00-00")
    b = joinpath(site, "meas_2025-04-01_07-00-05")
    _write_sample_metronix(a; fs = 64, n = 640, start_dt = DateTime(2025, 4, 1, 6, 0, 0))
    # the slow run records without a stop through the bursts, as a 128 Hz run
    # does
    slow = _write_sample_metronix(b; fs = 8, n = 2400, start_dt = DateTime(2025, 4, 1, 7, 0, 6))
    burst0 = _write_sample_metronix(b; fs = 2048, n = _BURST_N, run = "R000",
                                    start_dt = DateTime(2025, 4, 1, 7, 1, 0))
    burst1 = _write_sample_metronix(b; fs = 2048, n = _BURST_N, run = "R001",
                                    start_dt = DateTime(2025, 4, 1, 7, 3, 0))
    write(joinpath(b, "999_2025-04-01_07-05-00_2025-04-01_07-05-00_R002_2048H.xml"),
          "<?xml version='1.0'?><measurement><recording/></measurement>")
    write(joinpath(b, "Site_meas_2025-04-01_07-00-05.kml"), "<kml/>")
    return (; a, b, slow, burst0, burst1)
end

@testset "Metronix meas_ directory holding several runs" begin
    mktempdir() do root
        site = joinpath(root, "DF999")
        f = _write_sample_mixed_site(site)

        @test is_metronix_site(site)
        @test metronix_site_rates(site) == [8.0, 64.0, 2048.0]
        runs = metronix_site_runs(site)
        @test length(runs[8.0]) == 1 && length(runs[64.0]) == 1
        @test occursin("_R000_2048H", basename(runs[2048.0][1]))   # in order of start time
        @test occursin("_R001_2048H", basename(runs[2048.0][2]))
        _, empty_xmls = Timekeepers._metronix_site_index(site)
        @test basename.(empty_xmls) == ["999_2025-04-01_07-05-00_2025-04-01_07-05-00_R002_2048H.xml"]

        # one meas_ directory has an index as a separate site
        @test metronix_site_rates(f.b) == [8.0, 2048.0]

        # an XML identifies its run. Each .ats file of the run also identifies it
        r1 = read_metronix(runs[2048.0][2])
        @test sampling_rate(r1) == 2048.0
        @test maximum(abs.(r1.channels[:e1].data .- f.burst1["Ex"])) < 1e-5
        @test start_time(r1) == DateTime(2025, 4, 1, 7, 3, 0)
        r1b = read_metronix(joinpath(f.b, "999_V01_C02_R001_THx_BL_2048H.ats"))
        @test r1b.metadata[:metronix_xml_path] == r1.metadata[:metronix_xml_path]
        slow = read_metronix(only(runs[8.0]))
        @test maximum(abs.(slow.channels[:bz].data .- f.slow["Hz"])) < 1e-5
        @test length(slow.channels[:bz].data) == 2400

        @test_throws ErrorException read_metronix(f.b)        # which one of the three runs?
        @test_throws ErrorException read_metronix(only(empty_xmls))
        @test sampling_rate(read_metronix(f.a)) == 64.0       # one run: the directory is sufficient

        # above 1 kHz, the bursts are joined end to end, and all the samples stay
        ta = load_metronix_site(site; rate = 2048)
        @test length(timestamp(ta)) == 2 * _BURST_N
        @test issorted(timestamp(ta))
        @test maximum(abs.(values(ta[:e1]) .- vcat(f.burst0["Ex"], f.burst1["Ex"]))) < 1e-5
        @test Timekeepers._sample_rate_from_timearray(ta) == 2048.0
        # read from the rate directory DF999.2048, which the load made
        @test Timekeepers._ta_meta(ta)[:metronix_runs] == metronix_site_runs(site * ".2048")[2048.0]
        @test Timekeepers._ta_meta(ta)[:site_dir] == Timekeepers._norm_path(site)
        @test Timekeepers._ta_meta(load_metronix_site(site * ".8"))[:site_dir] == Timekeepers._norm_path(site)

        ta_b, _ = Timekeepers._load_site_any(f.b; rate = 8.0)
        @test Timekeepers._ta_meta(ta_b)[:site_dir] == Timekeepers._norm_path(site)
    end
end

# The rate directories put together again into one site: each file of the
# original site goes to all its copies
function _unsplit(roots)
    m = Dict{String, Vector{String}}()
    for root in roots, (rel, p) in _tree(root)
        push!(get!(m, rel, String[]), p)
    end
    return m
end

function _reproduces(site, roots)
    orig, m = _tree(site), _unsplit(roots)
    return Set(keys(orig)) == Set(keys(m)) &&
           all(read(p) == read(orig[k]) for (k, ps) in m for p in ps)
end

_snapshot(dir) = Dict(k => read(p) for (k, p) in _tree(dir))

# Write the name of the XML of each burst in its .ats headers, as an ADU-07
# does
function _name_xml_in_headers!(meas)
    for p in readdir(meas; join = true)
        occursin(r"_R000_.*_2048H\.ats$", p) || continue
        bytes = read(p)
        xml = only(filter(n -> occursin(r"_R000_2048H\.xml$", n), readdir(meas)))
        bytes[449:(448 + length(xml))] = codeunits(xml)
        bytes[(449 + length(xml)):512] .= 0x00
        write(p, bytes)
    end
end

@testset "Metronix site separated by sampling rate" begin
    mktempdir() do root
        site = joinpath(root, "DF999")
        f = _write_sample_mixed_site(site)
        @test !metronix_site_is_split(site)
        rate_dirs = split_metronix_site(site)
        @test rate_dirs == Timekeepers._norm_path(site) .* [".8", ".64", ".2048"]
        @test sort(readdir(root)) == ["DF999", "DF999.2048", "DF999.64", "DF999.8"]
        @test metronix_site_rates(site * ".2048") == [2048.0]
        @test metronix_site_rates(site * ".8") == [8.0]
        @test Timekeepers._metronix_split_source(site * ".2048") == Timekeepers._norm_path(site)
        @test Timekeepers._metronix_split_source(site) === nothing

        b = basename(f.b)
        names(r) = sort(readdir(joinpath(site * "." * r, b)))
        @test all(n -> occursin("_2048H", n) || endswith(n, ".kml"), names("2048"))
        @test all(n -> occursin("_8H", n) || endswith(n, ".kml"), names("8"))
        @test "999_2025-04-01_07-05-00_2025-04-01_07-05-00_R002_2048H.xml" in names("2048")  # did not record
        @test "Site_meas_2025-04-01_07-00-05.kml" in names("8")
        @test "Site_meas_2025-04-01_07-00-05.kml" in names("2048")
        @test _reproduces(site, rate_dirs)                    # each file, byte for byte

        # the split never copies over a rate directory that exists, because it
        # can hold cuts
        @test metronix_site_is_split(site)
        xml8 = joinpath(site * ".8", b, only(filter(n -> endswith(n, ".xml"), names("8"))))
        write(xml8, "edited")
        split_metronix_site(site)
        @test read(xml8, String) == "edited"
        @test_throws ErrorException split_metronix_site(site * ".2048")

        # a missing rate is made again, fully, after a split that stopped before
        # the end
        rm(site * ".2048"; recursive = true)
        mkpath(site * ".2048.partial/junk")
        @test !metronix_site_is_split(site)
        split_metronix_site(site)
        @test metronix_site_is_split(site)
        @test !ispath(site * ".2048.partial")
        @test _reproduces(site, [site * ".2048"]) == false   # only its rate
        @test sort(collect(keys(_tree(site * ".2048")))) ==
              sort([k for k in keys(_tree(site)) if occursin("2048H", k) || (endswith(k, ".kml") && startswith(k, b))])
    end
end

@testset "Metronix write with no masks changes nothing" begin
    mktempdir() do root
        site = joinpath(root, "DF999")
        _write_sample_mixed_site(site)
        before = _snapshot(site)
        @test write_metronix_site_masked(site) == String[]
        @test _snapshot(site) == before
        @test sort(readdir(root)) == ["DF999", "DF999.2048", "DF999.64", "DF999.8"]
        @test _reproduces(site, site .* [".8", ".64", ".2048"])

        @test_throws ErrorException write_metronix_site_masked(site; format = :MTH5)
        @test_throws ErrorException write_metronix_site_masked(site; format = :flat)
    end
end

@testset "Metronix write renumbers the runs after a cut" begin
    mktempdir() do root
        site = joinpath(root, "DF999")
        f = _write_sample_mixed_site(site)
        _name_xml_in_headers!(f.b)
        before = _snapshot(site)

        iv = (DateTime(2025, 4, 1, 7, 1, 5), DateTime(2025, 4, 1, 7, 1, 6))
        dirs = write_metronix_site_masked(site; rate_intervals = Dict(2048.0 => [iv]))
        @test dirs == [Timekeepers._norm_path(site * ".2048")]
        @test _snapshot(site) == before                   # the site never changes
        @test sort(readdir(root)) == ["DF999", "DF999.2048", "DF999.64", "DF999.8"]
        @test _reproduces(site, site .* [".8", ".64"]) == false   # 2048 Hz is not in them
        for r in ("8", "64")
            @test all(read(p) == read(joinpath(site, k)) for (k, p) in _tree(site * "." * r))
        end

        b = joinpath(site * ".2048", basename(f.b))
        names = readdir(b)
        # no change: burst R001, the scheduled XML of R002 and the .kml
        for n in names
            (occursin("_R001_", n) || occursin("_R002_", n) || endswith(n, ".kml")) || continue
            @test read(joinpath(f.b, n)) == read(joinpath(b, n))
        end
        # R000 keeps its number up to the cut. The remaining part is R003, after
        # R001 and R002
        @test count(n -> occursin(r"_R000_.*2048H\.ats$", n), names) == 5
        @test count(n -> occursin(r"_R003_.*2048H\.ats$", n), names) == 5
        xml0 = only(filter(n -> occursin(r"_R000_2048H\.xml$", n), names))
        xml3 = only(filter(n -> occursin(r"_R003_2048H\.xml$", n), names))
        @test xml0 == "999_2025-04-01_07-01-00_2025-04-01_07-01-05_R000_2048H.xml"
        @test xml3 == "999_2025-04-01_07-01-07_2025-04-01_07-01-20_R003_2048H.xml"
        text3 = read(joinpath(b, xml3), String)
        @test count("_R003_", text3) == 5 && !occursin("_R000_", text3)
        for n in filter(n -> occursin(r"_R00[03]_.*2048H\.ats$", n), names)
            @test Timekeepers._ats_xml_name(read(joinpath(b, n))) == (occursin("_R003_", n) ? xml3 : xml0)
        end
        # each run has its own XML, and the meas_ directory has its .kml
        @test all(r -> r.xml !== nothing, first(Timekeepers._metronix_site_index(site * ".2048")))
        @test "Site_meas_2025-04-01_07-00-05.kml" in names

        @test length(metronix_site_runs(site * ".2048")[2048.0]) == 3
        r3 = read_metronix(joinpath(b, xml3))
        @test start_time(r3) == DateTime(2025, 4, 1, 7, 1, 7)
        @test maximum(abs.(r3.channels[:e1].data .- f.burst0["Ex"][(7 * 2048 + 1):end])) < 1e-5
        # 07:01:05 to 07:01:06 is masked. 07:01:06 to 07:01:07 is also removed.
        # Thus, R003 starts on a whole second
        @test readlines(joinpath(site * ".2048", "mask.csv")) ==
              ["start_sample,end_sample,start_time,end_time",
               "$(5 * 2048 + 1),$(7 * 2048),2025-04-01T07:01:05.000,2025-04-01T07:01:06.999"]

        # the same masks cut nothing more the second time
        @test write_metronix_site_masked(site; rate_intervals = Dict(2048.0 => [iv])) == String[]

        # a meas_ directory with all its runs masked is removed
        whole = (DateTime(2025, 4, 1, 5, 0, 0), DateTime(2025, 4, 1, 6, 30, 0))
        dirs = @test_logs (:warn, r"was masked; removing it") match_mode = :any write_metronix_site_masked(site; intervals = [whole])
        @test dirs == [Timekeepers._norm_path(site * ".64")]
        @test readdir(site * ".64") == ["README.md", "mask.csv"]
        @test readlines(joinpath(site * ".64", "mask.csv"))[2] == "1,640,2025-04-01T06:00:00.000,2025-04-01T06:00:09.984"
    end
end

@testset "Metronix write cuts each rate by its own masks" begin
    mktempdir() do root
        site = joinpath(root, "DF999")
        f = _write_sample_mixed_site(site)
        split_metronix_site(site)
        # masked at 2048 Hz. The 8 Hz run records through the same second
        iv = (DateTime(2025, 4, 1, 7, 1, 0, 500), DateTime(2025, 4, 1, 7, 1, 0, 200))
        # the 0.2 s before the cut is too short to keep as a run. The cut ignores
        # it
        dirs = @test_logs (:warn, r"too short to store") match_mode = :any write_metronix_site_masked(site * ".2048"; rate_intervals = Dict(2048.0 => [iv]))
        @test dirs == [Timekeepers._norm_path(site * ".2048")]
        b = joinpath(site * ".2048", basename(f.b))

        # the 8 Hz run records through the cut, but it keeps all its samples
        @test all(read(p) == read(joinpath(site, k)) for (k, p) in _tree(site * ".8"))

        # burst R000 keeps only the part from the next whole second after the
        # cut, with its own number
        xml = only(filter(n -> occursin(r"_R000_2048H\.xml$", n), readdir(b)))
        @test xml == "999_2025-04-01_07-01-01_2025-04-01_07-01-20_R000_2048H.xml"
        seg = read_metronix(joinpath(b, xml))
        @test sampling_rate(seg) == 2048.0
        @test maximum(abs.(seg.channels[:e1].data .- f.burst0["Ex"][2049:end])) < 1e-5
        @test length(metronix_site_runs(site * ".2048")[2048.0]) == 2

        # the 0.2 s before the cut is also removed, with the remaining part of
        # that second
        @test readlines(joinpath(site * ".2048", "mask.csv"))[2:end] == ["1,2048,2025-04-01T07:01:00.000,2025-04-01T07:01:00.999"]
        @test occursin("2048 Hz runs of `DF999`", read(joinpath(site * ".2048", "README.md"), String))
    end
end

@testset "plot decimation breaks the line at a gap" begin
    secs = vcat(collect(0.0:0.001:9.999), collect(80_000.0:0.001:80_009.999))
    col = sin.(1:length(secs))
    masked = falses(length(secs))
    xs, yc, ym = Float64[], Float32[], Float32[]
    Timekeepers._decimate_minmax!(xs, yc, ym, secs, col, masked, 1, length(secs), 100; gap_s = 1.0)
    breaks = findall(i -> isnan(yc[i]) && isnan(ym[i]), eachindex(yc))
    @test length(breaks) == 1
    @test xs[breaks[1]] == 9.999
    Timekeepers._decimate_minmax!(xs, yc, ym, secs, col, masked, 1, length(secs), 100)
    @test !any(isnan, yc)                                   # no gap_s: no break
end

@testset "Metronix write lists what it cut in mask.csv" begin
    mktempdir() do root
        site = joinpath(root, "RK137")
        meas = joinpath(site, "meas_2025-04-01_07-00-05")
        _write_sample_metronix(meas; fs = 8, n = 80)

        # 07:00:08 to 07:00:10 is masked. The part after it starts on the next
        # whole second
        iv = (DateTime(2025, 4, 1, 7, 0, 8), DateTime(2025, 4, 1, 7, 0, 10))
        dirs = write_metronix_site_masked(site; intervals = [iv])
        @test dirs == [Timekeepers._norm_path(site * ".8")]
        csv = joinpath(only(dirs), "mask.csv")
        @test readlines(csv) == ["start_sample,end_sample,start_time,end_time",
                                 "17,40,2025-04-01T07:00:08.000,2025-04-01T07:00:10.875"]
        readme = read(joinpath(only(dirs), "README.md"), String)
        @test length(split(strip(readme), '\n')) <= 4
        @test occursin("`RK137`", readme)

        # a subsequent cut of the part after it counts samples from the run as
        # recorded
        iv2 = (DateTime(2025, 4, 1, 7, 0, 13), DateTime(2025, 4, 1, 7, 0, 14))
        write_metronix_site_masked(site; intervals = [iv])          # cuts nothing more
        write_metronix_site_masked(site; intervals = [iv2])
        @test readlines(csv)[2:end] == ["17,40,2025-04-01T07:00:08.000,2025-04-01T07:00:10.875",
                                        "57,72,2025-04-01T07:00:13.000,2025-04-01T07:00:14.875"]
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
        # Three runs of 4 s each. A 6 s gap between files 1 and 2. No gap
        # between files 2 and 3
        _write_sample_lemi424(joinpath(dir, "run_a.txt"); n = 4, start = DateTime(2020, 1, 1, 0, 0, 0))
        _write_sample_lemi424(joinpath(dir, "run_b.txt"); n = 4, start = DateTime(2020, 1, 1, 0, 0, 10))
        _write_sample_lemi424(joinpath(dir, "run_c.txt"); n = 4, start = DateTime(2020, 1, 1, 0, 0, 14))
        # The load must ignore a file that is not data
        write(joinpath(dir, "notes.md"), "ignore me")

        ta, fmt = Timekeepers._load_site_directory(dir)
        @test fmt == :lemi424
        times = Timekeepers._ta_timestamps(ta)
        vals = Timekeepers._ta_values(ta)
        @test first(times) == DateTime(2020, 1, 1, 0, 0, 0)
        @test last(times) == DateTime(2020, 1, 1, 0, 0, 17)
        @test size(vals, 1) == 18
        @test Timekeepers._ta_colnames(ta) == [:bx, :by, :bz, :e1, :e2]
        # Samples 5..10 are in the gap (indices 5..10 -> seconds 4..9, included)
        @test all(isnan, vals[5:10, :])
        # The first and the last data points must be finite
        @test all(isfinite, vals[1, :])
        @test all(isfinite, vals[end, :])
        meta = Timekeepers._ta_meta(ta)
        @test meta[:n_files] == 3
        @test meta[:sample_rate] == 1.0
        @test meta[:site] == basename(dir)

        # Write the joined SITENAME.txt automatically, and make sure that a new
        # scan ignores it
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

        # a record with only electric channels and without metadata: the writer
        # writes 0 for the missing channels
        times = [DateTime(2020, 1, 1) + Second(i - 1) for i in 1:6]
        e_only = TimeArray(times, [Float64(i) for i in 1:6, _ in 1:2], [:e1, :e2])
        write_geomag(out_path, e_only)
        back = load_geomag(out_path)
        @test values(back[:e1]) ≈ 1:6
        @test all(iszero, values(back[:bz]))
    end
end

@testset "recorded channels" begin
    @test Timekeepers._channel_recorded([1.0, 2.0, NaN])
    @test !Timekeepers._channel_recorded([NaN, NaN])
    @test !Timekeepers._channel_recorded(zeros(5))           # input that is not connected
    @test Timekeepers._channel_recorded([3.0])               # one sample, finite

    # a channel of only NaN is not recorded. It is not bad data, and it masks
    # nothing
    vals = [1.0 NaN; 2.0 NaN; NaN NaN; 4.0 NaN]
    mask = TimekeeperMask(TimeArray([DateTime(2020, 1, 1) + Second(i) for i in 1:4], vals, [:e1, :bz]))
    Timekeepers._auto_mask_nan!(mask, vals)
    @test mask.masked == [false, false, true, false]
end

@testset "decade ticks" begin
    ticks = Timekeepers._decade_ticks(3)
    @test first(ticks(0.3, 3.0e4)) == [1.0, 100.0, 1.0e4]    # each second decade after 3
    @test first(ticks(2.0, 700.0)) == [10.0, 100.0]
    @test first(ticks(0.2, 0.9)) == [0.2, 0.5]              # in one decade: 1-2-5
    @test all(l -> l isa AbstractString, last(ticks(3.1, 3.3)))  # plain numbers, not 10^0.5
end

@testset "loading window log wrapping" begin
    @test Timekeepers._wrap_words("short line", 20) == ["short line"]
    @test Timekeepers._wrap_words("one two three four", 9) == ["one two", "three", "four"]
    @test Timekeepers._wrap_words("a" ^ 12, 5) == ["aaaaa", "aaaaa", "aa"]
end

@testset "spectral workspace" begin
    n = 512
    x = [sin(2π * (i - 1) / 32) for i in 1:n]
    y = [cos(2π * (i - 1) / 32) for i in 1:n]
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
end

@testset "startup status" begin
    ta = _small_timearray()
    app = TKApp(ta; size = (700, 420))
    @test occursin("Timekeepers ready", app.status_label.text[])
end

function _wait_until(cond; timeout = 60)
    t0 = time()
    while !cond()
        time() - t0 > timeout && error("timed out waiting")
        sleep(0.05)
    end
end

@testset "typed plot window" begin
    @test Timekeepers._parse_window_count("12") == 12
    @test Timekeepers._parse_window_count("0") === nothing
    @test Timekeepers._parse_window_count("1.5") === nothing
    @test Timekeepers._parse_window_count("") === nothing

    n = 7200
    times = [DateTime(2020, 1, 1) + Second(i - 1) for i in 1:n]
    app = TKApp(TimeArray(times, randn(n, 2), [:bx, :by]); size = (900, 600))
    @test app.window_seconds[] == Inf                        # the full record
    app.window_box.stored_string[] = "30"
    @test app.window_seconds[] == Inf                        # All ignores the count
    app.window_menu.i_selected[] = 3                         # hours
    @test app.window_seconds[] == 30 * 3600.0
    app.window_menu.i_selected[] = 2                         # minutes
    @test app.window_seconds[] == 30 * 60.0
    app.window_menu.i_selected[] = 5                         # All
    @test app.window_seconds[] == Inf
end

@testset "spectra split at a gap between runs" begin
    mktempdir() do root
        site = joinpath(root, "DF999")
        _write_sample_mixed_site(site)
        ta = load_metronix_site(site; rate = 2048)            # two bursts, 2 min apart
        app = TKApp(ta; size = (900, 600))
        Timekeepers._switch_view!(app, :time_spectra)
        _wait_until(() -> length(app.psd_axes) == 5)

        # the gap ends a stretch, as a mask does
        segs, n_runs = Timekeepers._visible_good_index_segments(app, Timekeepers._visible_x_window(app)...)
        @test n_runs == 2
        @test segs == [(1, _BURST_N), (_BURST_N + 1, 2 * _BURST_N)]

        Timekeepers._compute_psd_for_window!(app)             # All: the average of the two runs
        @test !isempty(app.psd_values[1][])
        @test occursin("in 2 runs", app.psd_header[])
        nfft, _ = Timekeepers._current_nfft(app)
        @test occursin("averaged over $(Timekeepers._welch_segment_count([_BURST_N, _BURST_N], nfft)) segments",
                       app.psd_header[])

        # the samples of one run give the same spectrum as the average of the
        # two runs, if the two runs hold the same signal
        segs_one = [view(app.raw_values, 1:_BURST_N, 1)]
        f1, p1, _ = Timekeepers._welch_psd_segments(segs_one, 2048.0; nfft = nfft)
        @test app.psd_values[1][] ≈ p1[2:end]

    end
    @test Timekeepers._welch_segment_count([100, 4096, 1000], 1024) == 7
end

@testset "channel switches" begin
    n = 600
    times = [DateTime(2020, 1, 1) + Second(i - 1) for i in 1:n]
    vals = cumsum(sin.((1:n) ./ 7) .+ (1:5)'; dims = 1)
    vals[:, 3] .= NaN                                        # bz not recorded
    vals[:, 5] .= 0.0                                        # e2 not connected
    app = TKApp(TimeArray(times, vals, [:bx, :by, :bz, :e1, :e2]); size = (900, 600))
    @test app.channel_on == [true, true, false, true, false]
    @test masked_samples(app.mask) == 0

    app.channel_boxes[1].checked[] = false
    @test !app.channel_on[1]
    @test !app.axes[1].yticklabelsvisible[]

    # the traces stay. The spectra come after the new layout. The channels
    # that are off stay empty
    Timekeepers._switch_view!(app, :time_spectra)
    @test length(app.axes) == 5
    _wait_until(() -> length(app.psd_axes) == 5)
    @test app.channel_on == [false, true, false, true, false]
    @test isempty(app.psd_values[1][]) && !isempty(app.psd_values[2][])

    app.channel_boxes[1].checked[] = true
    @test !isempty(app.psd_values[1][])
end

@testset "Metronix mixed-rate site" begin
    mktempdir() do root
        site = joinpath(root, "RK999")
        _write_sample_metronix(joinpath(site, "meas_2025-04-01_07-00-05"); n = 400, fs = 8)
        _write_sample_metronix(joinpath(site, "meas_2025-04-01_09-00-00"); n = 800, fs = 16,
                               start_dt = DateTime(2025, 4, 1, 9, 0, 0))
        @test metronix_site_rates(site) == [8.0, 16.0]
        @test_throws ErrorException load_metronix_site(site)          # which rate?
        @test Timekeepers._sample_rate_from_timearray(load_metronix_site(site; rate = 16)) == 16.0

        ta8, fmt = Timekeepers._load_site_any(site; rate = 8.0)
        app = TKApp(ta8; size = (900, 600))
        Timekeepers._apply_loaded_data!(app, ta8, fmt, site; site_rates = [8.0, 16.0], rate = 8.0)
        @test length(app.rate_menu.options[]) == 2 && app.rate_index == 1
        t = timestamp(app.data)
        mask_interval!(app.mask, t[50], t[80])
        intervals_8 = copy(app.mask.intervals)

        # one rate in memory: when you select 16 Hz, it replaces the 8 Hz record,
        # and the 8 Hz record keeps its mask
        Timekeepers._select_rate!(app, 2)
        _wait_until(() -> app.rate_index == 2)
        @test Timekeepers._sample_rate_from_timearray(app.data) == 16.0
        @test app.rate_intervals[8.0] == intervals_8

        Timekeepers._select_rate!(app, 1)
        _wait_until(() -> app.rate_index == 1)
        @test app.mask.intervals == intervals_8
    end
end

# Write a position into each .ats header in `dir`, as an ADU does
function _set_ats_position!(dir::AbstractString, lat::Real, lon::Real)
    for (d, _, files) in walkdir(dir), f in files
        endswith(f, ".ats") || continue
        open(joinpath(d, f), "r+") do io
            seek(io, 96)
            write(io, Int32(round(lat * 3.6e6)), Int32(round(lon * 3.6e6)), Int32(15000))
        end
    end
end

@testset "survey scan, base and remote sites" begin
    mktempdir() do root
        t0 = DateTime(2025, 4, 1, 7, 0, 0)
        two_hours = 8 * 7200
        for (name, start, fs, lat, lon) in (("siteA", t0, 8, 48.60, 7.60),
                                            ("siteB", t0 + Hour(1), 8, 48.61, 7.62),
                                            ("siteC", t0 + Minute(30), 8, 48.90, 7.90),
                                            ("siteD", t0, 16, 48.605, 7.61))
            meas = joinpath(root, "campaign", name, "meas_" * Dates.format(start, "yyyy-mm-dd_HH-MM-SS"))
            _write_sample_metronix(meas; n = fs == 8 ? two_hours : 16 * 7200, fs = fs, start_dt = start)
            _set_ats_position!(meas, lat, lon)
        end
        # a split rate directory next to its site contains the same data again
        cp(joinpath(root, "campaign", "siteA"), joinpath(root, "campaign", "siteA.8"))
        lemi = mkpath(joinpath(root, "campaign", "siteE"))
        _write_sample_lemi424(joinpath(lemi, "siteE_001.txt"); n = 600, start = DateTime(2020, 1, 1))
        # a telluric site next to siteA, which records with it: only electric
        # channels
        tel = joinpath(root, "campaign", "siteF", "meas_2025-04-01_07-00-00")
        _write_sample_metronix(tel; n = two_hours, fs = 8, start_dt = t0)
        foreach(f -> occursin("_TH", f) && rm(joinpath(tel, f)), readdir(tel))
        _set_ats_position!(tel, 48.601, 7.601)
        _write_sample_geomag(joinpath(mkpath(joinpath(root, "campaign", "siteG")), "siteG.txt"))

        s = scan_survey(root)
        @test [x.name for x in s] == ["siteA", "siteB", "siteC", "siteD", "siteE", "siteF", "siteG"]
        @test length(scan_survey(root; include_split = true)) == 8
        a, b, c, dsite, e, f, g = s.sites
        @test site_components(a) == [:e1, :e2, :bx, :by, :bz] && has_magnetic(a)
        @test site_components(f) == [:e1, :e2] && !has_magnetic(f)
        @test site_components(e) == [:e1, :e2, :bx, :by, :bz]                # LEMI-424
        @test g.format == :geomag && site_components(g) == [:e1, :e2, :bx, :by, :bz]
        @test isapprox(g.latitude, 60 + 35 / 60 + 14.4 / 3600; atol = 1e-6)
        @test only(g.runs).sample_rate ≈ 10
        # an input that is not connected gives zeros at the two ends: not
        # recorded
        @test Timekeepers._text_components(["2020 01 01 00 00 00 1 2 3 10 11 0 0 0 0",
                                            "2020 01 01 00 00 01 1 2 3 10 11 0.000 0.0 0 0"], :lemi424) ==
              [:bx, :by, :bz]
        @test isempty(site_references(s, "siteF"; min_overlap_hours = 0.5).base) == false   # F uses the field of A
        @test all(x -> x.site != "siteF", site_references(s, "siteA"; min_overlap_hours = 0.5).base)
        @test a.format == :metronix && e.format == :lemi424
        @test a.latitude ≈ 48.60 atol = 1e-6
        @test e.latitude ≈ 60 + 22 / 60 atol = 1e-6
        @test e.longitude ≈ 24 + 56 / 60 atol = 1e-6
        @test only(a.runs).stop - only(a.runs).start == Hour(2)
        @test only(e.runs).sample_rate == 1.0 && only(e.runs).n_samples == 600
        @test survey_rates(s) == [1.0, 8.0, 10.0, 16.0]

        @test overlap_seconds(a, b) ≈ 3600
        @test overlap_seconds(a, c; rate = 8) ≈ 5400
        @test overlap_seconds(a, dsite) == 0                       # same time, other rate
        @test overlap_seconds(a, dsite; rate = :all) ≈ 7200          # ... but :all makes a pair
        @test overlap_seconds(a, b; rate = :all) ≈ 3600
        @test common_window(s, "siteA", ["siteB", "siteD"]; rate = :all) == [(t0 + Hour(1), t0 + Hour(2))]
        @test "siteD" in [x.site for x in site_references(s, "siteA"; rate = :all).base]
        @test overlap_seconds(a, e) == 0
        @test overlap_intervals(a, b) == [(t0 + Hour(1), t0 + Hour(2))]
        m = overlap_matrix(s; rate = 8.0)
        @test m[1, 1] ≈ 2 && m[1, 2] ≈ 1 && m[2, 1] ≈ 1 && m[1, 4] == 0
        @test 1.5 < site_distance(a, b) < 2.2
        @test site_distance(a, c) > 30

        refs = site_references(s, "siteA"; base_km = 5, min_overlap_hours = 0.5)
        @test [x.site for x in refs.base] == ["siteB"]
        @test [x.site for x in refs.remote] == ["siteC"]
        @test refs.remote[1].overlap_fraction ≈ 0.75
        @test isempty(site_references(s, "siteA"; min_overlap_hours = 1.5).base)
        @test isempty(site_references(s, "siteA"; base_km = 100, min_overlap_hours = 0.5).remote)
        gap = site_references(s, "siteA"; base_km = 5, remote_km = 50, min_overlap_hours = 0.5)
        @test [x.site for x in gap.base] == ["siteB"] && isempty(gap.remote)   # siteC is between
        @test site_references(s, a; exclude = ["siteB"], min_overlap_hours = 0.5).base[1].excluded

        plan = reference_plan(s; base_km = 5, min_overlap_hours = 0.5)
        @test plan[1].base == ["siteB"] && plan[1].remote == ["siteC"]
        @test isempty(plan[5].base) && isempty(plan[5].remote)
        plan2 = reference_plan(s; base_km = 5, min_overlap_hours = 0.5,
                               exclude = Dict("siteA" => Set(["siteC"])))
        @test plan2[1].base == ["siteB"] && isempty(plan2[1].remote)

        # siteA records 07:00-09:00, siteB 08:00-10:00, siteC 07:30-09:30
        @test common_window(s, "siteA", ["siteB", "siteC"]) == [(t0 + Hour(1), t0 + Hour(2))]
        @test common_window(s, "siteA", ["siteC"]) == [(t0 + Minute(30), t0 + Hour(2))]
        @test isempty(common_window(s, "siteA", String[]))
        @test isempty(common_window(s, "siteA", ["siteD"]))                  # other rate
        @test plan[1].common == [(t0 + Hour(1), t0 + Hour(2))]

        @test plan[1].base_hours ≈ [1.0] && plan[1].remote_hours ≈ [1.5]

        path = write_reference_plan(joinpath(root, "plan.txt"), s; rate = 8.0, min_overlap_hours = 0.5)
        lines = readlines(path)
        @test length(lines) == 8                                        # one header line, one row for each site
        @test split(lines[1]) == ["site", "base", "overlap", "(h)", "remote", "overlap", "(h)"]
        @test split(lines[2]) == ["siteA", "siteB", "1.00", "siteC", "1.50"]
        @test split(lines[6]) == ["siteE", "-", "-", "-", "-"]
        @test first(findfirst("siteC", lines[2])) == first(findfirst("remote", lines[1]))   # aligned
        back = read_reference_plan(path)
        @test length(back) == 7
        @test back[2].site == "siteB" && back[2].base == ["siteA"] && back[2].remote == ["siteC"]
        @test back[2].base_hours ≈ [1.0] && back[2].remote_hours ≈ [1.5]
        @test isempty(back[5].base) && isempty(back[5].remote_hours)
        wide = write_reference_plan(joinpath(root, "wide.txt"),
            [(site = "Sarıçam", base = ["a", "b"], base_hours = [12.345, 1.0],
              remote = String[], remote_hours = Float64[]),
             (site = "x", base = String[], base_hours = Float64[],
              remote = ["Sarıçam"], remote_hours = [3.0])])
        r1, r2 = read_reference_plan(wide)
        @test r1.site == "Sarıçam" && r1.base == ["a", "b"] && r1.base_hours ≈ [12.35, 1.0]
        @test isempty(r1.remote) && r2.remote == ["Sarıçam"] && r2.remote_hours ≈ [3.0]

        @test TKDash(s; size = (1200, 800)).rate === :all            # different instruments: all rates
        dash = TKDash(s; size = (1200, 800), min_overlap_hours = 0.5, rate = 8.0)
        @test dash.rate == 8.0 && dash.focus == 0                 # opens on the overview
        @test dash.status.text[] == Timekeepers.DASH_HINT          # the instruction, in grey
        @test length(dash.charts) == 1 && dash.chart_rows == [1:7]
        Timekeepers._focus!(dash, 3)                               # siteC: base chart above remote chart
        @test dash.site_menu.i_selected[] == 4
        @test length(dash.charts) == 2
        @test dash.chart_rows == [[3], [3, 2, 1]]                  # remote tie: nearer first
        @test dash.common == [(t0 + Hour(1), t0 + Hour(2))]              # siteC with siteA and siteB
        # hover: the row of siteB in the remote chart of siteC, on its run
        x = Dates.value(t0 + Hour(2) - dash.t0) / 3.6e6                     # 09:00 on the chart
        h = Timekeepers._hover(dash, dash.charts[2], Timekeepers.Point2f(x, 2.0))
        @test h !== nothing && occursin("siteB · remote", h[2]) && occursin("run ", h[2])
        @test Timekeepers._hover(dash, dash.charts[2], Timekeepers.Point2f(x, 2.5)) === nothing
        @test startswith(Timekeepers._site_tip(dash, 1), "siteA · remote\n")
        Timekeepers._focus!(dash, 0)
        @test length(dash.charts) == 1 && isempty(dash.common)
        @test startswith(Timekeepers._site_tip(dash, 1), "siteA\nown 2.0h")
        @test Timekeepers._bar_label(42732) == "[11.9h]"
        @test Timekeepers._bar_label(5400, 34200) == "[1.5h/9.5h]"
        Timekeepers._focus!(dash, 1)
        @test dash.chart_rows[1] == [1, 2]
        Timekeepers._toggle!(dash, 2)                              # drop siteB from siteA
        @test reference_plan(dash)[1].base == String[]
        @test dash.chart_rows[1] == [1]                            # ... and from the charts
        @test Timekeepers._site_roles(dash)[2][2] == Timekeepers.DASH_IDLE   # hollow on the map
        Timekeepers._toggle!(dash, 2)
        @test reference_plan(dash)[1].base == ["siteB"]
        Timekeepers._toggle!(dash, 2)
        Timekeepers._toggle!(dash, 3)                              # drop the two, then restore them
        @test length(dash.exclude["siteA"]) == 2
        Timekeepers._restore!(dash)
        @test !haskey(dash.exclude, "siteA") && dash.chart_rows == [[1, 2], [1, 3]]
        @test occursin("Restored 2 dropped sites for siteA", dash.status.text[])
        Timekeepers._toggle!(dash, 5)                              # siteE is neither: no change
        @test !haskey(dash.exclude, "siteE") && reference_plan(dash)[1].remote == ["siteC"]
        Timekeepers.limits!(dash.map_axis, 7.5, 7.7, 48.55, 48.65)            # zoom in
        zoomed = dash.map_axis.targetlimits[]
        Timekeepers._focus!(dash, 2)                                         # a new site keeps the zoom
        @test dash.map_axis.targetlimits[] == zoomed
        Timekeepers._reset_zoom!(dash)
        @test dash.map_axis.targetlimits[].origin[1] < 7.6 - 0.01
        @test dash.zoom_button !== nothing                         # below the map
        dash.map_open = false
        Timekeepers._layout_map!(dash)
        @test dash.map_axis === nothing && dash.zoom_button === nothing
        @test isempty(TKDash(Timekeepers.Survey(root, Timekeepers.SurveySite[])).chart_rows)
    end
end

@testset "survey site with a gap between runs" begin
    mktempdir() do root
        t0 = DateTime(2025, 4, 1, 7, 0, 0)
        # siteX records 07:00-08:00 and 09:00-10:00 in two meas_ directories.
        # siteY records 07:30-09:30 without a stop
        for (name, start, lat) in (("siteX", t0, 48.60), ("siteX", t0 + Hour(2), 48.60),
                                   ("siteY", t0 + Minute(30), 48.61))
            meas = joinpath(root, name, "meas_" * Dates.format(start, "yyyy-mm-dd_HH-MM-SS"))
            _write_sample_metronix(meas; n = name == "siteX" ? 8 * 3600 : 8 * 7200, fs = 8, start_dt = start)
            _set_ats_position!(meas, lat, 7.60)
        end
        s = scan_survey(root)
        x, y = s.sites
        @test length(x.runs) == 2 && recording_seconds(x) ≈ 7200
        @test recording_intervals(x) == [(t0, t0 + Hour(1)), (t0 + Hour(2), t0 + Hour(3))]
        @test overlap_intervals(x, y) == [(t0 + Minute(30), t0 + Hour(1)), (t0 + Hour(2), t0 + Hour(2) + Minute(30))]
        @test overlap_seconds(x, y) ≈ 3600
        refs = site_references(s, "siteX")
        @test only(refs.base).site == "siteY" && only(refs.base).overlap_fraction ≈ 0.5
        @test length(common_window(s, "siteX", ["siteY"])) == 2

        dash = TKDash(s; size = (1200, 800))
        Timekeepers._focus!(dash, 1)
        @test dash.chart_rows[1] == [1, 2] && length(dash.common) == 2
        # a hover on siteX in its gap gives the site but no run
        gap = Dates.value(t0 + Minute(90) - dash.t0) / 3.6e6
        h = Timekeepers._hover(dash, dash.charts[1], Timekeepers.Point2f(gap, 1.0))
        @test h !== nothing && startswith(h[2], "siteX") && !occursin("run ", h[2])
        inrun = Dates.value(t0 + Minute(150) - dash.t0) / 3.6e6
        @test occursin("run ", Timekeepers._hover(dash, dash.charts[1], Timekeepers.Point2f(inrun, 1.0))[2])
        plan = only(r for r in reference_plan(dash) if r.site == "siteX")
        @test plan.base == ["siteY"] && plan.base_hours ≈ [1.0]
    end
end
