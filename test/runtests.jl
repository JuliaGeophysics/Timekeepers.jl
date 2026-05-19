using Dates
using Bonito
using Test
using Timekeepers
using TimeSeries

const LEMI_SAMPLE_TEXT = """
2020 10 04 00 00 00 23772.512   238.148 41845.187  33.88  25.87   142.134   -45.060   213.787     8.224 12.78 2199.0 3404.83963 N 10712.84474 W 12 2 0
2020 10 04 00 00 01 23772.514   238.159 41845.202  33.89  25.87   142.014   -45.172   213.662     8.112 12.78 2199.0 3404.83965 N 10712.84471 W 12 2 0
2020 10 04 00 00 02 23772.522   238.149 41845.191  33.90  25.87   142.033   -45.156   213.688     8.136 12.78 2199.1 3404.83964 N 10712.84469 W 12 2 0
2020 10 04 00 00 03 23772.537   238.143 41845.210  33.89  25.87   142.023   -45.122   213.679     8.177 12.78 2199.2 3404.83965 N 10712.84467 W 12 2 0
"""

function lemi_sample_file()
    path = tempname() * ".txt"
    write(path, LEMI_SAMPLE_TEXT)
    return path
end

@testset "LEMI-424 reader" begin
    LEMI_SAMPLE = lemi_sample_file()
    run = read_lemi424(LEMI_SAMPLE; site = "LEMI090")
    @test run.site == "LEMI090"
    @test run.instrument == "LEMI-424"
    @test run.source_format == :lemi424
    @test sampling_rate(run) == 1.0
    @test length(run.channels[:bx].data) == 4
    @test start_time(run) == DateTime(2020, 10, 4, 0, 0, 0)
    @test isapprox(run.channels[:bx].data[1], 23772.512; atol = 1e-8)
    @test isapprox(run.metadata[:latitude], 34.08066; atol = 1e-4)
    @test isapprox(run.metadata[:longitude], -107.21406; atol = 1e-4)
end

@testset "LEMI-424 TimeArray loader" begin
    LEMI_SAMPLE = lemi_sample_file()
    full = read_lemi424(LEMI_SAMPLE; site = "LEMI090")
    tv = load_lemi424(LEMI_SAMPLE; site = "LEMI090")
    @test tv isa TimeArray
    @test timestamp(tv)[1] == DateTime(2020, 10, 4, 0, 0, 0)
    @test colnames(tv) == [:bx, :by, :bz, :e1, :e2]
    @test values(tv)[:, 1] == full.channels[:bx].data
    @test values(tv)[:, 2] == full.channels[:by].data
    @test values(tv)[:, 3] == full.channels[:bz].data
    @test values(tv)[:, 4] == full.channels[:e1].data
    @test values(tv)[:, 5] == full.channels[:e2].data
end

@testset "TimeArray conversion" begin
    LEMI_SAMPLE = lemi_sample_file()
    tv = to_timearray(read_lemi424(LEMI_SAMPLE; site = "LEMI090"); components = [:bx, :by, :bz, :e1, :e2])
    @test tv isa TimeArray
    @test length(timestamp(tv)) == 4
    @test colnames(tv) == [:bx, :by, :bz, :e1, :e2]
    @test values(tv)[1, 1] ≈ 23772.512
    run = from_timearray(tv; site = "LEMI090", instrument = "LEMI-424")
    @test run.channels[:e1].data == values(tv)[:, 4]
end

@testset "All-component masking" begin
    LEMI_SAMPLE = lemi_sample_file()
    tv = load_lemi424(LEMI_SAMPLE; site = "LEMI090")
    mask = TimekeeperMask(tv)
    mask_interval!(mask, DateTime(2020, 10, 4, 0, 0, 1), DateTime(2020, 10, 4, 0, 0, 2))
    @test masked_samples(mask) == 2
    @test sample_weights(mask) == [1.0, 0.0, 0.0, 1.0]
    cleaned = cleaned_timearray(tv, mask)
    @test all(isnan, values(cleaned)[2, :])
    @test all(isnan, values(cleaned)[3, :])
    @test values(cleaned)[1, 1] == values(tv)[1, 1]
    dropped = cleaned_timearray(tv, mask; mode = :drop)
    @test length(timestamp(dropped)) == 2
    segments = good_segments(tv, mask)
    @test length(segments) == 2
    tmp = tempname() * ".csv"
    write_mask(tmp, mask)
    loaded = read_mask(tmp, tv)
    @test loaded.masked == mask.masked
    cleaned_path = tempname() * ".csv"
    write_cleaned(cleaned_path, tv, mask)
    @test occursin("NaN", read(cleaned_path, String))
    unmask_interval!(mask, DateTime(2020, 10, 4, 0, 0, 2), DateTime(2020, 10, 4, 0, 0, 2))
    @test masked_samples(mask) == 1
    clear_mask!(mask)
    @test masked_samples(mask) == 0
end

@testset "TKApp" begin
    LEMI_SAMPLE = lemi_sample_file()
    tv = load_lemi424(LEMI_SAMPLE; site = "LEMI090")
    app = TKApp(tv; window_samples = 2)
    empty_app = TKApp(; window_samples = 1)
    path_app = TKApp(LEMI_SAMPLE; site = "LEMI090", window_samples = 2)
    @test app isa TKApp
    @test empty_app isa TKApp
    @test app.app isa Bonito.App
    @test app.app.handler(nothing, nothing) !== nothing
    @test empty_app.app.handler(nothing, nothing) !== nothing
    @test path_app isa TKApp
end

@testset "LEMI-424 writer" begin
    LEMI_SAMPLE = lemi_sample_file()
    run = read_lemi424(LEMI_SAMPLE; site = "LEMI090")
    tmp = tempname() * ".txt"
    write_lemi424(tmp, run)
    reread = read_lemi424(tmp; site = "LEMI090-copy")
    @test length(reread.channels[:bx].data) == length(run.channels[:bx].data)
    @test reread.channels[:e2].data[3] ≈ run.channels[:e2].data[3]
end
