using DataFrames
using JSON3
using Statistics: median
using Test

using DustyCult

@testset "relative-time fitting metadata" begin
    bandpass = Dict("g" => Bandpass(477.0), "V" => Bandpass(545.0))
    times = [2.4600000e6, 2.4600010e6, 2.4600020e6]
    lightcurve = lightcurve_from_relative_fluxes(
        times,
        [1.0, 0.92, 0.99],
        [0.02, 0.02, 0.02],
        ["g", "V", "g"],
        bandpass,
    )

    time_zero = DustyCult._fit_time_zero(lightcurve, (; t0_center=2.4600010e6))
    centered = DustyCult._centered_prior_kwargs((; t0_center=2.4600010e6, t0_width=7.0), time_zero)

    @test time_zero == 2.4600010e6
    @test centered.t0_center == 0.0
    @test centered.t0_width == 7.0
    @test DustyCult._fit_time_zero(lightcurve, (;)) == median(times)
end

@testset "absolute output artifacts" begin
    source = DataFrame(
        t0=[-0.25, 0.25],
        log_v=log.([1.0, 1.2]),
        b=[0.0, 0.1],
        log_tau0=log.([0.3, 0.4]),
        log_lambda0=log.([510.0, 520.0]),
        alpha=[0.0, 0.1],
        log_sigma_y=log.([0.2, 0.25]),
        log_sigma_x_plus=log.([0.3, 0.35]),
        log_sigma_x_minus=log.([0.4, 0.45]),
    )
    samples = fit_samples_table(source; scale=:physical, t0_offset=2.4600010e6)
    sampled = fit_samples_table(source; scale=:sampled, t0_offset=2.4600010e6)

    @test samples.t0 == [2.46000075e6, 2.46000125e6]
    @test sampled.t0 == [2.46000075e6, 2.46000125e6]
    @test samples.v ≈ [1.0, 1.2]

    bandpass = Dict("g" => Bandpass(477.0), "V" => Bandpass(545.0))
    lightcurve = lightcurve_from_relative_fluxes(
        [2.4600000e6, 2.4600010e6],
        [1.0, 0.95],
        [0.02, 0.02],
        ["g", "V"],
        bandpass,
    )
    intervals = (;
        lower95=[0.9, 0.91],
        lower68=[0.94, 0.945],
        median=[0.99, 0.95],
        upper68=[1.01, 0.98],
        upper95=[1.03, 1.0],
    )
    table = predictive_intervals_table((; lightcurve, predictive_intervals=intervals))

    @test table.time == lightcurve.time

    manifest_path = tempname()
    write_fit_manifest(
        manifest_path,
        Dict("prior_kwargs" => Dict("t0_center" => 2.4600010e6)),
        Dict("samples" => "samples.parquet");
        fit=(; time_zero=2.4600010e6, fit_time_mode=:relative),
    )
    manifest = JSON3.read(read(manifest_path, String))

    @test manifest.fit_time.mode == "relative"
    @test manifest.fit_time.time_zero == 2.4600010e6
end
