using Random

"""
Fit a `LightCurve` with the dust transit model and return posterior products.
"""
function fit_dust_transit(lightcurve::LightCurve, grid, star::StarParams;
    priors = nothing,
    prior_kwargs = (;),
    n_predictive_draws = 200,
    rng = Random.default_rng(),
    kwargs...,
)
    inputs = fit_inputs(lightcurve)
    fit_priors = priors === nothing ? default_dust_priors(lightcurve; prior_kwargs...) : priors
    chain = sample_dust_transit(
        grid,
        inputs.times,
        inputs.wavelengths,
        inputs.observed,
        inputs.errors,
        star,
        fit_priors;
        kwargs...,
    )
    samples = fit_samples_table(chain; scale=:physical)
    sampled_samples = fit_samples_table(chain; scale=:sampled)
    predictive_fluxes = n_predictive_draws === nothing ? nothing : posterior_predictive_fluxes(
        chain,
        grid,
        inputs.times,
        inputs.wavelengths,
        star;
        n_draws=n_predictive_draws,
        rng=rng,
    )
    predictive_intervals = predictive_fluxes === nothing ? nothing : posterior_predictive_intervals(predictive_fluxes)

    (;
        lightcurve,
        inputs,
        priors=fit_priors,
        chain,
        samples,
        sampled_samples,
        predictive_fluxes,
        predictive_intervals,
    )
end

"""
Read a relative-flux light curve, fit it, and return posterior products.
"""
function fit_relative_flux_lightcurve(path, bandpass, grid, star::StarParams;
    time_col = :time,
    relative_flux_col = :relative_flux,
    relative_flux_error_col = :relative_flux_error,
    band_col = :band,
    kwargs...,
)
    lightcurve = read_relative_flux_lightcurve(
        path,
        bandpass;
        time_col=time_col,
        relative_flux_col=relative_flux_col,
        relative_flux_error_col=relative_flux_error_col,
        band_col=band_col,
    )

    fit_dust_transit(lightcurve, grid, star; kwargs...)
end

"""
Read a delta-magnitude light curve, fit it, and return posterior products.
"""
function fit_delta_magnitude_lightcurve(path, bandpass, grid, star::StarParams;
    time_col = :time,
    delta_mag_col = :delta_mag,
    mag_error_col = :mag_error,
    band_col = :band,
    kwargs...,
)
    lightcurve = read_delta_magnitude_lightcurve(
        path,
        bandpass;
        time_col=time_col,
        delta_mag_col=delta_mag_col,
        mag_error_col=mag_error_col,
        band_col=band_col,
    )

    fit_dust_transit(lightcurve, grid, star; kwargs...)
end

"""
Read a magnitude light curve, fit it, and return posterior products.
"""
function fit_magnitude_lightcurve(path, bandpass, reference_magnitude, grid, star::StarParams;
    time_col = :time,
    mag_col = :mag,
    mag_error_col = :mag_error,
    band_col = :band,
    kwargs...,
)
    lightcurve = read_magnitude_lightcurve(
        path,
        bandpass,
        reference_magnitude;
        time_col=time_col,
        mag_col=mag_col,
        mag_error_col=mag_error_col,
        band_col=band_col,
    )

    fit_dust_transit(lightcurve, grid, star; kwargs...)
end
