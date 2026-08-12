using Random
using Statistics: median

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
    time_zero = _fit_time_zero(lightcurve, prior_kwargs)
    fit_times = inputs.times .- time_zero
    fit_prior_kwargs = _centered_prior_kwargs(prior_kwargs, time_zero)
    fit_priors = priors === nothing ? default_dust_priors(
        fit_times,
        inputs.wavelengths,
        inputs.observed;
        fit_prior_kwargs...,
    ) : priors
    chain = sample_dust_transit(
        grid,
        fit_times,
        inputs.wavelengths,
        inputs.observed,
        inputs.errors,
        star,
        fit_priors;
        kwargs...,
    )
    samples = fit_samples_table(chain; scale=:physical, t0_offset=time_zero)
    sampled_samples = fit_samples_table(chain; scale=:sampled, t0_offset=time_zero)
    predictive_fluxes = n_predictive_draws === nothing ? nothing : posterior_predictive_fluxes(
        chain,
        grid,
        fit_times,
        inputs.wavelengths,
        star;
        n_draws=n_predictive_draws,
        rng=rng,
    )
    predictive_intervals = predictive_fluxes === nothing ? nothing : posterior_predictive_intervals(predictive_fluxes)

    (;
        lightcurve,
        inputs,
        fit_inputs=(;
            times=fit_times,
            wavelengths=inputs.wavelengths,
            observed=inputs.observed,
            errors=inputs.errors,
        ),
        time_zero,
        fit_time_mode=:relative,
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

function _fit_time_zero(lightcurve::LightCurve, prior_kwargs)
    t0_center = _prior_kwarg(prior_kwargs, :t0_center, nothing)
    t0_center === nothing || return Float64(t0_center)

    median(Float64.(lightcurve.time))
end

function _centered_prior_kwargs(prior_kwargs, time_zero)
    prior_kwargs === nothing && return (;)

    values = Dict{Symbol,Any}()
    for (key, value) in pairs(prior_kwargs)
        symbol = Symbol(key)
        values[symbol] = symbol === :t0_center && value !== nothing ? Float64(value) - time_zero : value
    end

    (; values...)
end

function _prior_kwarg(prior_kwargs, key::Symbol, default)
    prior_kwargs === nothing && return default

    try
        if haskey(prior_kwargs, key)
            return prior_kwargs[key]
        end
    catch
    end

    string_key = string(key)
    try
        if haskey(prior_kwargs, string_key)
            return prior_kwargs[string_key]
        end
    catch
    end

    if prior_kwargs isa NamedTuple && hasproperty(prior_kwargs, key)
        return getproperty(prior_kwargs, key)
    end

    default
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
