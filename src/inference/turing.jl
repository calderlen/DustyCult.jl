using ADTypes: AutoForwardDiff
using Distributions: Normal

import ForwardDiff
import Turing

function chi2(observed, modeled, errors)
    total = 0.0

    for (obs, mod, err) in zip(observed, modeled, errors)
        total += ((obs - mod) / err)^2
    end

    total
end

function reduced_chi2(observed, modeled, errors, n_parameters)
    dof = length(observed) - n_parameters

    chi2(observed, modeled, errors) / dof
end

function dust_model_fluxes(grid, times, wavelengths, star::StarParams, dust::DustOccultationParams)
    [
        dust_transit_flux(grid, t, lambda, star, dust)
        for (t, lambda) in zip(times, wavelengths)
    ]
end

function dust_model_fluxes(grid, times, wavelength::Real, star::StarParams, dust::DustOccultationParams)
    [dust_transit_flux(grid, t, wavelength, star, dust) for t in times]
end

Turing.@model function dust_transit_model(grid, times, wavelengths, observed, errors, star, priors)
    t0 ~ priors.t0
    log_v ~ priors.log_v
    b ~ priors.b
    log_tau0 ~ priors.log_tau0
    log_lambda0 ~ priors.log_lambda0
    alpha ~ priors.alpha
    log_sigma_y ~ priors.log_sigma_y
    log_sigma_x_plus ~ priors.log_sigma_x_plus
    log_sigma_x_minus ~ priors.log_sigma_x_minus

    dust = DustOccultationParams(
        t0 = t0,
        v = exp(log_v),
        b = b,
        tau0 = exp(log_tau0),
        lambda0 = exp(log_lambda0),
        alpha = alpha,
        sigma_y = exp(log_sigma_y),
        sigma_x_plus = exp(log_sigma_x_plus),
        sigma_x_minus = exp(log_sigma_x_minus),
    )

    modeled = dust_model_fluxes(grid, times, wavelengths, star, dust)

    for i in eachindex(observed)
        observed[i] ~ Normal(modeled[i], errors[i])
    end
end

function sample_dust_transit(grid, times, wavelengths, observed, errors, star::StarParams, priors::DustPriors; n_samples = 1000, n_adapt = 1000, n_chains = 4, parallel = Turing.MCMCSerial(), progress = true, kwargs...)
    model = dust_transit_model(grid, times, wavelengths, observed, errors, star, priors)
    sampler = Turing.NUTS(; adtype=AutoForwardDiff())

    Turing.sample(
        model,
        sampler,
        parallel,
        n_samples,
        n_chains;
        num_warmup=n_adapt,
        discard_initial=n_adapt,
        progress=progress,
        kwargs...,
    )
end
