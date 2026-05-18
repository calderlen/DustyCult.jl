module DustyCult

export StarParams
export DustOccultationParams
export mu, limb_darkening, I, stellar_grid
export gaussian, gaussian2d, two_gaussian2d
export powerlaw_opacity, optical_depth, transmission
export occulter_center
export F0, dust_transit_flux
export lnL, chi2, reduced_chi2
export dust_model_fluxes
export DustPriors, default_dust_priors
export dust_transit_model, sample_dust_transit
export fit_dust_transit
export fit_relative_flux_lightcurve, fit_delta_magnitude_lightcurve
export fit_magnitude_lightcurve
export valid_quadratic_limb_darkening, valid_positive_scales
export Bandpass, LightCurve
export delta_mag_to_relative_flux, relative_flux_to_delta_mag
export mag_error_to_relative_flux_error
export lightcurve_from_relative_fluxes, read_relative_flux_lightcurve
export lightcurve_from_delta_magnitudes, read_delta_magnitude_lightcurve
export lightcurve_from_magnitudes, read_magnitude_lightcurve
export fit_inputs
export DEFAULT_FIT_LABELS
export fit_samples_table, posterior_predictive_fluxes
export posterior_predictive_intervals
export read_fit_config, run_fit_config, fit_lightcurve_from_config
export build_bandpass, build_star_params, build_stellar_grid
export write_fit_outputs, write_fit_manifest
export lightcurve_output_table, predictive_intervals_table, predictive_fluxes_table
export parse_fit_cli_args, fit_lightcurve_cli

include("star/parameters.jl")
include("star/intensity.jl")
include("star/grid.jl")

include("occulter/shape.jl")
include("occulter/opacity.jl")
include("occulter/trajectory.jl")

include("model.jl")
include("lightcurve.jl")
include("inference/likelihood.jl")
include("inference/priors.jl")
include("inference/turing.jl")
include("posterior.jl")
include("fit.jl")
include("bridge.jl")

end
