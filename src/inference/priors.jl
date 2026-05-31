using Distributions: Normal
using Statistics: median

Base.@kwdef struct DustPriors
    t0
    log_v
    b
    log_tau0
    log_lambda0
    alpha
    log_sigma_y
    log_sigma_x_plus
    log_sigma_x_minus
end

function default_dust_priors(lightcurve::LightCurve; kwargs...)
    default_dust_priors(lightcurve.time, lightcurve.wavelength, lightcurve.relative_flux; kwargs...)
end

"""
The priors are centered on the deepest observed point, the median wavelength, and a half-depth dip-duration estimate. Positive model scales are represented as log-space normal priors to match `dust_transit_model`.
"""
function default_dust_priors(times, wavelengths, observed;
    t0_center = nothing, # time when the dust occulter center is at x=0
    t0_width = nothing,
    log_v_width = 1.0,
    b_center = 0.0, # impact parameter
    b_width = 0.5,
    log_tau0_width = 1.5,
    log_lambda0_width = 0.25,
    alpha_center = 0.0, # dust opacity color slope, so alpha > 0 means more opaque at shorter wavelengths
    alpha_width = 2.0,
    log_sigma_width = 0.75, 
    min_depth = 1e-3, # optical depth floor
    min_duration = 1e-6, #time floor
)
    time = collect(times)
    flux = collect(observed)
    wavelength = wavelengths isa Real ? [wavelengths] : collect(wavelengths)

    length(time) == length(flux) || throw(ArgumentError("times and observed must have the same length"))
    isempty(time) && throw(ArgumentError("times and observed cannot be empty"))
    isempty(wavelength) && throw(ArgumentError("wavelengths cannot be empty"))

    time_span = maximum(time) - minimum(time)
    time_scale = max(time_span, min_duration)
    baseline = median(flux)
    dip_index = argmin(flux)
    depth = max(baseline - flux[dip_index], min_depth)
    duration = max(_dip_duration_guess(time, flux, baseline, depth), min_duration)

    t0_guess = t0_center === nothing ? time[dip_index] : Float64(t0_center)
    v_guess = 2.0 / duration
    tau0_guess = depth
    lambda0_guess = median(wavelength)
    sigma_guess = 0.25

    lambda0_guess > 0 || throw(ArgumentError("wavelengths must be positive"))

    DustPriors(
        t0 = Normal(t0_guess, something(t0_width, time_scale / 4)),
        log_v = Normal(log(v_guess), log_v_width),
        b = Normal(b_center, b_width),
        log_tau0 = Normal(log(tau0_guess), log_tau0_width),
        log_lambda0 = Normal(log(lambda0_guess), log_lambda0_width),
        alpha = Normal(alpha_center, alpha_width),
        log_sigma_y = Normal(log(sigma_guess), log_sigma_width),
        log_sigma_x_plus = Normal(log(sigma_guess), log_sigma_width),
        log_sigma_x_minus = Normal(log(sigma_guess), log_sigma_width),
    )
end

function valid_quadratic_limb_darkening(u1, u2)
    u1 >= 0 && u1 + u2 <= 1 && u1 + 2u2 >= 0
end

function valid_positive_scales(values...)
    all(>(0), values)
end

function _dip_duration_guess(time, flux, baseline, depth)
    threshold = baseline - depth / 2
    in_dip = findall(<=(threshold), flux)

    if length(in_dip) >= 2
        duration = maximum(time[in_dip]) - minimum(time[in_dip])
        duration > 0 && return duration
    end

    time_span = maximum(time) - minimum(time)
    time_span > 0 ? time_span / 5 : 1.0
end
