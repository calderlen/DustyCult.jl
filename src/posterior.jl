using DataFrames
using MCMCChains
using Random
using Statistics: quantile

const SAMPLED_FIT_PARAMETERS = (
    :t0,
    :log_v,
    :b,
    :log_tau0,
    :log_lambda0,
    :alpha,
    :log_sigma_y,
    :log_sigma_x_plus,
    :log_sigma_x_minus,
)

const PHYSICAL_FIT_PARAMETERS = (
    :t0,
    :v,
    :b,
    :tau0,
    :lambda0,
    :alpha,
    :sigma_y,
    :sigma_x_plus,
    :sigma_x_minus,
)

const PHYSICAL_TO_SAMPLED_LOG = Dict(
    :v => :log_v,
    :tau0 => :log_tau0,
    :lambda0 => :log_lambda0,
    :sigma_y => :log_sigma_y,
    :sigma_x_plus => :log_sigma_x_plus,
    :sigma_x_minus => :log_sigma_x_minus,
)

const DEFAULT_FIT_LABELS = Dict{Symbol,String}(
    :t0 => "t0",
    :v => "v",
    :b => "b",
    :tau0 => "tau0",
    :lambda0 => "lambda0",
    :alpha => "alpha",
    :sigma_y => "sigma_y",
    :sigma_x_plus => "sigma_x_plus",
    :sigma_x_minus => "sigma_x_minus",
    :log_v => "log_v",
    :log_tau0 => "log_tau0",
    :log_lambda0 => "log_lambda0",
    :log_sigma_y => "log_sigma_y",
    :log_sigma_x_plus => "log_sigma_x_plus",
    :log_sigma_x_minus => "log_sigma_x_minus",
)

"""
Convert an MCMC chain into table. By default the sampled log-space scale parameters are exponentiated into the physical dust model parameters used by `DustOccultationParams`.
"""
function fit_samples_table(chain; scale=:physical, parameters=nothing, t0_offset=0.0)
    scale = Symbol(scale)
    scale in (:physical, :sampled) || throw(ArgumentError("scale must be :physical or :sampled"))
    requested = _fit_sample_parameters(scale, parameters)
    source = DataFrame(chain)
    samples = DataFrame()

    for parameter in requested
        if scale === :physical
            _insert_physical_parameter!(samples, source, parameter)
        else
            samples[!, parameter] = source[!, parameter]
        end
    end

    _offset_t0!(samples, t0_offset)

    samples
end

"""
Draw posterior samples from `chain` and evaluate model fluxes for each draw. Rows correspond to posterior draws and columns correspond to `times`.
"""
function posterior_predictive_fluxes(chain, grid, times, wavelengths, star; n_draws=200, rng=Random.default_rng())
    sample_table = fit_samples_table(chain; scale=:physical)
    n_samples = nrow(sample_table)
    draw_count = min(n_draws, n_samples)
    draw_indices = draw_count == n_samples ? collect(1:n_samples) : Random.randperm(rng, n_samples)[1:draw_count]

    first_dust = _dust_params_from_sample(sample_table, first(draw_indices))
    first_prediction = collect(dust_model_fluxes(grid, times, wavelengths, star, first_dust))
    predictions = Matrix{Float64}(undef, length(draw_indices), length(first_prediction))
    predictions[1, :] = first_prediction

    for row in 2:length(draw_indices)
        draw_index = draw_indices[row]
        dust = _dust_params_from_sample(sample_table, draw_index)
        predictions[row, :] = dust_model_fluxes(grid, times, wavelengths, star, dust)
    end

    predictions
end

"""
Summarize posterior predictive flux draws into 95%, 68%, and median intervals. `predictions` should have posterior draws in rows and prediction points in columns.
"""
function posterior_predictive_intervals(predictions)
    qs = [quantile(view(predictions, :, index), [0.025, 0.16, 0.5, 0.84, 0.975]) for index in axes(predictions, 2)]

    (;
        lower95=[q[1] for q in qs],
        lower68=[q[2] for q in qs],
        median=[q[3] for q in qs],
        upper68=[q[4] for q in qs],
        upper95=[q[5] for q in qs],
    )
end

function _fit_sample_parameters(scale, parameters)
    if parameters === nothing
        return collect(scale === :physical ? PHYSICAL_FIT_PARAMETERS : SAMPLED_FIT_PARAMETERS)
    elseif parameters isa Symbol || parameters isa AbstractString
        return [Symbol(parameters)]
    else
        return Symbol.(collect(parameters))
    end
end

function _insert_physical_parameter!(target, source, parameter)
    if haskey(PHYSICAL_TO_SAMPLED_LOG, parameter)
        sampled_parameter = PHYSICAL_TO_SAMPLED_LOG[parameter]
        target[!, parameter] = exp.(source[!, sampled_parameter])
    else
        target[!, parameter] = source[!, parameter]
    end

    target
end

function _offset_t0!(samples, t0_offset)
    t0_offset = Float64(t0_offset)
    if t0_offset != 0.0 && hasproperty(samples, :t0)
        samples[!, :t0] = samples[!, :t0] .+ t0_offset
    end

    samples
end

function _dust_params_from_sample(sample_table, row)
    DustOccultationParams(
        t0=sample_table[row, :t0],
        v=sample_table[row, :v],
        b=sample_table[row, :b],
        tau0=sample_table[row, :tau0],
        lambda0=sample_table[row, :lambda0],
        alpha=sample_table[row, :alpha],
        sigma_y=sample_table[row, :sigma_y],
        sigma_x_plus=sample_table[row, :sigma_x_plus],
        sigma_x_minus=sample_table[row, :sigma_x_minus],
    )
end
