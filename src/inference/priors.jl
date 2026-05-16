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

function valid_quadratic_limb_darkening(u1, u2)
    u1 >= 0 && u1 + u2 <= 1 && u1 + 2u2 >= 0
end

function valid_positive_scales(values...)
    all(>(0), values)
end
