Base.@kwdef struct DustOccultationParams{T}
    t0::T = 0.0
    v::T = 1.0
    b::T = 0.0
    tau0::T = 1.0
    lambda0::T = 1.0
    alpha::T = 0.0
    sigma_y::T = 0.25
    sigma_x_plus::T = 0.25
    sigma_x_minus::T = 0.25
end

function F0(grid, star::StarParams)
    intensities = I.(grid.x, grid.y, Ref(star))
    sum(intensities .* grid.weights)
end

function dust_transit_flux(grid, t, lambda, star::StarParams, dust::DustOccultationParams)
    intensities = I.(grid.x, grid.y, Ref(star))

    x_center, y_center = occulter_center(t, dust.t0, dust.v, dust.b)
    x_occ = grid.x .- x_center
    y_occ = grid.y .- y_center

    shape = two_gaussian2d.(
        x_occ,
        y_occ,
        dust.sigma_y,
        dust.sigma_x_plus,
        dust.sigma_x_minus,
    )

    tau = optical_depth.(shape, lambda, dust.tau0, dust.lambda0, dust.alpha)
    sum(intensities .* transmission.(tau) .* grid.weights) / F0(grid, star)
end
