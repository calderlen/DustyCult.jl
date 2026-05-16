function mu(x, y, R)
    r2 = x^2 + y^2
    r2 > R^2 && return 0.0
    sqrt(max(0.0, 1 - r2 / R^2))
end

function limb_darkening(x, y, R, u1, u2)
    x^2 + y^2 > R^2 && return 0.0
    mu_local = mu(x, y, R)
    1 - u1*(1-mu_local) - u2*(1-mu_local)^2
end

function I(x, y, star::StarParams)
    star.I0 * limb_darkening(x, y, star.R, star.u1, star.u2)
end
