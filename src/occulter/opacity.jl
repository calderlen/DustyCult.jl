function powerlaw_opacity(lambda, tau0, lambda0, alpha)
    tau0 * (lambda / lambda0)^(-alpha)
end

function optical_depth(shape, lambda, tau0, lambda0, alpha)
    shape * powerlaw_opacity(lambda, tau0, lambda0, alpha)
end

transmission(tau) = exp(-tau)
