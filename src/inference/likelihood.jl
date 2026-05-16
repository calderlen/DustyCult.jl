lnL(observed, modeled, errors) = -0.5 * sum(((obs - mod) / err)^2 + log(2*pi * err^2) for (obs, mod, err) in zip(observed, modeled, errors))
