#!/usr/bin/env julia

import Pkg

Pkg.activate(joinpath(@__DIR__, ".."))

using DustyCult

exit(DustyCult.fit_lightcurve_cli(ARGS))
