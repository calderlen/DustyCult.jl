# DustyCult

julia tools for modeling transient stellar dimming from circumstellar dust occultation

code structure
- `src/star/`: stellar parameters, limb darkening, and stellar grids
- `src/occulter/`: dust-screen shapes, opacity laws, and sky-plane motion
- `src/model.jl`: forward flux models
- `src/lightcurve.jl`: light-curve containers, magnitude conversions, and data readers
- `src/inference/`: likelihoods, priors, and Turing sampling support
- `src/posterior.jl`: posterior sample tables and posterior predictive summaries

From the repository root:

```julia
using Pkg
Pkg.activate(".")
using DustyCult
```
