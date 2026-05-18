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

## Fit CLI

Run a light-curve fit from JSON and write Python-readable Parquet artifacts:

```sh
julia --project=. -e 'using Pkg; Pkg.resolve(); Pkg.instantiate()'
julia --project=. scripts/fit_lightcurve.jl --config run.json --out output/run_001
```

The JSON config declares the light-curve file, bandpass wavelengths, stellar
parameters, grid size, sampler settings, priors, and posterior predictive draw
count. Supported light-curve formats are `relative_flux`, `delta_magnitude`,
and `magnitude`; magnitude configs must include `reference_magnitude`.

The output directory contains:
- `samples.parquet`
- `sampled_samples.parquet`
- `lightcurve.parquet`
- `predictive_intervals.parquet`
- `predictive_fluxes.parquet`
- `manifest.json`

Python applications should call the CLI directly with `subprocess` and manage
their own workflow-specific state around the resulting artifacts.
