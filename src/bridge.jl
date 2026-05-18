using DataFrames
using Dates
using JSON3
using Parquet2
using Random

const FIT_OUTPUT_SCHEMA_VERSION = "1.0"
struct MissingConfigSentinel end
const _MISSING_CONFIG = MissingConfigSentinel()

"""
Read a DustyCult fit config from a JSON file.
"""
read_fit_config(path) = JSON3.read(read(path, String))

"""
Run a DustyCult light-curve fit from a JSON config file and write Parquet artifacts.
"""
function run_fit_config(config_path::AbstractString, output_dir::AbstractString)
    config = read_fit_config(config_path)
    fit = fit_lightcurve_from_config(config; config_path=config_path)
    paths = write_fit_outputs(fit, output_dir, config; config_path=config_path)

    (; fit, output_dir=abspath(output_dir), paths)
end

"""
Run a DustyCult light-curve fit from an already parsed config object.
"""
function fit_lightcurve_from_config(config; config_path=nothing)
    lightcurve_config = _config_get(config, "lightcurve")
    lightcurve_path = _resolve_config_path(_config_get(lightcurve_config, "path"), config_path)
    lightcurve_format = Symbol(string(_config_get(lightcurve_config, "format")))
    columns = _config_get(lightcurve_config, "columns", (;))
    bandpass = build_bandpass(_config_get(config, "bandpass"))
    star = build_star_params(_config_get(config, "star", (;)))
    grid = build_stellar_grid(_config_get(config, "grid", (;)), star)
    prior_kwargs = _named_tuple(_config_get(config, "prior_kwargs", (;)))
    sampling_kwargs = _sampling_kwargs(_config_get(config, "sampling", (;)))
    seed = _config_get(_config_get(config, "sampling", (;)), "seed", nothing)
    n_predictive_draws = _n_predictive_draws(config)
    rng = seed === nothing ? Random.default_rng() : Random.MersenneTwister(seed)

    if seed !== nothing
        Random.seed!(seed)
    end

    fit_kwargs = (;
        prior_kwargs,
        n_predictive_draws,
        rng,
        sampling_kwargs...,
    )

    if lightcurve_format === :relative_flux
        fit_relative_flux_lightcurve(
            lightcurve_path,
            bandpass,
            grid,
            star;
            time_col=_column_name(columns, "time", "time"),
            relative_flux_col=_column_name(columns, "relative_flux", "relative_flux"),
            relative_flux_error_col=_column_name(columns, "relative_flux_error", "relative_flux_error"),
            band_col=_column_name(columns, "band", "band"),
            fit_kwargs...,
        )
    elseif lightcurve_format === :delta_magnitude
        fit_delta_magnitude_lightcurve(
            lightcurve_path,
            bandpass,
            grid,
            star;
            time_col=_column_name(columns, "time", "time"),
            delta_mag_col=_column_name(columns, "delta_mag", "delta_mag"),
            mag_error_col=_column_name(columns, "mag_error", "mag_error"),
            band_col=_column_name(columns, "band", "band"),
            fit_kwargs...,
        )
    elseif lightcurve_format === :magnitude
        reference_magnitude = _reference_magnitude(config, lightcurve_config)

        fit_magnitude_lightcurve(
            lightcurve_path,
            bandpass,
            reference_magnitude,
            grid,
            star;
            time_col=_column_name(columns, "time", "time"),
            mag_col=_column_name(columns, "mag", "mag"),
            mag_error_col=_column_name(columns, "mag_error", "mag_error"),
            band_col=_column_name(columns, "band", "band"),
            fit_kwargs...,
        )
    else
        throw(ArgumentError("unsupported lightcurve.format: $(lightcurve_format)"))
    end
end

"""
Build a bandpass dictionary from config mapping band names to wavelengths.
"""
function build_bandpass(config)
    isempty(_config_pairs(config)) && throw(ArgumentError("bandpass config cannot be empty"))

    Dict(
        string(name) => Bandpass(Float64(_config_get(entry, "wavelength")))
        for (name, entry) in _config_pairs(config)
    )
end

"""
Build stellar parameters from config, using `StarParams` defaults for omitted fields.
"""
function build_star_params(config)
    StarParams(
        R=Float64(_config_get(config, "R", 1.0)),
        I0=Float64(_config_get(config, "I0", 1.0)),
        u1=Float64(_config_get(config, "u1", 0.0)),
        u2=Float64(_config_get(config, "u2", 0.0)),
    )
end

"""
Build the stellar grid from config and the selected star.
"""
function build_stellar_grid(config, star::StarParams)
    stellar_grid(star.R; n=Int(_config_get(config, "n", 101)))
end

"""
Write Python-readable fit artifacts into `output_dir`.
"""
function write_fit_outputs(fit, output_dir::AbstractString, config; config_path=nothing)
    mkpath(output_dir)

    paths = Dict{String,Any}()
    paths["samples"] = _write_parquet(joinpath(output_dir, "samples.parquet"), _draw_table(fit.samples))
    paths["sampled_samples"] = _write_parquet(joinpath(output_dir, "sampled_samples.parquet"), _draw_table(fit.sampled_samples))
    paths["lightcurve"] = _write_parquet(joinpath(output_dir, "lightcurve.parquet"), lightcurve_output_table(fit.lightcurve))

    if fit.predictive_intervals !== nothing
        paths["predictive_intervals"] = _write_parquet(
            joinpath(output_dir, "predictive_intervals.parquet"),
            predictive_intervals_table(fit),
        )
    else
        paths["predictive_intervals"] = nothing
    end

    if fit.predictive_fluxes !== nothing
        paths["predictive_fluxes"] = _write_parquet(
            joinpath(output_dir, "predictive_fluxes.parquet"),
            predictive_fluxes_table(fit),
        )
    else
        paths["predictive_fluxes"] = nothing
    end

    manifest_path = joinpath(output_dir, "manifest.json")
    write_fit_manifest(manifest_path, config, paths; config_path=config_path, output_dir=output_dir)
    paths["manifest"] = abspath(manifest_path)

    paths
end

function lightcurve_output_table(lightcurve::LightCurve)
    DataFrame(
        point_id=collect(1:length(lightcurve.time)),
        time=lightcurve.time,
        band=lightcurve.band,
        wavelength=lightcurve.wavelength,
        relative_flux=lightcurve.relative_flux,
        relative_flux_error=lightcurve.relative_flux_error,
    )
end

function predictive_intervals_table(fit)
    table = DataFrame(
        point_id=collect(1:length(fit.lightcurve.time)),
        time=fit.lightcurve.time,
        band=fit.lightcurve.band,
        wavelength=fit.lightcurve.wavelength,
        observed=fit.lightcurve.relative_flux,
        error=fit.lightcurve.relative_flux_error,
    )

    for name in propertynames(fit.predictive_intervals)
        table[!, name] = getproperty(fit.predictive_intervals, name)
    end

    table
end

function predictive_fluxes_table(fit)
    predictions = fit.predictive_fluxes
    n_draws, n_points = size(predictions)
    row_count = n_draws * n_points
    table = DataFrame(
        draw_id=Vector{Int}(undef, row_count),
        point_id=Vector{Int}(undef, row_count),
        time=Vector{Float64}(undef, row_count),
        band=Vector{String}(undef, row_count),
        wavelength=Vector{Float64}(undef, row_count),
        predictive_flux=Vector{Float64}(undef, row_count),
    )

    row = 1
    for draw in 1:n_draws
        for point in 1:n_points
            table.draw_id[row] = draw
            table.point_id[row] = point
            table.time[row] = fit.lightcurve.time[point]
            table.band[row] = fit.lightcurve.band[point]
            table.wavelength[row] = fit.lightcurve.wavelength[point]
            table.predictive_flux[row] = predictions[draw, point]
            row += 1
        end
    end

    table
end

function write_fit_manifest(path, config, artifact_paths; config_path=nothing, output_dir=nothing)
    manifest = Dict{String,Any}(
        "schema_version" => FIT_OUTPUT_SCHEMA_VERSION,
        "package" => Dict{String,Any}(
            "name" => "DustyCult",
            "version" => _package_version(),
        ),
        "created_at" => string(Dates.now(Dates.UTC)) * "Z",
        "config_path" => config_path === nothing ? nothing : abspath(config_path),
        "output_dir" => output_dir === nothing ? dirname(abspath(path)) : abspath(output_dir),
        "config" => config,
        "artifacts" => artifact_paths,
    )

    open(path, "w") do io
        JSON3.pretty(io, manifest)
        println(io)
    end

    abspath(path)
end

function parse_fit_cli_args(args)
    config_path = nothing
    output_dir = nothing
    help = false
    index = 1

    while index <= length(args)
        arg = args[index]

        if arg in ("-h", "--help")
            help = true
        elseif arg == "--config"
            index += 1
            index <= length(args) || throw(ArgumentError("--config requires a value"))
            config_path = args[index]
        elseif startswith(arg, "--config=")
            config_path = arg[length("--config=")+1:end]
        elseif arg == "--out"
            index += 1
            index <= length(args) || throw(ArgumentError("--out requires a value"))
            output_dir = args[index]
        elseif startswith(arg, "--out=")
            output_dir = arg[length("--out=")+1:end]
        else
            throw(ArgumentError("unknown argument: $arg"))
        end

        index += 1
    end

    if !help
        config_path === nothing && throw(ArgumentError("--config is required"))
        output_dir === nothing && throw(ArgumentError("--out is required"))
    end

    (; config=config_path, out=output_dir, help)
end

function fit_lightcurve_cli(args=ARGS)
    try
        parsed = parse_fit_cli_args(args)

        if parsed.help
            print(_fit_lightcurve_usage())
            return 0
        end

        result = run_fit_config(parsed.config, parsed.out)
        println("Wrote DustyCult fit outputs to $(result.output_dir)")
        return 0
    catch err
        showerror(stderr, err)
        println(stderr)
        println(stderr, _fit_lightcurve_usage())
        return 1
    end
end

function _draw_table(samples)
    table = copy(samples)
    insertcols!(table, 1, :draw_id => collect(1:nrow(table)))

    table
end

function _write_parquet(path, table)
    Parquet2.writefile(path, table)

    abspath(path)
end

function _n_predictive_draws(config)
    predictive_config = _config_get(config, "posterior_predictive", (;))
    predictive_config === nothing && return nothing

    _config_get(predictive_config, "n_draws", 200)
end

function _sampling_kwargs(config)
    kwargs = Dict{Symbol,Any}()

    for name in ("n_samples", "n_adapt", "n_chains", "progress")
        if _config_has(config, name)
            kwargs[Symbol(name)] = _config_get(config, name)
        end
    end

    (; kwargs...)
end

function _reference_magnitude(config, lightcurve_config)
    reference_magnitude = if _config_has(lightcurve_config, "reference_magnitude")
        _config_get(lightcurve_config, "reference_magnitude")
    elseif _config_has(config, "reference_magnitude")
        _config_get(config, "reference_magnitude")
    else
        throw(ArgumentError("magnitude light curves require reference_magnitude"))
    end

    Dict(string(name) => Float64(value) for (name, value) in _config_pairs(reference_magnitude))
end

function _named_tuple(config)
    pairs = _config_pairs(config)
    isempty(pairs) && return (;)

    keys = Tuple(Symbol(name) for (name, _) in pairs)
    values = Tuple(value for (_, value) in pairs)

    NamedTuple{keys}(values)
end

function _column_name(columns, key, default)
    string(_config_get(columns, key, default))
end

function _config_get(config, key::AbstractString, default=_MISSING_CONFIG)
    if _config_has(config, key)
        return _config_value(config, key)
    end

    default === _MISSING_CONFIG && throw(ArgumentError("missing required config key: $key"))

    default
end

function _config_has(config, key::AbstractString)
    for candidate in (key, Symbol(key))
        try
            haskey(config, candidate) && return true
        catch
        end
    end

    config isa NamedTuple && return hasproperty(config, Symbol(key))

    false
end

function _config_value(config, key::AbstractString)
    for candidate in (key, Symbol(key))
        try
            haskey(config, candidate) && return config[candidate]
        catch
        end
    end

    config isa NamedTuple && return getproperty(config, Symbol(key))

    throw(ArgumentError("missing required config key: $key"))
end

function _config_pairs(config)
    config === nothing && return Pair{String,Any}[]

    if config isa AbstractDict
        return [string(key) => value for (key, value) in pairs(config)]
    elseif config isa NamedTuple
        return [string(key) => getproperty(config, key) for key in keys(config)]
    else
        return [string(key) => value for (key, value) in pairs(config)]
    end
end

function _resolve_config_path(path, config_path)
    path_string = string(path)
    resolved = if isabspath(path_string) || config_path === nothing
        path_string
    else
        joinpath(dirname(abspath(config_path)), path_string)
    end

    abspath(resolved)
end

function _package_version()
    version = Base.pkgversion(@__MODULE__)

    version === nothing ? nothing : string(version)
end

function _fit_lightcurve_usage()
    """
    Usage:
      julia --project=. scripts/fit_lightcurve.jl --config run.json --out output/run_001

    Required:
      --config PATH    JSON fit configuration
      --out PATH       Output directory for Parquet artifacts and manifest.json
    """
end
