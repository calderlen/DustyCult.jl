struct Bandpass
    wavelength
end

struct LightCurve
    time
    relative_flux
    relative_flux_error
    band
    wavelength
end

delta_mag_to_relative_flux(delta_mag::Real) = 10.0^(-0.4 * delta_mag)

relative_flux_to_delta_mag(relative_flux::Real) = -2.5 * log10(relative_flux)

mag_error_to_relative_flux_error(mag_error::Real, relative_flux::Real) = log(10) / 2.5 * relative_flux * mag_error

function lightcurve_from_delta_magnitudes(time, delta_mag, mag_error, band, bandpass)
    time = collect(time)
    delta_mag = collect(delta_mag)
    mag_error = collect(mag_error)
    band = string.(collect(band))

    relative_flux = delta_mag_to_relative_flux.(delta_mag)
    relative_flux_error = mag_error_to_relative_flux_error.(mag_error, relative_flux)
    wavelength = [bandpass[name].wavelength for name in band]

    LightCurve(time, relative_flux, relative_flux_error, band, wavelength)
end

function lightcurve_from_relative_fluxes(time, relative_flux, relative_flux_error, band, bandpass)
    time = collect(time)
    relative_flux = collect(relative_flux)
    relative_flux_error = collect(relative_flux_error)
    band = string.(collect(band))
    wavelength = [bandpass[name].wavelength for name in band]

    LightCurve(time, relative_flux, relative_flux_error, band, wavelength)
end

function lightcurve_from_magnitudes(time, mag, mag_error, band, bandpass, reference_magnitude)
    mag = collect(mag)
    band = string.(collect(band))
    delta_mag = [
        mag[i] - reference_magnitude[band[i]]
        for i in eachindex(mag)
    ]

    lightcurve_from_delta_magnitudes(time, delta_mag, mag_error, band, bandpass)
end

"""
Read a delimited light-curve file containing relative fluxes and flux errors.
"""
function read_relative_flux_lightcurve(path, bandpass; time_col = :time, relative_flux_col = :relative_flux, relative_flux_error_col = :relative_flux_error, band_col = :band)
    row = _read_delimited_rows(path)
    header = row[1]
    data = row[2:end]

    time_idx = _column_index(header, time_col)
    relative_flux_idx = _column_index(header, relative_flux_col)
    relative_flux_error_idx = _column_index(header, relative_flux_error_col)
    band_idx = _column_index(header, band_col)

    time = Float64[]
    relative_flux = Float64[]
    relative_flux_error = Float64[]
    band = String[]

    for entry in data
        push!(time, parse(Float64, entry[time_idx]))
        push!(relative_flux, parse(Float64, entry[relative_flux_idx]))
        push!(relative_flux_error, parse(Float64, entry[relative_flux_error_idx]))
        push!(band, entry[band_idx])
    end

    lightcurve_from_relative_fluxes(time, relative_flux, relative_flux_error, band, bandpass)
end

function read_delta_magnitude_lightcurve(path, bandpass; time_col = :time, delta_mag_col = :delta_mag, mag_error_col = :mag_error, band_col = :band)
    row = _read_delimited_rows(path)
    header = row[1]
    data = row[2:end]

    time_idx = _column_index(header, time_col)
    delta_mag_idx = _column_index(header, delta_mag_col)
    mag_error_idx = _column_index(header, mag_error_col)
    band_idx = _column_index(header, band_col)

    time = Float64[]
    delta_mag = Float64[]
    mag_error = Float64[]
    band = String[]

    for entry in data
        push!(time, parse(Float64, entry[time_idx]))
        push!(delta_mag, parse(Float64, entry[delta_mag_idx]))
        push!(mag_error, parse(Float64, entry[mag_error_idx]))
        push!(band, entry[band_idx])
    end

    lightcurve_from_delta_magnitudes(time, delta_mag, mag_error, band, bandpass)
end

function read_magnitude_lightcurve(path, bandpass, reference_magnitude; time_col = :time, mag_col = :mag, mag_error_col = :mag_error, band_col = :band)
    row = _read_delimited_rows(path)
    header = row[1]
    data = row[2:end]

    time_idx = _column_index(header, time_col)
    mag_idx = _column_index(header, mag_col)
    mag_error_idx = _column_index(header, mag_error_col)
    band_idx = _column_index(header, band_col)

    time = Float64[]
    mag = Float64[]
    mag_error = Float64[]
    band = String[]

    for entry in data
        push!(time, parse(Float64, entry[time_idx]))
        push!(mag, parse(Float64, entry[mag_idx]))
        push!(mag_error, parse(Float64, entry[mag_error_idx]))
        push!(band, entry[band_idx])
    end

    lightcurve_from_magnitudes(time, mag, mag_error, band, bandpass, reference_magnitude)
end

function fit_inputs(lightcurve::LightCurve)
    (;
        times = lightcurve.time,
        wavelengths = lightcurve.wavelength,
        observed = lightcurve.relative_flux,
        errors = lightcurve.relative_flux_error,
    )
end

function _read_delimited_rows(path)
    raw_line = readlines(path)
    line = [
        strip(entry)
        for entry in raw_line
        if !isempty(strip(entry)) && !startswith(strip(entry), "#")
    ]
    delimiter = _detect_delimiter(first(line))

    [_split_delimited_line(entry, delimiter) for entry in line]
end

function _detect_delimiter(line)
    occursin(",", line) && return ','
    occursin("\t", line) && return '\t'

    nothing
end

function _split_delimited_line(line, delimiter)
    if delimiter === nothing
        split(strip(line))
    else
        strip.(split(line, delimiter))
    end
end

function _column_index(header, column)
    column_name = string(column)
    index = findfirst(==(column_name), header)

    index === nothing && throw(ArgumentError("missing required column $column_name; available columns: $(join(header, ", "))"))

    index
end
