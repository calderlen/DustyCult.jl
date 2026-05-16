gaussian(x, sigma) = exp(-0.5 * (x / sigma)^2)

function gaussian2d(x, y, sigma_x, sigma_y)
    gaussian(x, sigma_x) * gaussian(y, sigma_y)
end

function two_gaussian2d(x, y, sigma_y, sigma_x_1, sigma_x_2)
    0.5 * (gaussian2d(x, y, sigma_x_1, sigma_y) + gaussian2d(x, y, sigma_x_2, sigma_y))
end
