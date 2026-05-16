function stellar_grid(R; n::Integer=101)
    xs = collect(range(-R, R; length=n))
    ys = collect(range(-R, R; length=n))
    dx = xs[2] - xs[1]
    dy = ys[2] - ys[1]

    x = [xi for _ in ys, xi in xs]
    y = [yi for yi in ys, _ in xs]
    mask = x.^2 .+ y.^2 .<= R^2
    weights = fill(dx * dy, n, n)

    (; x, y, mask, weights)
end
