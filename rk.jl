#defines the step for RK1
function step_rk1(f, t, y, h)
    k1 = f(t, y)
    return y .+ h .* k1
end

#main solver for RK1
function solve_rk1(f, t_span, y0, h)
    t0, tf = t_span
    N = ceil(Int, (tf - t0) / h)
    t_grid = range(t0, length=N+1, step=h)
    Y = zeros(length(t_grid), length(y0))
    Y[1,:] = y0
    for n in 1:N
        Y[n+1,:] = step_rk1(f, t_grid[n], Y[n,:], h)
    end
    return t_grid, Y
end

#defines the step for RK2
function step_rk2(f, t, y, h)
    k1 = f(t, y)
    k2 = f(t + h, y .+ h .* k1)
    return y .+ (h / 2) .* (k1 .+ k2)
end

#main solver for RK2
function solve_rk2(f, t_span, y0, h)
    t0, tf = t_span
    N = ceil(Int, (tf - t0) / h)
    t_grid = range(t0, length=N+1, step=h)
    Y = zeros(length(t_grid), length(y0))
    Y[1,:] = y0
    for n in 1:N
        Y[n+1,:] = step_rk2(f, t_grid[n], Y[n,:], h)
    end
    return t_grid, Y
end

#defines the step for RK3
function step_rk3(f, t, y, h)
    k1 = f(t, y)
    k2 = f(t + h / 2, y .+ (h / 2) .* k1)
    k3 = f(t + h, y .+ h .* (-k1 .+ 2 * k2))
    return y .+ (h / 6) .* (k1 .+ 4 * k2 .+ k3)
end

#main solver for RK3
function solve_rk3(f, t_span, y0, h)
    t0, tf = t_span
    N = ceil(Int, (tf - t0) / h)
    t_grid = range(t0, length=N+1, step=h)
    Y = zeros(length(t_grid), length(y0))
    Y[1,:] = y0
    for n in 1:N
        Y[n+1,:] = step_rk3(f, t_grid[n], Y[n,:], h)
    end
    return t_grid, Y
end

#defines the step for RK4
function step_rk4(f, t, y, h)
    k1 = f(t, y)
    k2 = f(t + h / 2, y .+ (h / 2) .* k1)
    k3 = f(t + h / 2, y .+ (h / 2) .* k2)
    k4 = f(t + h, y .+ h .* k3)
    return y .+ (h / 6) .* (k1 .+ 2 * k2 .+ 2 * k3 .+ k4)
end

#main solver for RK4
function solve_rk4(f, t_span, y0, h)
    t0, tf = t_span
    N = ceil(Int, (tf - t0) / h)
    t_grid = range(t0, length=N+1, step=h)
    Y = zeros(length(t_grid), length(y0))
    Y[1,:] = y0
    for n in 1:N
        Y[n+1,:] = step_rk4(f, t_grid[n], Y[n,:], h)
    end
    return t_grid, Y
end

#defines the step for RK5
function step_rk5(f, t, y, h)
    k1 = f(t, y)
    k2 = f(t + h / 4, y .+ (h / 4) .* k1)
    k3 = f(t + h / 4, y .+ (h / 8) .* (k1 .+ k2))
    k4 = f(t + h / 2, y .+ (h / 2) .* k3)
    k5 = f(t + 3 * h / 4, y .+ (h / 16) .* (3 * k1 .+ 9 * k4))
    k6 = f(t + h, y .+ (h / 7) .* (2 * k1 .+ 3 * k2 .+ 4 * k4 .- 12 * k3))
    return y .+ (h / 90) .* (7 * k1 .+ 32 * k3 .+ 12 * k4 .+ 32 * k5 .+ 7 * k6)
end

#main solver for RK5
function solve_rk5(f, t_span, y0, h)
    t0, tf = t_span
    N = ceil(Int, (tf - t0) / h)
    t_grid = range(t0, length=N+1, step=h)
    Y = zeros(length(t_grid), length(y0))
    Y[1,:] = y0
    for n in 1:N
        Y[n+1,:] = step_rk5(f, t_grid[n], Y[n,:], h)
    end
    return t_grid, Y
end