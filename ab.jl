#defines the step for AB2
function step_ab2(f, t_n, y_n, y_prev, f_n, f_prev, h)
    return y_n .+ h .* ((3/2)*f_n .- (1/2)*f_prev)
end

#main solver for AB2
function solve_ab2(f, t_span, y0, h)
    t0, tf = t_span
    N = ceil(Int, (tf - t0) / h)
    t_grid = range(t0, length=N+1, step=h)
    Y = zeros(length(t_grid), length(y0))
    F = Vector{Vector{Float64}}(undef, N+1)
    
    Y[1,:] = y0
    F[1] = f(t_grid[1], Y[1,:])
    
    #uses RK2 for first step
    k1 = f(t_grid[1], Y[1,:])
    k2 = f(t_grid[1] + h, Y[1,:] .+ h .* k1)
    Y[2,:] = Y[1,:] .+ 0.5 * h .* (k1 .+ k2)
    F[2] = f(t_grid[2], Y[2,:])
    
    for n in 2:N
        Y[n+1,:] = step_ab2(f, t_grid[n], Y[n,:], Y[n-1,:], F[n], F[n-1], h)
        F[n+1] = f(t_grid[n+1], Y[n+1,:])
    end
    
    return t_grid, Y
end


#defines the step for AB3
function step_ab3(f, t_n, y_n, y_prev, y_prev2, f_n, f_prev, f_prev2, h)
    return y_n .+ h .* ((23/12)*f_n .- (16/12)*f_prev .+ (5/12)*f_prev2)
end

#main solver for AB3
function solve_ab3(f, t_span, y0, h)
    t0, tf = t_span
    N = ceil(Int, (tf - t0) / h)
    t_grid = range(t0, length=N+1, step=h)
    Y = zeros(length(t_grid), length(y0))
    F = Vector{Vector{Float64}}(undef, N+1)
    
    Y[1,:] = y0
    F[1] = f(t_grid[1], Y[1,:])
    
    #uses RK3 for first two steps
    k1 = f(t_grid[1], Y[1,:])
    k2 = f(t_grid[1] + h/2, Y[1,:] .+ (h/2) .* k1)
    k3 = f(t_grid[1] + h, Y[1,:] .+ h .* (-k1 .+ 2*k2))
    Y[2,:] = Y[1,:] .+ (h/6) .* (k1 .+ 4*k2 .+ k3)
    F[2] = f(t_grid[2], Y[2,:])
    
    k1 = f(t_grid[2], Y[2,:])
    k2 = f(t_grid[2] + h/2, Y[2,:] .+ (h/2) .* k1)
    k3 = f(t_grid[2] + h, Y[2,:] .+ h .* (-k1 .+ 2*k2))
    Y[3,:] = Y[2,:] .+ (h/6) .* (k1 .+ 4*k2 .+ k3)
    F[3] = f(t_grid[3], Y[3,:])
    
    for n in 3:N
        Y[n+1,:] = step_ab3(f, t_grid[n], Y[n,:], Y[n-1,:], Y[n-2,:], F[n], F[n-1], F[n-2], h)
        F[n+1] = f(t_grid[n+1], Y[n+1,:])
    end
    
    return t_grid, Y
end


#defines the step for AB4
function step_ab4(f, t_n, y_n, y_prev, y_prev2, y_prev3, f_n, f_prev, f_prev2, f_prev3, h)
    return y_n .+ h .* ((55/24)*f_n .- (59/24)*f_prev .+ (37/24)*f_prev2 .- (9/24)*f_prev3)
end

#main solver for AB4
function solve_ab4(f, t_span, y0, h)
    t0, tf = t_span
    N = ceil(Int, (tf - t0) / h)
    t_grid = range(t0, length=N+1, step=h)
    Y = zeros(length(t_grid), length(y0))
    F = Vector{Vector{Float64}}(undef, N+1)
    
    Y[1,:] = y0
    F[1] = f(t_grid[1], Y[1,:])
    
    #uses RK4 for first three steps
    for i in 1:3
        k1 = f(t_grid[i], Y[i,:])
        k2 = f(t_grid[i] + h/2, Y[i,:] .+ (h/2) .* k1)
        k3 = f(t_grid[i] + h/2, Y[i,:] .+ (h/2) .* k2)
        k4 = f(t_grid[i] + h, Y[i,:] .+ h .* k3)
        Y[i+1,:] = Y[i,:] .+ (h/6) .* (k1 .+ 2*k2 .+ 2*k3 .+ k4)
        F[i+1] = f(t_grid[i+1], Y[i+1,:])
    end
    
    for n in 4:N
        Y[n+1,:] = step_ab4(f, t_grid[n], Y[n,:], Y[n-1,:], Y[n-2,:], Y[n-3,:],
                            F[n], F[n-1], F[n-2], F[n-3], h)
        F[n+1] = f(t_grid[n+1], Y[n+1,:])
    end
    
    return t_grid, Y
end


#defines the step for AB5
function step_ab5(f, t_n, y_n, y_prev, y_prev2, y_prev3, y_prev4, f_n, f_prev, f_prev2, f_prev3, f_prev4, h)
    return y_n .+ h .* ((1901/720)*f_n .- (2774/720)*f_prev .+ (2616/720)*f_prev2
                       .- (1274/720)*f_prev3 .+ (251/720)*f_prev4)
end

#main solver for AB5
function solve_ab5(f, t_span, y0, h)
    t0, tf = t_span
    N = ceil(Int, (tf - t0) / h)
    t_grid = range(t0, length=N+1, step=h)
    Y = zeros(length(t_grid), length(y0))
    F = Vector{Vector{Float64}}(undef, N+1)
    
    Y[1,:] = y0
    F[1] = f(t_grid[1], Y[1,:])
    
    #uses RK4 for first four steps
    for i in 1:4
        k1 = f(t_grid[i], Y[i,:])
        k2 = f(t_grid[i] + h/2, Y[i,:] .+ (h/2) .* k1)
        k3 = f(t_grid[i] + h/2, Y[i,:] .+ (h/2) .* k2)
        k4 = f(t_grid[i] + h, Y[i,:] .+ h .* k3)
        Y[i+1,:] = Y[i,:] .+ (h/6) .* (k1 .+ 2*k2 .+ 2*k3 .+ k4)
        F[i+1] = f(t_grid[i+1], Y[i+1,:])
    end
    
    for n in 5:N
        Y[n+1,:] = step_ab5(f, t_grid[n], Y[n,:], Y[n-1,:], Y[n-2,:], Y[n-3,:], Y[n-4,:],
                            F[n], F[n-1], F[n-2], F[n-3], F[n-4], h)
        F[n+1] = f(t_grid[n+1], Y[n+1,:])
    end
    
    return t_grid, Y
end