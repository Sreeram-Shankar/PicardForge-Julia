#defines the SDIRK step with Gauss-Seidel relaxation
function step_sdirk(f, t, y, h, A, b, c; sweeps=12, tol=1e-10)
    s = length(b)
    n = length(y)

    #defines the initial guesses for all stages
    Y = [copy(y) for _ in 1:s]

    #implements Gauss-Seidel relaxation
    for _ in 1:sweeps
        Y_old = deepcopy(Y)
        for i in 1:s
            rhs = zeros(n)
            for j in 1:s
                rhs .+= A[i,j] * f(t + c[j]*h, Y[j])
            end
            Y[i] = y .+ h .* rhs
        end
        diff_norm = sqrt(sum(sum((Y[i] - Y_old[i]).^2) for i in 1:s))
        if diff_norm < tol
            break
        end
    end

    #computes the final state update
    K = [f(t + c[i]*h, Y[i]) for i in 1:s]
    y_next = y .+ h .* sum(b[i] .* K[i] for i in 1:s)
    return y_next
end


#solves the nonlinear system of equations with a Gauss-Seidel relaxation SDIRK2
function solve_sdirk2(f, t_span, y0, h; sweeps=12, tol=1e-10)
    gamma = 1.0 - 1.0/sqrt(2.0)
    A = [
        gamma  0.0;
        1.0 - gamma  gamma
    ]
    b = [1.0 - gamma, gamma]
    c = [gamma, 1.0]

    t0, tf = t_span
    N = ceil(Int, (tf - t0)/h)
    t_grid = range(t0, length=N+1, step=h)
    Y = zeros(length(t_grid), length(y0))
    Y[1,:] = y0

    for n in 1:N
        Y[n+1,:] = step_sdirk(f, t_grid[n], Y[n,:], h, A, b, c; sweeps=sweeps, tol=tol)
    end
    return t_grid, Y
end


#solves the nonlinear system of equations with a Gauss-Seidel relaxation SDIRK3
function solve_sdirk3(f, t_span, y0, h; sweeps=12, tol=1e-10)
    gamma = 0.435866521508459
    A = [
        gamma           0.0               0.0;
        0.2820667395    gamma             0.0;
        1.208496649    -0.644363171       gamma
    ]
    b = [1.208496649, -0.644363171, gamma]
    c = [gamma, 0.7179332605, 1.0]

    t0, tf = t_span
    N = ceil(Int, (tf - t0)/h)
    t_grid = range(t0, length=N+1, step=h)
    Y = zeros(length(t_grid), length(y0))
    Y[1,:] = y0

    for n in 1:N
        Y[n+1,:] = step_sdirk(f, t_grid[n], Y[n,:], h, A, b, c; sweeps=sweeps, tol=tol)
    end
    return t_grid, Y
end


#solves the nonlinear system of equations with a Gauss-Seidel relaxation SDIRK4
function solve_sdirk4(f, t_span, y0, h; sweeps=12, tol=1e-10)
    gamma = 0.572816062482135
    A = [
        gamma           0.0           0.0             0.0;
       -0.6557110092    gamma         0.0             0.0;
        0.757184241     0.237758128   gamma           0.0;
        0.155416858     0.701913790   0.142669351     gamma
    ]
    b = [0.155416858, 0.701913790, 0.142669351, gamma]
    c = [gamma, 0.344, 0.995, 1.0]

    t0, tf = t_span
    N = ceil(Int, (tf - t0)/h)
    t_grid = range(t0, length=N+1, step=h)
    Y = zeros(length(t_grid), length(y0))
    Y[1,:] = y0

    for n in 1:N
        Y[n+1,:] = step_sdirk(f, t_grid[n], Y[n,:], h, A, b, c; sweeps=sweeps, tol=tol)
    end
    return t_grid, Y
end