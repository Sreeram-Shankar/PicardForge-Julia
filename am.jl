#solves the nonlinear system of equations with a Gauss-Seidel relaxation AM2
function solve_am2(f, t_span, y0, h; sweeps=12, tol=1e-10)
    t0, tf = t_span
    N = ceil(Int, (tf - t0)/h)
    t_grid = range(t0, length=N+1, step=h)
    Y = zeros(length(t_grid), length(y0))
    F = zeros(length(t_grid), length(y0))
    Y[1,:] = y0
    F[1,:] = f(t_grid[1], Y[1,:])

    #uses BE bootstrap
    Y[2,:] = Y[1,:] .+ h .* f(t_grid[2], Y[1,:])
    F[2,:] = f(t_grid[2], Y[2,:])

    for n in 2:N
        t_next = t_grid[n+1]
        f_n = F[n,:]
        y = copy(Y[n,:])

        #implements Gauss-Seidel relaxation
        for _ in 1:sweeps
            y_old = copy(y)
            f_next = f(t_next, y)
            y = Y[n,:] .+ h .* (0.5*f_next .+ 0.5*f_n)
            if norm(y - y_old) < tol
                break
            end
        end

        Y[n+1,:] = y
        F[n+1,:] = f_next
    end

    return t_grid, Y
end

#solves the nonlinear system of equations with a Gauss-Seidel relaxation AM3
function solve_am3(f, t_span, y0, h; sweeps=12, tol=1e-10)
    #uses AM2 to bootstrap 1 step
    _, Y = solve_am2(f, (t_span[1], t_span[1]+2*h), y0, h; sweeps=sweeps, tol=tol)
    t0, tf = t_span
    N = ceil(Int, (tf-t0)/h)
    t_grid = range(t0, length=N+1, step=h)

    #compute initial F
    F = zeros(length(t_grid), length(y0))
    for i in 1:3
        F[i,:] = f(t_grid[i], Y[i,:])
    end

    #defines the main AM3 solver
    for n in 3:N
        t_next = t_grid[n+1]
        f_n = F[n,:]
        f_nm1 = F[n-1,:]
        y = copy(Y[n,:])

        #implements Gauss-Seidel relaxation
        for _ in 1:sweeps
            y_old = copy(y)
            f_next = f(t_next, y)
            y = Y[n,:] .+ h .* ((5/12)*f_next .+ (2/3)*f_n .- (1/12)*f_nm1)
            if norm(y - y_old) < tol
                break
            end
        end

        Y[n+1,:] = y
        F[n+1,:] = f_next
    end

    return t_grid, Y
end

#solves the nonlinear system of equations with a Gauss-Seidel relaxation AM4
function solve_am4(f, t_span, y0, h; sweeps=12, tol=1e-10)
    #uses AM3 to bootstrap 2 steps
    _, Y = solve_am3(f, (t_span[1], t_span[1]+3*h), y0, h; sweeps=sweeps, tol=tol)
    t0, tf = t_span
    N = ceil(Int, (tf-t0)/h)
    t_grid = range(t0, length=N+1, step=h)

    #compute initial F
    F = zeros(length(t_grid), length(y0))
    for i in 1:4
        F[i,:] = f(t_grid[i], Y[i,:])
    end

    #defines the main AM4 solver
    for n in 4:N
        t_next = t_grid[n+1]
        f_n = F[n,:]
        f_nm1 = F[n-1,:]
        f_nm2 = F[n-2,:]
        y = copy(Y[n,:])

        #implements Gauss-Seidel relaxation
        for _ in 1:sweeps
            y_old = copy(y)
            f_next = f(t_next, y)
            y = Y[n,:] .+ h .* ((3/8)*f_next .+ (19/24)*f_n .- (5/24)*f_nm1 .+ (1/24)*f_nm2)
            if norm(y - y_old) < tol
                break
            end
        end

        Y[n+1,:] = y
        F[n+1,:] = f_next
    end

    return t_grid, Y
end

#solves the nonlinear system of equations with a Gauss-Seidel relaxation AM5
function solve_am5(f, t_span, y0, h; sweeps=12, tol=1e-10)
    #uses AM4 to bootstrap 3 steps of history
    _, Y = solve_am4(f, (t_span[1], t_span[1] + 4*h), y0, h; sweeps=sweeps, tol=tol)
    t0, tf = t_span
    N = ceil(Int, (tf - t0) / h)
    t_grid = range(t0, length=N+1, step=h)

    #computes initial F
    F = zeros(length(t_grid), length(y0))
    for i in 1:5
        F[i,:] = f(t_grid[i], Y[i,:])
    end

    #defines the main AM5 solver
    for n in 5:N
        t_next = t_grid[n+1]
        f_n = F[n,:]
        f_nm1 = F[n-1,:]
        f_nm2 = F[n-2,:]
        f_nm3 = F[n-3,:]
        y = copy(Y[n,:])

        #implements Gauss-Seidel relaxation
        for _ in 1:sweeps
            y_old = copy(y)
            f_next = f(t_next, y)
            y = Y[n,:] .+ h .* (
                (251/720)*f_next .+
                (646/720)*f_n .-
                (264/720)*f_nm1 .+
                (106/720)*f_nm2 .-
                (19/720)*f_nm3
            )
            if norm(y - y_old) < tol
                break
            end
        end

        Y[n+1,:] = y
        F[n+1,:] = f_next
    end

    return t_grid, Y
end