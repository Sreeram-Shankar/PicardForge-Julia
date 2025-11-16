#solves the nonlinear system of equations with a Gauss-Seidel relaxation BE (BDF1)
function solve_be(f, t_span, y0, h; sweeps=12, tol=1e-10)
    t0, tf = t_span
    N = ceil(Int, (tf - t0)/h)
    t_grid = range(t0, length=N+1, step=h)
    Y = zeros(length(t_grid), length(y0))
    Y[1,:] = y0

    #defines the main BE solver
    for n in 1:N
        t = t_grid[n]
        t_next = t + h
        y = copy(Y[n,:])

        #implements Gauss-Seidel relaxation
        for _ in 1:sweeps
            y_old = copy(y)
            y = Y[n,:] .+ h .* f(t_next, y)
            if norm(y - y_old) < tol
                break
            end
        end

        Y[n+1,:] = y
    end

    return t_grid, Y
end


#solves the nonlinear system of equations with a Gauss-Seidel relaxation BDF2
function solve_bdf2(f, t_span, y0, h; sweeps=12, tol=1e-10)
    #uses backward euler to bootstrap
    _, Y_BE = solve_be(f, (t_span[1], t_span[1]+h), y0, h; sweeps=sweeps, tol=tol)
    y1 = Y_BE[end,:]

    #defines the main BDF2 solver
    t0, tf = t_span
    N = ceil(Int, (tf - t0)/h)
    t_grid = range(t0, length=N+1, step=h)
    Y = zeros(length(t_grid), length(y0))
    Y[1,:] = y0
    Y[2,:] = y1

    for n in 2:N
        y_prev = Y[n-1,:]
        y = copy(Y[n,:])
        t_next = t_grid[n+1]

        #defines the rhs
        rhs = (-4*Y[n,:] .+ y_prev) / (2*h)

        #implements Gauss-Seidel relaxation
        for _ in 1:sweeps
            y_old = copy(y)
            f_val = f(t_next, y)
            y = (rhs .+ f_val) / (3/(2*h))
            if norm(y - y_old) < tol
                break
            end
        end

        Y[n+1,:] = y
    end

    return t_grid, Y
end


#solves the nonlinear system of equations with a Gauss-Seidel relaxation BDF3
function solve_bdf3(f, t_span, y0, h; sweeps=12, tol=1e-10)
    #bootstrap using BE then BDF2
    _, Y_BE = solve_be(f, (t_span[1], t_span[1]+h), y0, h; sweeps=sweeps, tol=tol)
    y1 = Y_BE[end,:]
    _, Y_BDF2 = solve_bdf2(f, (t_span[1], t_span[1]+2*h), y0, h; sweeps=sweeps, tol=tol)
    y2 = Y_BDF2[end,:]

    #defines the main BDF3 solver
    t0, tf = t_span
    N = ceil(Int, (tf-t0)/h)
    t_grid = range(t0, length=N+1, step=h)

    Y = zeros(length(t_grid), length(y0))
    Y[1,:] = y0
    Y[2,:] = y1
    Y[3,:] = y2

    for n in 3:N
        y_prev3 = Y[n-3,:]
        y_prev2 = Y[n-2,:]
        y_prev1 = Y[n-1,:]
        y = copy(Y[n,:])
        t_next = t_grid[n+1]

        #defines the rhs
        rhs = (-11*Y[n,:] .+ 18*y_prev1 .- 9*y_prev2 .+ 2*y_prev3)/(6*h)

        #implements Gauss-Seidel relaxation
        for _ in 1:sweeps
            y_old = copy(y)
            f_val = f(t_next, y)
            y = (rhs .+ f_val) / (11/(6*h))
            if norm(y - y_old) < tol
                break
            end
        end

        Y[n+1,:] = y
    end

    return t_grid, Y
end


#solves the nonlinear system of equations with a Gauss-Seidel relaxation BDF4
function solve_bdf4(f, t_span, y0, h; sweeps=12, tol=1e-10)
    #bootstrap using BDF1–3
    _, Y_BE = solve_be(f, (t_span[1], t_span[1]+h), y0, h; sweeps=sweeps, tol=tol)
    y1 = Y_BE[end,:]
    _, Y_BDF2 = solve_bdf2(f, (t_span[1], t_span[1]+2*h), y0, h; sweeps=sweeps, tol=tol)
    y2 = Y_BDF2[end,:]
    _, Y_BDF3 = solve_bdf3(f, (t_span[1], t_span[1]+3*h), y0, h; sweeps=sweeps, tol=tol)
    y3 = Y_BDF3[end,:]

    #defines the main BDF4 solver
    t0, tf = t_span
    N = ceil(Int, (tf-t0)/h)
    t_grid = range(t0, length=N+1, step=h)

    Y = zeros(length(t_grid), length(y0))
    Y[1,:] = y0
    Y[2,:] = y1
    Y[3,:] = y2
    Y[4,:] = y3

    for n in 4:N
        y_prev4 = Y[n-4,:]
        y_prev3 = Y[n-3,:]
        y_prev2 = Y[n-2,:]
        y_prev1 = Y[n-1,:]
        y = copy(Y[n,:])
        t_next = t_grid[n+1]

        #defines the rhs
        rhs = (-25*Y[n,:] .+ 48*y_prev1 .- 36*y_prev2 .+ 16*y_prev3 .- 3*y_prev4)/(12*h)

        #implements Gauss-Seidel relaxation
        for _ in 1:sweeps
            y_old = copy(y)
            f_val = f(t_next, y)
            y = (rhs .+ f_val) / (25/(12*h))
            if norm(y - y_old) < tol
                break
            end
        end

        Y[n+1,:] = y
    end

    return t_grid, Y
end


#solves the nonlinear system of equations with a Gauss-Seidel relaxation BDF5
function solve_bdf5(f, t_span, y0, h; sweeps=12, tol=1e-10)
    #bootstrap up to BDF4
    _, Y_BE = solve_be(f, (t_span[1], t_span[1]+h), y0, h; sweeps=sweeps, tol=tol)
    y1 = Y_BE[end,:]
    _, Y_BDF2 = solve_bdf2(f, (t_span[1], t_span[1]+2*h), y0, h; sweeps=sweeps, tol=tol)
    y2 = Y_BDF2[end,:]
    _, Y_BDF3 = solve_bdf3(f, (t_span[1], t_span[1]+3*h), y0, h; sweeps=sweeps, tol=tol)
    y3 = Y_BDF3[end,:]
    _, Y_BDF4 = solve_bdf4(f, (t_span[1], t_span[1]+4*h), y0, h; sweeps=sweeps, tol=tol)
    y4 = Y_BDF4[end,:]

    #defines the main BDF5 solver
    t0, tf = t_span
    N = ceil(Int, (tf-t0)/h)
    t_grid = range(t0, length=N+1, step=h)

    Y = zeros(length(t_grid), length(y0))
    Y[1,:] = y0
    Y[2,:] = y1
    Y[3,:] = y2
    Y[4,:] = y3
    Y[5,:] = y4

    for n in 5:N
        y_prev5 = Y[n-5,:]
        y_prev4 = Y[n-4,:]
        y_prev3 = Y[n-3,:]
        y_prev2 = Y[n-2,:]
        y_prev1 = Y[n-1,:]
        y = copy(Y[n,:])
        t_next = t_grid[n+1]

        #defines the rhs
        rhs = (-137*Y[n,:] .+ 300*y_prev1 .- 300*y_prev2 .+ 200*y_prev3 .- 75*y_prev4 .+ 12*y_prev5)/(60*h)

        #implements Gauss-Seidel relaxation
        for _ in 1:sweeps
            y_old = copy(y)
            f_val = f(t_next, y)
            y = (rhs .+ f_val) / (137/(60*h))
            if norm(y - y_old) < tol
                break
            end
        end

        Y[n+1,:] = y
    end

    return t_grid, Y
end


#solves the nonlinear system of equations with a Gauss-Seidel relaxation BDF6
function solve_bdf6(f, t_span, y0, h; sweeps=12, tol=1e-10)
    #bootstrap up to BDF5
    _, Y_BE = solve_be(f, (t_span[1], t_span[1]+h), y0, h; sweeps=sweeps, tol=tol)
    y1 = Y_BE[end,:]
    _, Y_BDF2 = solve_bdf2(f, (t_span[1], t_span[1]+2*h), y0, h; sweeps=sweeps, tol=tol)
    y2 = Y_BDF2[end,:]
    _, Y_BDF3 = solve_bdf3(f, (t_span[1], t_span[1]+3*h), y0, h; sweeps=sweeps, tol=tol)
    y3 = Y_BDF3[end,:]
    _, Y_BDF4 = solve_bdf4(f, (t_span[1], t_span[1]+4*h), y0, h; sweeps=sweeps, tol=tol)
    y4 = Y_BDF4[end,:]
    _, Y_BDF5 = solve_bdf5(f, (t_span[1], t_span[1]+5*h), y0, h; sweeps=sweeps, tol=tol)
    y5 = Y_BDF5[end,:]

    #defines the main BDF6 solver
    t0, tf = t_span
    N = ceil(Int, (tf-t0)/h)
    t_grid = range(t0, length=N+1, step=h)

    Y = zeros(length(t_grid), length(y0))
    Y[1,:] = y0
    Y[2,:] = y1
    Y[3,:] = y2
    Y[4,:] = y3
    Y[5,:] = y4
    Y[6,:] = y5

    for n in 6:N
        y_prev6 = Y[n-6,:]
        y_prev5 = Y[n-5,:]
        y_prev4 = Y[n-4,:]
        y_prev3 = Y[n-3,:]
        y_prev2 = Y[n-2,:]
        y_prev1 = Y[n-1,:]
        y = copy(Y[n,:])
        t_next = t_grid[n+1]

        #defines the rhs
        rhs = (-147*Y[n,:] .+ 360*y_prev1 .- 450*y_prev2 .+ 400*y_prev3 .- 225*y_prev4 .+ 72*y_prev5 .- 10*y_prev6)/(60*h)

        #implements Gauss-Seidel relaxation
        for _ in 1:sweeps
            y_old = copy(y)
            f_val = f(t_next, y)
            y = (rhs .+ f_val) / (147/(60*h))
            if norm(y - y_old) < tol
                break
            end
        end

        Y[n+1,:] = y
    end

    return t_grid, Y
end