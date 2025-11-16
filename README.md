# 📦 PicardForge-Julia
### *A fixed-step ODE & PDE time-integration library featuring Picard–Gauss–Seidel implicit solvers and classical explicit methods — now in native Julia*

---

## ✨ Overview

**PicardForge-Julia** is a lightweight numerical solver library implementing a complete suite of **fixed-step ODE solvers**, designed especially for **semi-discretized PDEs** (heat equation, diffusion, conduction models, etc.).

Unlike traditional implicit solvers that rely on Newton iterations or Jacobian solves, **PicardForge-Julia uses Picard fixed-point iteration with Gauss–Seidel relaxation** to solve all implicit systems. This gives very efficient and stable integration for diffusion-type PDEs while maintaining full A-stability (and L-stability for Radau).

The library includes:

- **Explicit Runge–Kutta (RK1–RK6)**
- **Adams–Bashforth multistep (AB2–AB5)**
- **Adams–Moulton multistep (AM2–AM5)**
- **Backward Differentiation Formulas (BDF1–BDF6)**
- **SDIRK (2nd–4th order)**
- **Fully implicit Gauss, Radau, and Lobatto collocation IRK (s = 1–5)**

All implicit families are solved using **Picard–Gauss–Seidel**, giving a robust, Jacobian-free integration strategy ideal for large systems from PDE discretization.

---

## 🚀 Features

### ✔ Full suite of classic numerical integrators

| Family | Methods | Notes |
|-------|---------|-------|
| **Explicit RK** | RK1–RK6 | Classical explicit Butcher tables |
| **Adams–Bashforth** | AB2–AB5 | Explicit multistep |
| **Adams–Moulton** | AM2–AM5 | Implicit multistep (Picard solved) |
| **BDF** | BDF1–BDF6 | All implicit, Picard–GS iteration |
| **SDIRK** | SDIRK2–SDIRK4 | Diagonally implicit RK |
| **Gauss–Legendre IRK** | s = 1–5 | A-stable, symplectic |
| **Radau IIA IRK** | s = 2–5 | L-stable, great for stiff PDEs |
| **Lobatto IIIC IRK** | s = 2–5 | Symmetric, stiffly accurate |

### ✔ Picard–Gauss–Seidel nonlinear iteration

- No Jacobian matrices  
- No Newton factorization  
- Stage-by-stage relaxation  
- Matrix-free  
- Very effective for diffusion-dominated PDEs  

### ✔ Static Butcher tableaus included

All Gauss, Radau, and Lobatto tables (s = 1–5) are embedded directly in `irk.jl`, generated using high precision (from MethodForge-style quadrature construction) and stored as `const` arrays.

### ✔ PDE-ready design

Well-suited for time integration of systems obtained from:

- radial diffusion / conduction  
- thermal evolution models  
- finite-difference / finite-volume parabolic PDEs  
- generic semi-discretized systems \( u_t = F(u) \)

---

## 📁 Repository Structure

PicardForge-Julia/
│
├── rk.jl # RK1–RK6 explicit solvers
├── ab.jl # AB2–AB5 explicit multistep
├── am.jl # AM2–AM5 implicit multistep (Picard–GS)
├── bdf.jl # BDF1–BDF6 implicit multistep (Picard–GS)
├── irk.jl # Gauss/Radau/Lobatto collocation IRK (Picard–GS)
├── sdirk.jl # SDIRK2–SDIRK4 implicit RK
└── PicardForge.jl # Unified module export

markdown
Copy code

Each solver file contains a family of functions of the form:

- `solve_rk1`, `solve_rk2`, …, `solve_rk6`
- `solve_ab2`, …, `solve_ab5`
- `solve_am2`, …, `solve_am5`
- `solve_bdf1` (or `solve_be`), …, `solve_bdf6`
- `solve_sdirk2`, `solve_sdirk3`, `solve_sdirk4`
- `solve_collocation` (for Gauss/Radau/Lobatto IRK)

and these are re-exported via `PicardForge.jl`.

---

## 🧠 How Picard–Gauss–Seidel nonlinear iteration works

All implicit collocation and multistep methods solve fixed-point equations of the form

\[
Y = y_n + h A F(Y)
\]

for stages \( Y = (Y_1, \dots, Y_s) \), where \( A \) is the Butcher matrix and \( F(Y) = (f(t_n + c_i h, Y_i))_{i=1}^s \).

Instead of Newton’s method (which forms Jacobians and solves linear systems), we apply a **Picard fixed-point iteration**:

\[
Y^{(k+1)} = y_n + h A F(Y^{(k)})
\]

This is implemented using **Gauss–Seidel stage relaxation** in Julia:

```julia
function step_collocation(f, t, y, h, A, b, c; sweeps=12, tol=1e-10)
    s = length(b)
    n = length(y)

    # initial stage guesses
    Y = [copy(y) for _ in 1:s]

    for _ in 1:sweeps
        Y_old = deepcopy(Y)

        for i in 1:s
            rhs = zeros(eltype(y), n)
            for j in 1:s
                rhs .+= A[i,j] * f(t + c[j]*h, Y[j])
            end
            Y[i] = y .+ h .* rhs
        end

        # global stage-difference norm
        if norm(vcat([Y[i] .- Y_old[i] for i in 1:s]...)) < tol
            break
        end
    end

    K = [f(t + c[i]*h, Y[i]) for i in 1:s]
    y_next = y .+ h .* sum(b[i] .* K[i] for i in 1:s)
    return y_next
end
