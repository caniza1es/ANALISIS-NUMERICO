function gauss_elimination(A, b)
    n = length(b)
    Ab = [A b]
    for i in 1:n-1
        for j in i+1:n
            if Ab[i, i] == 0
                error("0 pivot")
            end
            m = Ab[j, i] / Ab[i, i]
            Ab[j, :] = Ab[j, :] - m * Ab[i, :]
        end
    end
    
    return Ab
end

function back_substitution(Ab)
    n = size(Ab, 1)
    x = zeros(n)
    
    for i in n:-1:1
        x[i] = (Ab[i, end] - sum(Ab[i, i+1:end-1] .* x[i+1:end])) / Ab[i, i]
    end
    
    return x
end

function solve_gaussian(A, b)
    Ab = gauss_elimination(A, b)
    return back_substitution(Ab)
end

function partial_pivoting(A, b)
    n = length(b)
    Ab = [A b]
    
    for i in 1:n-1

        pivot_row = i + argmax(abs.(Ab[i:n, i])) - 1

        Ab[i, :], Ab[pivot_row, :] = Ab[pivot_row, :], Ab[i, :]
        for j in i+1:n
            if Ab[i, i] == 0
                error("0 pivot")
            end
            m = Ab[j, i] / Ab[i, i]
            Ab[j, :] = Ab[j, :] - m * Ab[i, :]
        end
    end
    
    return Ab
end

function solve_partial_pivoting(A, b)
    Ab = partial_pivoting(A, b)
    return back_substitution(Ab)
end

function scaled_pivoting(A, b)
    n = length(b)
    Ab = [A b]
    s = maximum(abs.(A), dims=2) 
    s = s[:, 1]  
    for i in 1:n-1
        pivot_row = i + argmax(abs.(Ab[i:n, i]) ./ s[i:n]) - 1
        Ab[i, :], Ab[pivot_row, :] = Ab[pivot_row, :], Ab[i, :]
        s[i], s[pivot_row] = s[pivot_row], s[i]
        for j in i+1:n
            if Ab[i, i] == 0
                error("0 pivot")
            end
            m = Ab[j, i] / Ab[i, i]
            Ab[j, :] = Ab[j, :] - m * Ab[i, :]
        end
    end
    
    return Ab
end

function solve_scaled_pivoting(A, b)
    Ab = scaled_pivoting(A, b)
    return back_substitution(Ab)
end

function lu_factorization(A)
    n = size(A, 1)
    L = Matrix{Float64}(I(n))  
    U = copy(A)  

    for i in 1:n-1
        for j in i+1:n
            if U[i, i] == 0
                error("0 pivot")
            end
            L[j, i] = U[j, i] / U[i, i]
            U[j, :] = U[j, :] - L[j, i] * U[i, :]
        end
    end

    return L, U
end

function I(n)
    return Diagonal(ones(n))
end

function Diagonal(v)
    n = length(v)
    D = zeros(n, n)
    for i in 1:n
        D[i, i] = v[i]
    end
    return D
end


function forward_substitution(L, b)
    n = length(b)
    y = zeros(n)

    for i in 1:n
        y[i] = (b[i] - sum(L[i, 1:i-1] .* y[1:i-1])) / L[i, i]
    end

    return y
end

function inverse_using_lu(A)
    n = size(A, 1)
    L, U = lu_factorization(A)
    I = Diagonal(ones(n))  
    invA = zeros(n, n)

    for i in 1:n
        e = I[:, i]
        y = forward_substitution(L, e)
        invA[:, i] = back_substitution([U y])
    end

    return invA
end

function solve_using_inverse(A, b)
    invA = inverse_using_lu(A)
    return invA * b
end




A = [3.0 -0.1 -0.2; 0.1 7.0 -0.3; 0.3 -0.2 10.0]
b = [7.85; -19.3; 71.4]


x_gaussian = solve_gaussian(A, b)
println("Solución usando eliminación gaussiana: ", x_gaussian)


x_partial = solve_partial_pivoting(A, b)
println("Solución usando pivoteo parcial: ", x_partial)


x_scaled = solve_scaled_pivoting(A, b)
println("Solución usando pivoteo escalado: ", x_scaled)



x_inverse = solve_using_inverse(A, b)
println("Solución usando la inversa: ", x_inverse)


#METODOS ITERATIVOS

function jacobi(A, b, x0, tol, max_iter)
    n = length(b)
    x = copy(x0)
    x_new = similar(x)

    iter = 0
    while iter < max_iter
        for i in 1:n
            sigma = 0.0
            for j in 1:n
                if j != i
                    sigma += A[i, j] * x[j]
                end
            end
            x_new[i] = (b[i] - sigma) / A[i, i]
        end

        if norm(x_new - x) < tol
            return x_new, iter
        end

        x .= x_new
        iter += 1
    end

    return x, iter
end

function gauss_seidel(A, b, x0, tol, max_iter)
    n = length(b)
    x = copy(x0)
    x_new = similar(x)

    iter = 0
    while iter < max_iter
        for i in 1:n
            sigma = 0.0
            for j in 1:n
                if j != i
                    sigma += A[i, j] * x_new[j]
                end
            end
            x_new[i] = (b[i] - sigma) / A[i, i]
        end

        if norm(x_new - x) < tol
            return x_new, iter
        end

        x .= x_new
        iter += 1
    end

    return x, iter
end

function iterative_refinement(A, b, x0, tol, max_iter)
    x, _ = jacobi(A, b, x0, tol, max_iter)
    r = b - A * x
    Ar, _ = jacobi(A, r, zeros(length(b)), tol, max_iter)
    x += Ar
    return x
end

function norm(x)
    return sqrt(sum(abs2, x))
end


A = [10.0 2.0 1.0; 1.0 5.0 1.0; 2.0 3.0 10.0]
b = [7.0, -8.0, 6.0]
x0 = [0.0, 0.0, 0.0]
tol = 1e-6
max_iter = 1000


x_jacobi, iter_jacobi = jacobi(A, b, x0, tol, max_iter)
println("Solución usando Jacobi después de $iter_jacobi iteraciones: ", x_jacobi)


x_gs, iter_gs = gauss_seidel(A, b, x0, tol, max_iter)
println("Solución usando Gauss-Seidel después de $iter_gs iteraciones: ", x_gs)


x_refined = iterative_refinement(A, b, x0, tol, max_iter)
println("Solución usando Refinamiento Iterativo: ", x_refined)



function linear_regression(x, y)
    n = length(x)
    X = hcat(ones(n), x)
    β = X \ y
    return β
end


function saturation_growth_regression(x, y)
    n = length(x)
    X = hcat(ones(n), x, x.^2)
    β = X \ log.(y)
    a, b, c = exp(β[1]), exp(β[2]), exp(β[3])
    return a, b, c
end


function power_regression(x, y)
    n = length(x)
    X = hcat(ones(n), log.(x))
    β = X \ log.(y)
    a, b = exp(β[1]), β[2]
    return a, b
end


function quadratic_regression(x, y)
    n = length(x)
    X = hcat(ones(n), x, x.^2)
    β = X \ y
    return β
end




using Plots


x = 1:10
y = [3, 7, 9, 12, 15, 18, 20, 21, 22, 23]


β_linear = linear_regression(x, y)
y_linear = β_linear[1] .+ β_linear[2] * x


a, b, c = saturation_growth_regression(x, y)
y_saturation_growth = a .* (1 .- exp.(-b .* x)) .+ c



a_power, b_power = power_regression(x, y)
y_power = a_power * x.^b_power


β_quadratic = quadratic_regression(x, y)
β0, β1, β2 = β_quadratic
y_quadratic = β0 .+ β1 * x .+ β2 * x.^2


scatter(x, y, label="Datos", legend=:topright)
plot!(x, y_linear, label="Ajuste Lineal", linestyle=:dash)
plot!(x, y_saturation_growth, label="Ajuste Tasa de Crecimiento de Saturación", linestyle=:dot)
plot!(x, y_power, label="Ajuste Potencial", linestyle=:dashdot)
plot!(x, y_quadratic, label="Ajuste Parabólico", linestyle=:dashdotdot)
xlabel!("X")
ylabel!("Y")
title!("Ajuste de Datos a Diferentes Modelos")



function legendre(n, x)
    if n == 0
        return 1.0
    elseif n == 1
        return x
    else
        return ((2 * n - 1) * x * legendre(n - 1, x) - (n - 1) * legendre(n - 2, x)) / n
    end
end


function legendre_polynomials(n, x)
    P = [legendre(i, xi) for xi in x, i in 0:n-1]
    return P
end


function legendre_least_squares_fit(x, y, n)
    P = legendre_polynomials(n, x)
    coefficients = P \ y
    return coefficients
end


x_test = collect(range(-1, stop=1, length=100))
y_test = exp.(x_test)


coefficients_test = legendre_least_squares_fit(x_test, y_test, 4)


plot(x_test, y_test, label="Exp(x)", xlabel="x", ylabel="f(x)", legend=:topleft)
plot!(x_test, legendre_polynomials(4, x_test) * coefficients_test, label="Legendre Fit")


function fourier_fit(x, y, num_terms)
    N = length(x)
    coeffs = zeros(2*num_terms)
    for i in 1:num_terms
        coeffs[2*i-1] = sum(y .* cos.(i .* x)) / N
        coeffs[2*i] = sum(y .* sin.(i .* x)) / N
    end
    return coeffs
end


function evaluate_fourier(x, coeffs)
    num_terms = length(coeffs) ÷ 2
    y = zeros(length(x))
    for i in 1:num_terms
        y .+= coeffs[2*i-1] * cos.(i .* x) .+ coeffs[2*i] * sin.(i .* x)
    end
    return y
end


x = range(0, stop=2π, length=100)
y_true = sin.(x) + 0.5 .* cos.(2 .* x) + 0.2 .* sin.(3 .* x)
noise = 0.1 .* randn(length(x))
data = y_true .+ noise

num_terms = 5
coeffs = fourier_fit(x, data, num_terms)
fitted_data = evaluate_fourier(x, coeffs)

plot(x, data, label="Data", xlabel="x", ylabel="f(x)", legend=:topleft)
plot!(x, fitted_data, label="Fourier Series Fit")
