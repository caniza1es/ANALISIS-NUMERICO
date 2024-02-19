include("raices.jl")

# Define the function whose root is to be found
f(x) = x^3 - 2x - 5

# Define the interval and error frontier
intervalo = (2.0, 3.0)
error = 1e-5

# Call the metodo_biseccion function
resultado = raices.biseccion(f, intervalo, error)

# Print the result
println("Root: ", resultado["raiz"])
println("Error: ", resultado["error"])
println("Iterations: ", resultado["iteraciones"])
println("valor evaluado ",f(resultado["raiz"]))

resultado = raices.falsaposicion(f, intervalo, error)
println("----------------------")
# Print the result
println("Root: ", resultado["raiz"])
println("Error: ", resultado["error"])
println("Iterations: ", resultado["iteraciones"])
println("valor evaluado ",f(resultado["raiz"]))

# Print the result
g(x) = 0.5 * (x + 2 / x)
resultado = raices.puntofijo(g, 1, error)


println("----------------------")
println("Root: ", resultado["raiz"])
println("Error: ", resultado["error"])
println("Iterations: ", resultado["iteraciones"])
println("valor evaluado ",g(resultado["raiz"]))