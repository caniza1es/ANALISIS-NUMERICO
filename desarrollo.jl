include("main.jl")

f_tan(x) = tan(x / 3)
f_exp_x(x) = exp(x) + x
f_cos(x) = cos(x)


h = 0.1  


derivada_tan_tres_puntos = derivada_tres_puntos(f_tan, 3, h)
derivada_tan_cinco_puntos = derivada_cinco_puntos(f_tan, 3, h)


derivada_exp_x_tres_puntos = derivada_tres_puntos(f_exp_x, 2, h)
derivada_exp_x_cinco_puntos = derivada_cinco_puntos(f_exp_x, 2, h)

h1_richardson = π / 3
h2_richardson = π / 6
derivada_cos_richardson = extrapolacion_richardson(f_cos, π/4, h1_richardson, h2_richardson)

println("\\section{Punto 1}")
println("\\subsection{Punto 1a}")
println("\\textbf{Punto 1a.i:} \\\\")
println("Usando la fórmula de tres puntos, la derivada de \\( y = \\tan\\left(\\frac{x}{3}\\right) \\) en \\( x = 3 \\) es aproximadamente \\( f'(3) \\approx $derivada_tan_tres_puntos \\). \\\\")
println("Usando la fórmula de cinco puntos, la derivada es aproximadamente \\( f'(3) \\approx $derivada_tan_cinco_puntos \\). \\\\")
println("\\textbf{Punto 1a.ii:} \\\\")
println("Usando la fórmula de tres puntos, la derivada de \\( y = e^x + x \\) en \\( x = 2 \\) es aproximadamente \\( f'(2) \\approx $derivada_exp_x_tres_puntos \\). \\\\")
println("Usando la fórmula de cinco puntos, la derivada es aproximadamente \\( f'(2) \\approx $derivada_exp_x_cinco_puntos \\). \\\\")
println("\\subsection{Punto 1b}")
println("La primera derivada de \\( y = \\cos(x) \\) en \\( x = \\frac{\\pi}{4} \\) usando la extrapolación de Richardson con pasos \\( h_1 = \\frac{\\pi}{3} \\) y \\( h_2 = \\frac{\\pi}{6} \\) es aproximadamente \\( f'\\left(\\frac{\\pi}{4}\\right) \\approx $derivada_cos_richardson \\). \\\\")


f_normal(x) = 1 / sqrt(2π) * exp(-x^2 / 2)

a = -10  
b = 10   
h = 0.1  

puntos_inflexion_normal = puntos_inflexion(f_normal, a, b, h)


println("\\subsection{Punto 1c}")
println("Los puntos de inflexión de la función de distribución normal, \\( f(x) = \\frac{1}{\\sqrt{2\\pi}} e^{-\\frac{x^2}{2}} \\), se encuentran aproximadamente en: \\\\")
println("\\[ \\text{Puntos de inflexión:} ")
for punto in puntos_inflexion_normal
    println("$punto, ")
end
println("\\]")


f(x) = 6 + 3 * cos(x)

a = 0
b = π / 2

solucion_analitica = 6 * (b - a) + 3 * sin(b) - 3 * sin(a)

trapecio_simple1 = trapecio_simple(a, b, f)

trapecio_compuesto_n2 = trapecio_compuesto(a, b, f, 2)
trapecio_compuesto_n4 = trapecio_compuesto(a, b, f, 4)

simpson_un_tercio_simple = simpson_un_tercio(a, b, f)

simpson_un_tercio_compuesto_n4 = simpson_un_tercio_compuesto(a, b, f, 4)

simpson_tres_octavos_simple = simpson_tres_octavos(a, b, f)

simpson_tres_octavos_compuesto_n6 = simpson_tres_octavos_compuesto(a, b, f, 6)


println("\\section{Punto 2}")
println("\\subsection{Punto 2a}")
println("\\textbf{Punto 2a.i:} La solución analítica de la integral es \\( $solucion_analitica \\). \\\\")
println("\\textbf{Punto 2a.ii:} Resultado con la regla del trapecio simple: \\( $trapecio_simple1\\). \\\\")
println("\\textbf{Punto 2a.iii:} Resultados con la regla del trapecio compuesto, n=2: \\( $trapecio_compuesto_n2 \\), n=4: \\( $trapecio_compuesto_n4 \\). \\\\")
println("\\textbf{Punto 2a.iv:} Resultado con la regla de Simpson 1/3 simple: \\( $simpson_un_tercio_simple \\). \\\\")
println("\\textbf{Punto 2a.v:} Resultado con la regla de Simpson 1/3 compuesta, n=4: \\( $simpson_un_tercio_compuesto_n4 \\). \\\\")
println("\\textbf{Punto 2a.vi:} Resultado con la regla de Simpson 3/8 simple: \\( $simpson_tres_octavos_simple \\). \\\\")
println("\\textbf{Punto 2a.vii:} Resultado con la regla de Simpson 3/8 compuesta, n=6: \\( $simpson_tres_octavos_compuesto_n6 \\). \\\\")


Q_t(t) = 9 + 5 * cos(0.4 * t)^2
c_t(t) = 5 * exp(-0.5 * t) + 2 * exp(0.15 * t)


integrando(t) = Q_t(t) * c_t(t)


tiempo_inicial = 2
tiempo_final = 8
tolerancia = 0.001  


masa_transportada, iteraciones = romberg_integracion(integrando, tiempo_inicial, tiempo_final, tolerancia)


println("\\section{Punto 2b}")
println("La cantidad de masa transportada por el tubo desde \\( t = $tiempo_inicial \\) minutos hasta \\( t = $tiempo_final \\) minutos se calcula usando la integración de Romberg con una tolerancia de $tolerancia. \\\\")
println("La masa transportada es aproximadamente \\( $masa_transportada \\) mg, calculada en \\( $iteraciones \\) iteraciones de Romberg. \\\\")


g = 9.81  
m = 80    
c = 10    


v_t(t) = (g * m / c) * (1 - exp(-(c / m) * t))


tiempo_inicial = 0
tiempo_final = 8
tolerancia = 0.01  


distancia_recorrida, iteraciones = romberg_integracion(v_t, tiempo_inicial, tiempo_final, tolerancia)


println("\\section{Punto 2c}")
println("La distancia recorrida por el paracaidista en los primeros \\( $tiempo_final \\) segundos se calcula usando la integración de Romberg con una tolerancia del 1%. \\\\")
println("La distancia recorrida es aproximadamente \\( $distancia_recorrida \\) metros, calculada en \\( $iteraciones \\) iteraciones de Romberg. \\\\")


function f(t, y, k=0.06)
    return -k * sqrt(y)
end





y0 = 3.0
t0 = 0.0  
tf = 100.0  
h = 0.25  


tiempo_vaciado = runge_kutta__4(y0, t0, tf, h)
println("El tanque se vacía en aproximadamente $(tiempo_vaciado) minutos.")



f(x, y) = (1 + 4x) * sqrt(y)

# Condiciones iniciales y parámetros de la simulación
x0 = 0.0
y0 = 1.0
h = 0.05
xn = 1.0
n = Int((xn - x0)/h)

# Aplicación de los métodos
euler_y = euler(f, x0, y0, h, n)
midpoint_y = midpoint(f, x0, y0, h, n)
modified_euler_y = modified_euler(f, x0, y0, h, n)
runge_kutta_y = runge_kutta(f, x0, y0, h, n)

print(euler_y," ",midpoint_y," ",modified_euler_y, " " ,runge_kutta_y)

using Plots

m = 20.0
k = 20.0
x0 = [1.0, 0.0] 
t0 = 0.0
tf = 15.0
h = 0.01 
n = Int((tf - t0) / h) 


function sistema_ode(t, y, c)
    x, v = y
    dxdt = v
    dvdt = -(k/m) * x - (c/m) * v
    return [dxdt, dvdt]
end


function runge_kutta_4(f, y0, t0, tf, h, c)
    n = Int((tf - t0) / h)
    t = t0
    y = y0
    ys = [y0[1]]
    ts = [t0]

    for i = 1:n
        k1 = f(t, y, c)
        k2 = f(t + h/2, y + h/2 .* k1, c)
        k3 = f(t + h/2, y + h/2 .* k2, c)
        k4 = f(t + h, y + h .* k3, c)
        y += h * (k1 + 2 .* k2 + 2 .* k3 + k4) / 6
        t += h
        push!(ys, y[1])
        push!(ts, t)
    end
    return ts, ys
end


cs = [5.0, 40.0, 200.0]


plot(title="Desplazamiento del Sistema Masa-Resorte-Amortiguador", xlabel="Tiempo (s)", ylabel="Desplazamiento (m)", legend=:bottomleft)

for c in cs
    t_vals, x_vals = runge_kutta_4((t, y, c) -> sistema_ode(t, y, c), x0, t0, tf, h, c)
    plot!(t_vals, x_vals, label="c = $c Ns/m")
end


savefig("oscilador_amortiguado.png")



k = 0.026
p_max = 12000.0


function modelo_logistico(p)
    return k * (1 - p / p_max) * p
end


data_years = [1950, 1960, 1970, 1980, 1990, 2000]
data_population = [2555, 3040, 3708, 4454, 5276, 6079]


p0 = data_population[1]
t0 = data_years[1]
tf = data_years[end]
h = 1.0 


t_euler, p_euler = eulerk(modelo_logistico, p0, t0, tf, h)


t_rk4, p_rk4 = rk4(modelo_logistico, p0, t0, tf, h)


using Interpolations
interp = LinearInterpolation(data_years, data_population)
t_interp = t0:h:tf
p_interp = interp.(t_interp)


plot(t_euler, p_euler, label="Euler", title="Simulación del Crecimiento de la Población", xlabel="Año", ylabel="Población (millones)")
plot!(t_rk4, p_rk4, label="RK4")
scatter!(data_years, data_population, label="Datos Históricos", color="black")


savefig("poblacion_modelo_logistico.png")


