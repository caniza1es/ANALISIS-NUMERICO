#DERIVACION

function derivada_tres_puntos(f, x, h)
    (f(x + h) - f(x - h)) / (2h)
end

function derivada_cinco_puntos(f, x, h)
    (-f(x + 2h) + 8f(x + h) - 8f(x - h) + f(x - 2h)) / (12h)
end

function extrapolacion_richardson(f, x, h1, h2)
    D1 = derivada_tres_puntos(f, x, h1)
    D2 = derivada_tres_puntos(f, x, h2)
    D2 + (D2 - D1) / ((h1/h2)^2 - 1)
end

function segunda_derivada_cinco_puntos(f, x, h)
    return (-f(x + 2*h) + 16*f(x + h) - 30*f(x) + 16*f(x - h) - f(x - 2*h)) / (12*h^2)
end

function puntos_inflexion(f, a, b, h)
    x = a
    puntos_inflexion = []
    
    prev_fpp = segunda_derivada_cinco_puntos(f, x, h)
    x += h

    while x <= b
        fpp = segunda_derivada_cinco_puntos(f, x, h)
        if prev_fpp * fpp < 0
            push!(puntos_inflexion, x - h/2) 
        end
        
        x += h
        prev_fpp = fpp
    end
    
    return puntos_inflexion
end

#INTEGRACION

function trapecio_simple(a, b, f)
    return (f(a) + f(b)) * (b - a) / 2
end

function trapecio_compuesto(a, b, f, n)
    h = (b - a) / n
    suma = f(a) + f(b)
    for i = 1:(n-1)
        suma += 2 * f(a + i*h)
    end
    return suma * h / 2
end

function simpson_un_tercio(a, b, f)
    h = (b - a) / 2
    return (f(a) + 4*f(a + h) + f(b)) * h / 3
end

function simpson_un_tercio_compuesto(a, b, f, n)
    h = (b - a) / n
    suma = f(a) + f(b)
    for i = 1:(n-1)
        coeficiente = 3 - (i % 2) * 2  
        suma += coeficiente * f(a + i*h)
    end
    return suma * h / 3
end

function simpson_tres_octavos(a, b, f)
    h = (b - a) / 3
    return (f(a) + 3*f(a + h) + 3*f(a + 2*h) + f(b)) * 3 * h / 8
end

function simpson_tres_octavos_compuesto(a, b, f, n)
    #n multiplo de 3
    h = (b - a) / n
    suma = f(a) + f(b)
    for i = 1:(n-1)
        coeficiente = if i % 3 == 0; 2; else; 3; end
        suma += coeficiente * f(a + i*h)
    end
    return suma * 3 * h / 8
end

function error_relativo(exacto, aproximado)
    return abs((exacto - aproximado) / exacto) * 100
end

function romberg_integracion(f, a, b, tol)
    R = zeros(Float64, 10, 10)  
    h = b - a
    R[1, 1] = (f(a) + f(b)) * h / 2  
    n = 1
    for i in 2:10
        h /= 2
        sum = 0.0
        for k in 1:n
            sum += f(a + (2k-1) * h)
        end
        R[i, 1] = 0.5 * R[i-1, 1] + h * sum  
        for j in 2:i
            R[i, j] = (4^(j-1) * R[i, j-1] - R[i-1, j-1]) / (4^(j-1) - 1)  
        end
        if i > 2 && abs(R[i, i] - R[i-1, i-1]) < tol
            return R[i, i], i  
        end
        n *= 2
    end
    error("max iteraciones")
end

#VALOR INICIAL


function euler(f, x0, y0, h, n)
    x = x0
    y = y0
    for i = 1:n
        y += h * f(x, y)
        x += h
    end
    return y
end

function midpoint(f, x0, y0, h, n)
    x = x0
    y = y0
    for i = 1:n
        k1 = f(x, y)
        k2 = f(x + h/2, y + h/2 * k1)
        y += h * k2
        x += h
    end
    return y
end

function modified_euler(f, x0, y0, h, n)
    x = x0
    y = y0
    for i = 1:n
        k1 = f(x, y)
        k2 = f(x + h, y + h * k1)
        y += h * (k1 + k2) / 2
        x += h
    end
    return y
end

function runge_kutta(f, x0, y0, h, n)
    x = x0
    y = y0
    for i = 1:n
        k1 = f(x, y)
        k2 = f(x + h/2, y + h/2 * k1)
        k3 = f(x + h/2, y + h/2 * k2)
        k4 = f(x + h, y + h * k3)
        y += h * (k1 + 2*k2 + 2*k3 + k4) / 6
        x += h
    end
    return y
end

function runge_kutta__4(y0, t0, tf, h)
    t = t0
    y = y0
    while y > 0
        k1 = h * f(t, y)
        k2 = h * f(t + 0.5*h, y + 0.5*k1)
        k3 = h * f(t + 0.5*h, y + 0.5*k2)
        k4 = h * f(t + h, y + k3)
        y += (k1 + 2*k2 + 2*k3 + k4) / 6
        t += h
        if y < 0.0001  #tanque vaci
            break
        end
    end
    return t
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

function eulerk(f, p0, t0, tf, h)
    times = t0:h:tf
    p = zeros(length(times))
    p[1] = p0
    for i in 2:length(times)
        p[i] = p[i-1] + h * f(p[i-1])
    end
    return times, p
end

function rk4(f, p0, t0, tf, h)
    times = t0:h:tf
    p = zeros(length(times))
    p[1] = p0
    for i in 1:(length(times)-1)
        t = times[i]
        p_i = p[i]
        k1 = h * f(p_i)
        k2 = h * f(p_i + 0.5*k1)
        k3 = h * f(p_i + 0.5*k2)
        k4 = h * f(p_i + k3)
        p[i+1] = p_i + (k1 + 2*k2 + 2*k3 + k4) / 6
    end
    return times, p
end
