module numerico
    export tres_puntos,cinco_puntos,regla_trapezoidal

    function tres_puntos(f, x, h)
        D1 = (f(x - h) - f(x + h)) / (2h)
        D2 = (f(x - h/2) - f(x + h/2)) / h
        D = D2 + (D2 - D1) / 3
        return D
    end

    function cinco_puntos(f, x, h)
        D1 = (-f(x + 2h) + 8f(x + h) - 8f(x - h) + f(x - 2h)) / (12h)
        D2 = (-f(x + h) + 8f(x + h/2) - 8f(x - h/2) + f(x - h)) / (6h)
        D = D2 + (D2 - D1) / 15
        return D
    end
    
    function regla_trapezoidal(f, a, b)
        h = b - a
        return (h / 2) * (f(a) + f(b))
    end

    function regla_de_simpson(f, a, b)
        h = (b - a) / 2
        x0 = a
        x1 = a + h
        x2 = b
        return (h / 3) * (f(x0) + 4*f(x1) + f(x2))
    end

    function tres_octavos(f, a, b)
        h = (b - a) / 3
        return (3h / 8) * (f(a) + 3*f(a + h) + 3*f(a + 2h) + f(b))
    end
    
    function regla_de_boole(f, a, b)
        h = (b - a) / 4
        return (2h / 45) * (7f(a) + 32f(a + h) + 12f(a + 2h) + 32f(a + 3h) + 7f(b))
    end
    
    

end