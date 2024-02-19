module raices

export biseccion,falsaposicion

    using ForwardDiff

    function biseccion(f, intervalo, err)
        a, b = intervalo
        iteraciones = 0
        
        if sign(f(a)) == sign(f(b))
            error("F(a) y F(b) no tienen signos opuestos")
        end
        
        while (b - a) / 2 > err
            c = (a + b) / 2
            if sign(f(a)) != sign(f(c))
                b = c
            elseif sign(f(b)) != sign(f(c))
                a = c
            end
            iteraciones += 1
        end

        return Dict("raiz" => (a + b) / 2, "error" => (b - a) / 2, "iteraciones" => iteraciones)
    end

    function falsaposicion(f, intervalo, error)
        a, b = intervalo
        iteraciones = 0  
        if sign(f(a)) == sign(f(b))
            error("F(a) y F(b) no tienen signos opuestos")
        end
        c = a 
        c_previo = a 
        error_relativo = Inf 
    
        while true
            c_previo = c
            c = b - f(b) * (b - a) / (f(b) - f(a))
            if iteraciones > 0
                error_relativo = abs(c - c_previo) / abs(c)
            end
            if  error_relativo < error
                break
            end
            if sign(f(a)) == sign(f(c))
                a = c
            else
                b = c
            end
            iteraciones += 1
        end
        return Dict("raiz" => c, "error" => error_relativo, "iteraciones" => iteraciones)
    end

    function puntofijo(g, x0, error)
        max_iter = 100  # Maximum number of iterations
        iteraciones = 0  # Iteration counter
        x = x0  # Initial guess
        x_nuevo = x  # Initialize x_nuevo with the value of x to ensure it's defined
    
        while true
            x_nuevo = g(x)  # Apply the transformation g to the current estimate x
            if abs(x_nuevo - x) < error  # Check if the difference is within the tolerance
                break  # Convergence criterion met
            end
            x = x_nuevo  # Update x to the new estimate for the next iteration
            iteraciones += 1
            if iteraciones >= max_iter  # Check if maximum iterations reached
                print("Maximum iterations exceeded")
                return Dict("raiz" => "No encontrada", "error" => "N/A", "iteraciones" => iteraciones)
            end
        end
        return Dict("raiz" => x, "error" => abs(x_nuevo - x), "iteraciones" => iteraciones)
    end
    
end 