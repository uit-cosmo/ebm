module North
    using PyCall
    using Jacobi
    using LinearAlgebra
    using DifferentialEquations
    using StatsBase
    
    A = 211.2 - 18.
    B = 1/0.32
    D = 0.38
    s1 = -0.796
    s2 = -0.482
    s22 = 0.147
    Tc = -10
    b0 = 0.38
    a0 = 0.697
    a2 = -0.0779
    Q = 334.4
    c = 13.2
    n = 500 
    Delta_theta = 1/n
    r_E = 6371
    X = (collect(1:1:n ) .-0.5) * Delta_theta
    theta = collect(0:1:n-1) * Delta_theta
    
    lam = (1 .- theta.^2)/(Delta_theta.^2)

    CO2_parameter = 6
    CO2_base = 280
    v = ones(n)
    T0 = 30 .- ((collect(1:1:n) .- 0.5) .* 1/n) .* 50


    function construct_laplacian()
        py"""
        import numpy as np 
        from scipy import sparse

        def construct_laplacian_py(n, lam): 
            #Laplacian Construction.
            
            # Representation of sparse matrix and right-hand side
            main  = np.zeros(n)
            lower = np.zeros(n)
            upper = np.zeros(n)

            # Precompute sparse matrix
            main[:] = -np.roll(lam, 1) - lam
            lower[:] = lam  
            upper[:] = lam  

            laplacian = sparse.diags(
                diagonals=[main, lower, upper],
                offsets=[0, -1, 1], shape=(n, n),
                format='csr')
            laplacian.toarray()
            
            # boundary conditions
            laplacian[0, 0] = -lam[0]
            laplacian[0, 1] = lam[0]
            laplacian[-1, -2] = lam[-2]
            laplacian[-1, -1] = -lam[-2]

            return laplacian.toarray()
        """
        x = py"construct_laplacian_py"(n, lam)
    end

    function US(x)
        (1 .+ tanh.(x))/2
    end

    function albedo(T)
        albedo = b0*US(Tc.-T) + (a0 .+ a2 .* legendre.(X, 2)) .* US(T.-Tc)
    end

    function solar_incidation(t)
        1 .+ s1 * cos(2 * pi * t) * legendre.(X, 1) + 
            (s2 + s22 + cos(4 * pi * t)) * legendre.(X, 2)
    end

    function forcing_CO2(t)
        CO2ppm = 280 
        CO2 = CO2_parameter * log((CO2ppm*1.02^t)/CO2_base)
    end

    function define_equation(dTdt, T, p, t) 
        seasonality, albedo_td = p 
        if albedo_td == true
            a = albedo(T)
        else
            a = albedo(T0) 
        end
        if seasonality == true
            S = solar_incidation(t)
        else 
            S = 1 .+ s2 * legendre.(X, 2)
        end

        dTdt .= (-A*v - B*T + D * (laplacian * T) + Q*S.*a + forcing_CO2(t)*v)/c
    end

    function g(dTdt, T,p,t)
        dTdt .= 0.5
    end

    function calculate_SIA(sol)
        X = sol
        #ind = mapslices(t -> replace(u -> findfirst(x -> x < -10, t), Nothing => 0), X, dims = 1)
        ind = Array{Int16}(undef, size(X)[2]) #zeros(size(X)[2])
        for i = 1:size(X)[2]
            if isnothing(findfirst(x -> x <- 10, X[:, i])) 
                ind[i] = n
            else 
                ind[i] = findfirst(x -> x <- 10, X[:, i])
            end
           
        end
        area = 2 .* pi .* r_E^2 .* (1 .- theta[ind])
    end

    function calculate_global_average_T(sol) 
        weight_area = 2 .* pi .* r_E^2 .* (1 .- theta)
        global_temp = mean(sol, weights(weight_area), dims = 1)
    end



    function initialize_problem(p, noise, tspan)
        if noise == true 
            W = WienerProcess(0.0,0.0,0.0)
            prob = SDEProblem(define_equation, g, T0, tspan, p, noise = W)
        else 
            prob = ODEProblem{isinplace}(define_equation, T0, tspan, p)   
        end
    end
    const laplacian = construct_laplacian()

end