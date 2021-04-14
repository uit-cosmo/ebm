module North
    #######################################################################
    #This module contains the differential equations and all calculations #
    #including sea ice area calculation, global Temperature calculation...#
    #######################################################################
    using PyCall
    using Jacobi
    using LinearAlgebra
    using DifferentialEquations
    using StatsBase
    using Sundials
    using CSV
    using Tables
    using Interpolations
    using Plots
    using DataFrames
    using JLD2 
    using DelimitedFiles

    runs = 2
    A = 211.2 - 19#18.
    B = 1/0.26
    D = 0.38
    s1 = -0.796
    s2 = -0.482
    s22 = 0.147
    Tc = -10
    b0 = 0.38
    a0 = 0.60
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
    CO2ppm = 280 
    v = ones(n)
    T0 = 30 .- ((collect(1:1:n) .- 0.5) .* 1/n) .* 50
    directory = @__DIR__

    function construct_laplacian()
        """
        Laplacian Construction in Python. 
        """
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
                format="csr")
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
        """
        Function used for albedo calculation.
        """
        (1 .+ tanh.(x))/2
    end

    function albedo(T)
        """
        Function to define the albedo function.
        """
        albedo = b0*US(Tc.-T) + (a0 .+ a2 .* legendre.(X, 2)) .* US(T.-Tc)
    end

    function solar_incidation(t)
        """ 
        Calculates the seasonal incident solar radiation at latitude X, taken from North & Coakley 1979.
        
        Parameters:
        
        t: Time vector in years.    
        """
        1 .+ s1 * cos(2 * pi * t) * legendre.(X, 1) + 
            (s2 + s22 + cos(4 * pi * t)) * legendre.(X, 2) .+ (X .- 0.5).^3/4 #last term is lapse-rate feedback
    end

    function forcing_CO2(t, historical_forcing)
        """
        Formula for CO2 forcing. If historical_forcing == true,
        then the historical forcing is used instead of increasing CO2.

        Parameters:
        
        t: Time vector in years.  
        historical_forcing: boolean 


        Reference:
        Myhre, G., E.J. Highwood, K.P. Shine, and F. Stordal, 1998: 
        New estimates of radiative forcing due to well mixed greenhouse gases. 
        Geophysical Research Letters, 25, 2715-2718.
        """
        if historical_forcing == true 
            f = interpolate_historical_forcing(t) 
        else
            f = CO2_parameter * log((CO2ppm*1.01^t)/CO2_base)
        end
        f
    end

    function define_equation(dTdt, T, p, t) 
        """
        ODE/SDE to solve with/without seasonality and 
        with/without albedo feedback.
        """
        albedo_td, seasonality, historical_forcing = p 
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
        dTdt .= (-A*v - B*T + D * (laplacian * T) + Q*S.*a + 
                    forcing_CO2(t, historical_forcing)*v)/c

    end

    function g(dTdt, T,p,t)
        """
        Noise term in SDE. 
        """
        dTdt .= 0.05
    end

    function calculate_SIA(sol)
        """
        Calculates the sea ice area on the northern hemisphere. 
        """
        #ind = mapslices(t -> replace(u -> findfirst(x -> x < -10, t), Nothing => 0), X, dims = 1)
        ind = Array{Int16}(undef, size(sol)[2]) #zeros(size(X)[2])
        for i = 1:size(sol)[2]
            if isnothing(findfirst(x -> x < Tc, sol[:, i])) 
                ind[i] = n
            else 
                ind[i] = findfirst(x -> x < Tc, sol[:, i])
            end
        end
        area = 2 .* pi .* r_E^2 .* (1 .- theta[ind])
    end

    function calculate_global_average_T(sol) 
        """
        Calculating the average global temperature for every year, weighted by the areas of the latitude bands. 
        """
        weight_area = 2 .* pi .* r_E^2 .* (1 .- theta)
        global_temp = mean(sol, weights(weight_area), dims = 1)
    end

    function initialize_problem(p, noise, tspan)
        """
        Initialize the problem as ODE or SDE depending on input parameter noise.
        The equation gets solved for 20 years to get an initial guess for T,
        before running it for the whole tspan.  

        Parameters:
        p: tuple containing the booleans (albedo_td, seasonality, historical_forcing)
        noise: boolean if to solve ODE or SDE
        tspan: tuple of time span

        """

        if noise == true 
            W = WienerProcess(0.0,0.0,0.0)
            prob = SDEProblem(define_equation, g, T0, (0.0,20), p, noise = W)   
            sol =  solve(prob, ImplicitRKMil(), save_everystep=false, dt = 0.01, progress = true,
                progress_steps = 1)
            #prob = ODEProblem(define_equation, T0, (0.0,10), p)   
            #sol =  solve(prob, CVODE_BDF(), save_everystep=false, dt = 0.01, progress = true,
            #    progress_steps = 1)
            T0_new = sol[:, end]
            prob = SDEProblem(define_equation, g, T0_new, tspan, p, noise = W)
        else 
            prob = ODEProblem(define_equation, T0, (0.0,20), p)   
            sol =  solve(prob, CVODE_BDF(), save_everystep=false, dt = 0.01, progress = true,
                progress_steps = 1)
            T0_new = sol[:, end]
            prob = ODEProblem(define_equation, T0_new, tspan, p)   
        end
    end

    function solve_problem(p, noise, tspan, ensemble)
        """
        Solve the 1D (latitudinal) EBM on the northern hemisphere. 
        Returns the solution array.
        
        Parameters:
        p: tuple containing the booleans (albedo_td, seasonality, historical_forcing)
        noise: boolean if to solve ODE or SDE
        tspan: tuple of time span
        ensemble: boolean if ensemble run or single run
        
        """
        prob = initialize_problem(p, noise, tspan)
        if noise == true
            if ensemble == true 
                ensembleprob = EnsembleProblem(prob) 
                prob_no_noise = initialize_problem(p, false, tspan)
                # A run without noise is saved for analysis 
                sol_no_noise = solve(prob_no_noise, CVODE_BDF(), saveat = 1/12, progress = true,
                progress_steps = 1, dt=0.01)
                save_ensemble(sol_no_noise, "temp_no_noise") 

                sol = solve(ensembleprob,ImplicitRKMil(), saveat = 1/12, progress = true,
                progress_steps = 1, adaptive=false,dt=0.01, trajectories=runs)
            else
                @time begin
                sol = solve(prob, ImplicitRKMil(), saveat = 1/12, progress = true,
                            progress_steps = 1, adaptive=false,dt=0.01)
                end
            end
        else
            @time begin 
                # Rodas4(), or TRBDF2() seem to be the fastest solvers 
                # CVODE_BDF() from Sundials package is 10x faster: https://github.com/SciML/Sundials.jl
                # radau() from ODEInterfaceDiffEq package also very fast
                sol =  solve(prob, CVODE_BDF(), saveat = 1/12, dt = 0.01,progress = true,
                progress_steps = 1)
            end
        end
    end

    function load_historical_forcing() 
        """
        Function to load historical forcing data. 
        """
        forcing_data =  directory * "/input/Imbalance.txt"
        forcing = CSV.File(forcing_data,  skipto=2,  delim=" ", 
        ignorerepeated=true, drop = [1], limit = 131) |> Tables.matrix
        forcing_sum = sum(forcing, dims = 2)
    end

    function interpolate_historical_forcing(x) 
        """
        Function to interpolate historical forcing data. 
        """
        forcing =  load_historical_forcing()  
        itp = LinearInterpolation(collect(0:1:size(forcing, 1)-1), 
                                forcing[:, 1])
        itp(x)
    end

    function save_global_temp(sol) 
        """
        Function to save global temperature as CSV (single run).
        """
        CSV.write(directory * "/output/glob_temp_m.txt", convert(DataFrame, calculate_global_average_T(sol)'))
    end

    function save_ensemble(sol, filename) 
        """
        Function to save ensemble run as .jld (hdf5). 
        150 runs with 250 years (12 steps per year) ~ 1.8 GB.
        """
        @save directory * "/output/" * filename * "$(runs)r.jld" sol
    end


    function save_ensemble_csv(sol)
        """
        Function to save every ensemble run as .csv. 
        1 run with 250 years (12 steps per year) ~ 28 MB.
        """
        for i= 1:1:runs
            writedlm(directory * "/output/" * "run$(i).csv",  sol[:, :, i], ',')
        end
    end

    const laplacian = construct_laplacian()

end