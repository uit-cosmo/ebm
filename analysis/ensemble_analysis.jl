module EnsembleAnalysis 
    using PyCall
    using JLD2
    using DifferentialEquations
    using StatsBase, HypothesisTests
    include("../North.jl")
    include("../Plotting_functions.jl")
    using .North
    using .PlottingFunctions
    using FileIO

    directory = @__DIR__
    directory = directory[1:end-8]
    sol_no_noise = load(directory * "/output/" * "temp_no_noise10r.jld", "sol")
    JLD2.@load directory * "/output/" * "ensemble_sol_noalbedo10r.jld" sol

   
    function rolling_mean(array, window_size)
        """
        Function to calculate rolling mean using pandas from Python.
        """
        py"""
        import numpy as np 
        import pandas as pd

        def rolling_mean_py(array, window_size):
            mean = pd.Series(array).rolling(window = window_size).mean()
            return mean 
        """
        x = py"rolling_mean_py"(array, window_size)
    end

    function variance_autocorrelation(function_string, variable, window_size)
        """
        Function to calculate rolling variance/autorcorrelation using pandas from Python.
        """
        py"""
        import numpy as np 
        import pandas as pd
        from scipy import sparse
        from scipy.ndimage import gaussian_filter1d

        def rolling_statistics(function_string, variable_month, window_size):
            window_size = int(window_size/12)
            detrended_month = variable_month#variable_month - gaussian_filter1d(variable_month, sigma =  10, mode = "nearest")
            if function_string == "autocorrelation":
                #detrended_month = detrended_month[0:250]
                s = pd.Series(detrended_month)

                indicator_month = s.rolling(window_size).apply(lambda x: x.autocorr(), raw = False)


            elif function_string == "variance":
                indicator_month = pd.Series(detrended_month).rolling(window = window_size).var() 
            return np.array(indicator_month)
        
        """
        x = py"rolling_statistics"(function_string, variable, window_size)
    end


    function analysis()
        month = 3
        window_size = 12*50
        glob_temp = Array{Float32}(undef, size(sol, 2), size(sol, 3))
        sia = Array{Float32}(undef, size(sol, 2), size(sol, 3))
        sia_detrended = similar(sia)
        #sia_rolling = Array{Float32}(undef, size(sol, 2), size(sol, 3))
        glob_temp_year = Array{Float32}(undef, convert(Int16, floor(size(sol, 2)/12)), size(sol, 3))
        for i = 1:1:(size(sol, 3))
            #sol_detrended[:,:,i] = sol[:,:,i] .- sol_no_noise
            glob_temp[:, i] = North.calculate_global_average_T(sol[:,:,i])
            glob_temp_year[:, i] = mean(reshape(glob_temp[1:end-1, i], 12, :), dims = 1)
            sia[:, i] = North.calculate_SIA(sol[:,:,i])
            sia_detrended[:, i] = sia[:, i] - North.calculate_SIA(sol_no_noise[:,:])
            #sia_rolling[:, i] .= rolling_mean(sia[:, i], 12)
            #sia_var[:, i] = variance_autocorrelation("autocorrelation", sia[:, i], window_size, month)
        end
        sia_month = month_array(sia_detrended, month)

        quantiles_var, mean_sia_var = mean_quantiles_indicator(sia_month, "var", window_size, month)
        quantiles_ac, mean_sia_ac= mean_quantiles_indicator(sia_month, "ac", window_size, month)
        quantiles_sia, mean_sia= mean_quantiles_indicator(month_array(sia, month), "none", window_size, month)
        PlottingFunctions.confidence_plot(quantiles_var, convert(Int8, window_size/12), mean_sia_var, 25, "Variance [km²]", "variance_month$(month)")
        PlottingFunctions.confidence_plot(quantiles_ac, convert(Int8, window_size/12), mean_sia_ac, 25, "AC(1)", "ac1_month$(month)")
        PlottingFunctions.confidence_plot(quantiles_sia, 0, mean_sia, 25, "SIA [km²]", "sia_month$(month)")
        #PlottingFunctions.line_plot_ensemble(sia_var, "var SIA", 25)
        #PlottingFunctions.line_plot_ensemble(glob_temp_year, "T", 25)
        #PlottingFunctions.line_plot_ensemble(sia_rolling, "SIA", 25)
    end
    
    function calculate_quantile(array, p, window_size)
        array = array[window_size:end, :]
        array = replace(array, NaN=>missing)
        quantiles = Array{Float32}(undef, size(array, 1), 2)
        #quantiles[:, 2] = mapslices(x -> quantile(filter(!isnan, x), 1 - p), array, dims = 2)

        for i=1:1:size(array, 1)
            if all(ismissing.((array[i, :]))) 
                quantiles[i, 2] = NaN
                quantiles[i, 1] = NaN
            else 
                quantiles[i, 2] = quantile(skipmissing(array[i, :]), 1 - p)
                quantiles[i, 1] = quantile(skipmissing(array[i, :]), p)
            end
        end
        #quantiles[:, 2] = mapslices(x -> quantile(skipmissing(x), 1 - p), array, dims = 2)#quantile(array[1,:], 1 - p)
        #quantiles[:, 1] = mapslices(x -> quantile(skipmissing(x), p), array, dims = 2)
        
        return quantiles
    end

    function month_array(array, month)
        array[collect(month:12:end), :]
    end
    
    function mean_quantiles_indicator(variable, indicator, window_size, month)
        indicator_array = Array{Float32}(undef, size(variable, 1), size(variable, 2))
        
        if indicator == "var"
            mean_indicator = Array{Float32}(undef, size(indicator_array[50:end, :],1))
            for i=1:1:(size(variable, 2))
                indicator_array[:, i] = variance_autocorrelation("variance", variable[:, i], window_size)
            end
        elseif indicator == "ac"
            mean_indicator = Array{Float32}(undef, size(indicator_array[50:end, :],1))
            for i=1:1:(size(variable, 2))
                indicator_array[:, i] = variance_autocorrelation("autocorrelation", variable[:, i], window_size)
            end
        elseif indicator == "none"
            mean_indicator = Array{Float32}(undef, size(indicator_array, 1))
            quantiles_indicator = calculate_quantile(variable, 0.1, 1)
            mean!(mean_indicator, variable[:, :])

            return quantiles_indicator, mean_indicator
        end
        quantiles_indicator = calculate_quantile(indicator_array, 0.1, 50)
        mean!(mean_indicator, indicator_array[50:end, :])

        return quantiles_indicator, mean_indicator
    end

    analysis()

end