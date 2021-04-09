module EnsembleAnalysis 
    using PyCall
    using JLD2
    using DifferentialEquations
    using StatsBase, HypothesisTests
    include("../North.jl")
    include("../Plotting_functions.jl")
    using .North
    using .PlottingFunctions
    
    directory = @__DIR__
    directory = directory[1:end-8]
    JLD2.@load directory * "/output/" * "ensemble_sol_100r.jld" sol


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

    function variance_autocorrelation(function_string, variable, window_size, month)
        """
        Function to calculate rolling variance/autorcorrelation using pandas from Python.
        """
        py"""
        import numpy as np 
        import pandas as pd
        from scipy import sparse
        from scipy.ndimage import gaussian_filter1d

        def month_array(month, array): 
            return array[np.arange(month, array.size, 12)]

        def rolling_statistics(function_string, variable, window_size, month):
            window_size = int(window_size/12)
            variable_month = month_array(month, variable)
            detrended_month = variable_month - gaussian_filter1d(variable_month, sigma =  10, mode = "nearest")
            if function_string == "autocorrelation":
                indicator_month = pd.Series(detrended_month).rolling(window = window_size).apply(lambda x: x.autocorr(), raw = False)
            elif function_string == "variance":
                indicator_month = pd.Series(detrended_month).rolling(window = window_size).var() 
            return np.array(indicator_month)
        
        """
        x = py"rolling_statistics"(function_string, variable, window_size, month)
    end


    function analysis()
        window_size = 12*50
        glob_temp = Array{Float32}(undef, size(sol, 2), size(sol, 3))
        sia = Array{Float32}(undef, size(sol, 2), size(sol, 3))
        sia_rolling = Array{Float32}(undef, size(sol, 2), size(sol, 3))
        sia_var = Array{Float32}(undef, convert(Int16, floor(size(sol, 2)/12)), size(sol, 3))
        glob_temp_year = Array{Float32}(undef, convert(Int16, floor(size(sol, 2)/12)), size(sol, 3))
        for i = 1:1:(size(sol, 3))

            glob_temp[:, i] = North.calculate_global_average_T(sol[:,:,i])
            glob_temp_year[:, i] = mean(reshape(glob_temp[1:end-1, i], 12, :), dims = 1)
            sia[:, i] = North.calculate_SIA(sol[:,:,i])
            sia_rolling[:, i] .= rolling_mean(sia[:, i], 12)
            sia_var[:, i] = variance_autocorrelation("variance", sia[:, i], window_size, 3)
        end
        quantiles_sia = calculate_quantile(sia, 0.1, 12)
        quantiles_var = calculate_quantile(sia_var, 0.1, 50)
        mean_sia_var = mean(sia_var[50:end, :], dims = 2)
        mean_sia = mean(sia, dims = 2)
        #PlottingFunctions.confidence_plot(quantiles_var, mean_sia_var, 25, "variance")
        PlottingFunctions.confidence_plot(quantiles_sia, mean_sia, 250, "sia")
        #PlottingFunctions.line_plot_ensemble(sia_var, "var SIA", 25)
        #PlottingFunctions.line_plot_ensemble(glob_temp_year, "T", 25)
        #PlottingFunctions.line_plot_ensemble(sia_rolling, "SIA", 25)
    end
    
    function calculate_quantile(array, p, window_size)
        array = array[window_size:end, :]
        quantiles = Array{Float32}(undef, size(array, 1), 2)
        quantiles[:, 2] = mapslices(x -> quantile(x, 1 - p), array, dims = 2)#quantile(array[1,:], 1 - p)
        quantiles[:, 1] = mapslices(x -> quantile(x, p), array, dims = 2)
        
        return quantiles
    end
    
    analysis()

end