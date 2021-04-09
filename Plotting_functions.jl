module PlottingFunctions 
    ###############################################################
    #This module contains the functions to plot the solution etc. #
    #PyCall/PyPlot/matplotlib should be installed                 #
    ###############################################################

    using PyCall
    import PyPlot; const plt = PyPlot
    pygui(true)
    using CSV
    using Tables
    rcParams = PyCall.PyDict(plt.matplotlib["rcParams"])
    rcParams["font.size"] = 10
    rcParams["axes.linewidth"] = 2
    rcParams["axes.labelsize"] = 10
    rcParams["axes.labelweight"] = "bold"
    rcParams["font.weight"] = "bold"

    function simple_line_plot(array, y_label, steps_years)
        fig, ax = plt.subplots()
        ax.set_xlabel("Modeltime [year]")
        xticklabels = [string(x) for x = 0:steps_years:size(array, 1)]
        ax.set_xticks(collect(0:steps_years:size(array, 1)))
        ax.set_xticklabels(xticklabels)
        ax.grid(true, ls = "--")
        ax.set_ylabel(y_label)
        ax.plot(array, lw = 2)
        plt.show()
    end

    function line_plot_ensemble(array, y_label, steps_years) 
        runs = size(array, 2)
        fig, ax = plt.subplots()
        ax.set_xlabel("Modeltime [year]")
        xticklabels = [string(x) for x = 0:steps_years:size(array, 1)]
        ax.set_xticks(collect(0:steps_years:size(array, 1)))
        ax.set_xticklabels(xticklabels)
        ax.grid(true, ls = "--")
        ax.set_ylabel(y_label)
        for i = 1:1:runs
            ax.plot(array[:, i], lw = 2)
        end
        plt.show()
    end

    function confidence_plot(quantiles, mean, steps_years, y_label)
        fig, ax = plt.subplots()
        ax.set_xlabel("Modeltime [year]")
        xticklabels = [string(x) for x = 0:steps_years:size(quantiles, 1)]
        ax.set_xticks(collect(0:steps_years:size(quantiles, 1)))
        ax.set_xticklabels(xticklabels)
        ax.grid(true, ls = "--")
        ax.set_ylabel(y_label)
        ax.fill_between(collect(1:1:size(quantiles, 1)), quantiles[:, 1], quantiles[:, 2], alpha = 0.5)
        ax.plot(collect(1:1:size(mean, 1)), mean, lw = 2)

        plt.show()
    end

    function plot_historical_simulation(array_sim)
        forcing = load_historical_forcing() 
        temp = load_historical_temperature()

        fig, ax = plt.subplots()
        ax2 = ax.twinx()
        ax.set_xlabel("Modeltime [year]")
        xticklabels = [string(x) for x = 0:10:size(array_sim, 1)]
        ax.set_xticks(collect(0:10:size(array_sim, 1)))
        ax.set_xticklabels(xticklabels)
        ax.grid(true, ls = "--")
        ax.set_ylabel("Temperature [K]")
        ax2.set_ylabel("Forcing [W/m²]")
        ax.plot(collect(size(temp, 1)-size(array_sim, 1):1:size(temp, 1)-1),array_sim, lw = 2, label = "sim. T")
        ax.plot(temp, lw = 2, ls = "--", label = "hist. T")
        ax2.plot(forcing, lw = 2, color = "C2", label = "forcing")
        fig.legend()

        plt.show()
    end


    function plot_sia_historical_simulation(array_sim) 
        sia = load_ice_area_cmip6() * 1e6
        sia = sia[30*12:1:160*12]

        fig, ax = plt.subplots()
        ax.set_xlabel("Modeltime [year]")
        xticklabels = [string(x) for x = 0:25:convert(Int16, trunc(size(array_sim, 1)/12))]
        ax.set_xticks(collect(0:25*12:size(array_sim, 1)))
        ax.set_xticklabels(xticklabels)
        ax.grid(true, ls = "--")
        ax.set_ylabel("SIA [km²]")
        ax.plot(collect(size(sia, 1)-size(array_sim, 1):1:size(sia, 1)-1),array_sim, lw = 2, label = "sim. SIA")
        ax.plot(sia, lw = 2, ls = "--", label = "hist. SIA")
        fig.legend()

        plt.show()
    end

    function load_ice_area_cmip6() 
        directory = @__DIR__
        sia_data =  directory * "/input/" * "ACCESS-CM2.r1i1p1f1" * ".txt"
        sia = CSV.File(sia_data, delim=" ") |> Tables.matrix
    end
    
    function load_historical_forcing() 
        directory = @__DIR__
        forcing_data =  directory * "/input/" * "Imbalance" * ".txt"
        forcing = CSV.File(forcing_data,  skipto=2,  delim=" ", 
        ignorerepeated=true, drop = [1], limit = 131) |> Tables.matrix
        forcing = sum(forcing, dims = 2)
    end

    function load_historical_temperature()
        directory = @__DIR__
        temp_data =  directory * "/input/" * "HadCrut3gl" * ".txt"
        temp = CSV.File(temp_data,  skipto=2) |> Tables.matrix
    end
end