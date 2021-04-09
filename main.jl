##########################################################
#This file calls all the functions to solve the system.  #
#You can specify the parameters:                         #
#albedo/seasonality/historical_forcing/noise/ensemble    #
##########################################################


include("./North.jl")
include("./Plotting_functions.jl")
using .North
using .PlottingFunctions
using DifferentialEquations
using Plots
using StatsBase
using Logging: global_logger
using TerminalLoggers: TerminalLogger
using LaTeXStrings
global_logger(TerminalLogger())

albedo_td = true
seasonality = true
historical_forcing =  false
noise = true
ensemble = true

if historical_forcing == true
    tspan = (0.0,130)
else
    tspan = (0.0,250)
end

println(@__DIR__)

sol_noise =  North.solve_problem((albedo_td, seasonality, historical_forcing), noise, tspan, ensemble)
if ensemble == true 
    North.save_ensemble(sol_noise)
end

#North.save_global_temp(sol_noise) 

@time begin
    sia = North.calculate_SIA(sol_noise)
end
@time begin
    glob_temp = North.calculate_global_average_T(sol_noise)
end
@time begin
    glob_temp_year = mean(reshape(glob_temp[:, 1:end-1], 12, :), dims = 1)
end

if historical_forcing == true 
    PlottingFunctions.plot_historical_simulation(glob_temp_year[1, 1:end] .- mean(glob_temp_year[1, 1:end]))
    PlottingFunctions.plot_sia_historical_simulation(sia[1*12:end]) 
else
    PlottingFunctions.simple_line_plot(glob_temp_year[1, :], 
        L"Global Temperature [$^\circ C$]", 10)
    PlottingFunctions.simple_line_plot(glob_temp_year[1, :] .- mean(glob_temp_year), 
        L"Global Temperature Anomaly [$^\circ C$]", 10)
    PlottingFunctions.simple_line_plot(sia, "SIA [km²]", 50)
end

@gif for i ∈ 1:12:size(sol_noise)[2]
    plot(sol_noise[:, i], ylims = (-20,30))
end