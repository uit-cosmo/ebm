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

tspan = (0.0,130)
albedo_td = true
seasonality = true
historical_forcing =  true
noise = true

println(@__DIR__)

sol_noise =  North.solve_problem((albedo_td, seasonality, historical_forcing), noise, tspan)
North.save_global_temp(sol_noise) 
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
    PlottingFunctions.plot_historical_simulation(glob_temp_year[1, :] .- mean(glob_temp_year))
    PlottingFunctions.plot_sia_historical_simulation(sia) 
end

@gif for i ∈ 1:12:size(sol_noise)[2]
    plot(sol_noise[:, i], ylims = (-20,30))
end
"""
PlottingFunctions.simple_line_plot(glob_temp_year[1, :], 
                L"Global Temperature [$^\circ C$]", 10)
PlottingFunctions.simple_line_plot(glob_temp_year[1, :] .- mean(glob_temp_year), 
                L"Global Temperature Anomaly [$^\circ C$]", 10)
PlottingFunctions.simple_line_plot(sia, "SIA [km²]", 50)


@gif for i ∈ 1:12:size(sol_noise)[2]
    plot(sol_noise[:, i], ylims = (-20,30))
end
"""