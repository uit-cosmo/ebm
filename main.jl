include("./North.jl")
using .North
using DifferentialEquations
using Plots
using StatsBase
using Logging: global_logger
using TerminalLoggers: TerminalLogger
global_logger(TerminalLogger())

tspan = (0.0,130)
albedo_td = false
seasonality = true
historical_forcing =  true
noise = false

println(@__DIR__)

sol_noise =  North.solve_problem((albedo_td, seasonality, historical_forcing), noise, tspan)


@time begin
    sia = North.calculate_SIA(sol_noise)
end
@time begin
    glob_temp = North.calculate_global_average_T(sol_noise)
end
@time begin
    glob_temp_year = mean(reshape(glob_temp[:, 1:end-1], 12, :), dims = 1)
end

println(size(glob_temp_year))
display(plot(glob_temp_year[1, :]))
display(plot(glob_temp_year[1, :] .- mean(glob_temp_year)))
display(plot(sia))
#println(sia)
@gif for i ∈ 1:12:size(sol_noise)[2]
    plot(sol_noise[:, i], ylims = (-20,30))
end