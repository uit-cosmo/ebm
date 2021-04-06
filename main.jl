include("./North.jl")
using .North
using DifferentialEquations
using Plots
using StatsBase
using Logging: global_logger
using TerminalLoggers: TerminalLogger
global_logger(TerminalLogger())

tspan = (0.0,41)
albedo_td = true
seasonality = true
noise = false

println(@__DIR__)

sol_noise =  North.solve_problem((albedo_td, seasonality), noise, tspan)

plt2 = plot(sol_noise[:, end]);
display(plt2)
@time begin
    sia = North.calculate_SIA(sol_noise)
end
@time begin
    glob_temp = North.calculate_global_average_T(sol_noise)
end
@time begin
    glob_temp_year = mean(reshape(glob_temp[:, 1:end-1], :, 12), dims = 2)
end

println(size(glob_temp_year))
display(plot(glob_temp_year[:, 1]))
display(plot(sia))
#println(sia)
@gif for i ∈ 1:12:size(sol_noise)[2]
    plot(sol_noise[:, i], ylims = (-20,30))
end