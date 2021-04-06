include("./North.jl")
using .North
using DifferentialEquations
using Plots
using StatsBase
using Logging: global_logger
using TerminalLoggers: TerminalLogger
global_logger(TerminalLogger())




tspan = (0.0,400.9)
albedo_td = true
seasonality = true
noise = true
p = (true, true)


prob = North.initialize_problem((albedo_td, seasonality), noise, tspan)
#sol = solve(prob, TRBDF2())
@time begin
    sol_noise = solve(prob, ImplicitRKMil(), saveat = 1/12, progress = true,
    progress_steps = 1, adaptive=false,dt=0.01)
end

plt2 = plot(sol_noise[:, end]);
display(plt2)
@time begin
    sia = North.calculate_SIA(sol_noise)
end
@time begin
    glob_temp = North.calculate_global_average_T(sol_noise)
end
@time begin
    glob_temp_year = mean(reshape(glob_temp, :, 12), dims = 2)
end
println(size(glob_temp_year))
display(plot(glob_temp_year[:, 1]))
display(plot(sia))
#println(sia)
@gif for i ∈ 1:12:size(sol_noise)[2]
    plot(sol_noise[:, i])
end