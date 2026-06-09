using DrWatson
using GLMakie
using JSON
using Statistics

@quickactivate "neural_field"
include(srcdir("parallel_delayed.jl"))

nodes_array = [2,4,8,16,32,64]
num_reps = 10
time_arr = []
mean_arr = []
max_arr_all = Float64[]
min_arr_all = Float64[]

for proc in nodes_array
    data = JSON.parsefile("data/exp_raw/mpi_exp/results_nx30ny30np$proc.json")
    rhs_times = [entry["rhs"] for entry in data]
    # get summed times (all rhs parts)
    summed_arr = []
    for num in 1:proc
        summed = sum.(rhs_times[num])
        push!(summed_arr, summed)
    end 
    # get max per run
    max_arr = []
    for num in 1:num_reps
        max_time = maximum(summed_arr[rank][num] for rank in 1:proc)
        push!(max_arr, max_time)
    end 
    mean_time = mean(max_arr)
    # statistics 
    push!(min_arr_all, minimum(max_arr))
    push!(max_arr_all, maximum(max_arr))
end 

mean_arr_round = round.(mean_arr, digits = 2)

fig = Figure()
ax = Axis(
    fig[1,1]; 
    title = "Strong Scaling", 
    xlabel = "Number of processors", 
    ylabel = "Wall clock time", 
    xticks = nodes_array, 
    #yticks = mean_arr_round,
    xscale = log2, 
    yscale = log2
)

band!(
    ax,
    nodes_array,
    min_arr_all,
    max_arr_all,
    color = (:purple, 0.20)
)

lines!(ax, nodes_array, mean_arr_round, color = :purple)
scatter!(ax, nodes_array, mean_arr_round, color = :purple)
save("plots/strong_scaling.png",fig)