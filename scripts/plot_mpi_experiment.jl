using DrWatson
using GLMakie
using JSON
using Statistics

@quickactivate "neural_field"
include(srcdir("parallel_delayed.jl"))

### plot Strong Scaling 

nodes_array = [2,4,8,16,32,64]
num_nodes = length(nodes_array)
num_reps = 10
max_min_arr = Float64[]
max_arr_all = Float64[]
min_arr = Float64[]

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
    push!(min_arr, minimum(max_arr))
    push!(max_arr_all, maximum(max_arr))
    push!(max_min_arr, minimum(max_arr))    
end 

perfect_arr = [max_min_arr[1]/(2^(n-1)) for n in 1:num_nodes]
max_min_arr_round = round.(max_min_arr, digits = 2)

fig = Figure()
ax = Axis(
    fig[1,1]; 
    title = "Strong Scaling", 
    xlabel = "Number of processors", 
    ylabel = "Wall clock time (s)", 
    xticks = nodes_array, 
    xscale = log2, 
    yscale = log2
)


lines!(ax, nodes_array, perfect_arr, color = :gray, label = "Perfect scaling")
scatter!(ax, nodes_array, perfect_arr, color = :gray)
lines!(ax, nodes_array, max_min_arr_round, color = :purple, label = "Actual scaling")
scatter!(ax, nodes_array, max_min_arr_round, color = :purple)
axislegend(ax)
save("plots/strong_scaling.png",fig)


### Plot efficiency 
procs_arr = [2,4,8,16,32,64] 
efficiency = []
# get reference value for T1
data_t1 = JSON.parsefile("data/exp_raw/mpi_exp/results_nx30ny30np1.json")
t1 = maximum(sum.(data[1]["rhs"]))
num_reps = 10

for proc in procs_arr
    data_eff = JSON.parsefile("data/exp_raw/mpi_exp/results_nx30ny30np$proc.json")
    rhs_times = [entry["rhs"] for entry in data_eff]
    summed_arr = []
    for num in 1:proc
        summed = sum.(rhs_times[num])
        push!(summed_arr, summed)
    end
    slowest_worker = []
    for num in 1:num_reps
        slowest = maximum(summed_arr[rank][num] for rank in 1:proc)
        push!(slowest_worker, slowest)
    end 
    tp = minimum(slowest_worker)
    eff = t1/(proc*tp)
    push!(efficiency,eff)
end     

fig = Figure()
ax = Axis(
    fig[1,1]; 
    title = "Efficiency", 
    xlabel = "Number of processors", 
    ylabel = "Efficiency", 
    xticks = procs_arr, 
    xscale = log2
)

perfect_eff = [1 for _ in procs_arr]

lines!(ax,procs_arr, efficiency, color = :green)
scatter!(ax,procs_arr, efficiency, color = :green)
lines!(ax,procs_arr, perfect_eff, color = :gray)
save("plots/efficiency.png",fig)

#### Plot Computation/Comunication plot 
procs_arr = [2,4,8,16,32,64] 
comm_arr_full = []
comp_arr_full = []
num_reps = 10

for proc in procs_arr
    data_eff = JSON.parsefile("data/exp_raw/mpi_exp/results_nx30ny30np$proc.json")
    rhs_times = [entry["rhs"] for entry in data_eff]
    # extract communication
    communication_arr = [Float64[] for _ in 1:proc] 
    for i in 1:proc
        for rep in 1:num_reps
        comm = rhs_times[i][rep][7]
        push!(communication_arr[i], comm)
        end
    end 
    comp_arr = []
    for num in 1:proc
        summed = sum.(rhs_times[num]).-communication_arr[num]
        push!(comp_arr, summed)
    end
    slowest_worker_comp = []
    slowest_worker_comm = []
    for num in 1:num_reps
        slowest_comp = maximum(comp_arr[rank][num] for rank in 1:proc)
        slowest_comm = maximum(communication_arr[rank][num] for rank in 1:proc)
        push!(slowest_worker_comp, slowest_comp)
        push!(slowest_worker_comm, slowest_comm)
    end
    mean_comp = mean(slowest_worker_comp)
    mean_comm = mean(slowest_worker_comm)
    push!(comm_arr_full, mean_comm)
    push!(comp_arr_full, mean_comp)
end

total = comm_arr_full .+ comp_arr_full
ratio_comm = comm_arr_full ./total 
ratio_comp = comp_arr_full./total
total_1 = [1 for _ in total]

fig = Figure(size = (1000, 450))

x = 1:length(procs_arr)

ax1 = Axis(
    fig[1, 1],
    xlabel = "Number of processors",
    ylabel = "Ratio of Compunication and Computation",
    title = "Computation vs communication time",
    xticks = (x, string.(procs_arr))
)


barplot!(
    ax1,
    x,
    ratio_comp,
    color = :steelblue,
    label = "Computation"
)

barplot!(
    ax1,
    x,
    total_1,
    fillto = ratio_comp,
    color = :orange,
    label = "Communication"
)

axislegend(ax1, position = :rt)

ax2 = Axis(
    fig[1, 2], 
    xlabel = "Number of processors",
    ylabel = "Time [s]",
    title = "Computation vs communication time",
    xticks = procs_arr, 
    xscale = log2
)

lines!(ax2, procs_arr, comm_arr_full, label = "Communication")
scatter!(ax2, procs_arr, comm_arr_full)
lines!(ax2, procs_arr, comp_arr_full, label = "Computation")
scatter!(ax2, procs_arr, comp_arr_full)

axislegend(ax2, position = :rt)

save("plots/communication_vs_computation.png",fig)