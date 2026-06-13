using DrWatson
using GLMakie
using JSON
using Statistics

@quickactivate "neural_field"
include(srcdir("parallel_delayed.jl"))
include(srcdir("data_processing.jl"))

############################
# Strong Scaling plot 
############################

function plot_strong_scaling_nodes()
    proc_list = [1,2,4,8,16,32,64]
    nxny_arr = [30,60]
    num_reps = 10
    num_layers = 2
    best_time_arr_combined = []

    for nx in nxny_arr
        best_time_arr = get_strong_scaling_data(proc_list, num_reps, nx, nx, num_layers)
        best_time_arr_round = round.(best_time_arr, digits = 2)
        push!(best_time_arr_combined, best_time_arr_round)
    end 

    fig = Figure()
    ax = Axis(
        fig[1,1]; 
        title = "Strong Scaling - Nodes", 
        xlabel = "Number of processors", 
        ylabel = "Wall clock time (s)", 
        xticks = proc_list, 
        xscale = log2, 
        yscale = log2
    )

    for entry in eachindex(nxny_arr)
        nx = nxny_arr[entry]
        time_arr = best_time_arr_combined[entry]

        perfect_arr = time_arr[1] ./ proc_list
        num_el = nx*nx
        lines!(ax, proc_list, time_arr, label = "Actual scaling, num_elements =$num_el")
        scatter!(ax, proc_list, time_arr)

        perfect_label = entry == firstindex(nxny_arr) ? "Perfect scaling" : nothing
        lines!(ax, proc_list, perfect_arr, linestyle = :dash, color = :gray,label = perfect_label)
    end

    axislegend(ax)
    save("plots/strong_scaling_nodes.png",fig)
end

function plot_strong_scaling_layers()
    proc_list = [1,2,4,8,16,32,64]
    nx = ny = 50
    num_reps = 10
    num_layers_arr = [2,4,8]
    best_time_arr_combined = []

    for layer in num_layers_arr
        best_time_arr = get_strong_scaling_data(proc_list, num_reps, nx, ny, layer)
        best_time_arr_round = round.(best_time_arr, digits = 2)
        push!(best_time_arr_combined, best_time_arr_round)
    end 

    fig = Figure()
    ax = Axis(
        fig[1,1]; 
        title = "Strong Scaling - Layers", 
        xlabel = "Number of processors", 
        ylabel = "Wall clock time (s)", 
        xticks = proc_list, 
        xscale = log2, 
        yscale = log2
    )

    for entry in eachindex(num_layers_arr)
        layer = num_layers_arr[entry]
        time_arr = best_time_arr_combined[entry]

        perfect_arr = time_arr[1] ./ proc_list

        lines!(ax, proc_list, time_arr, label = "Actual scaling, $(layer) layers")
        scatter!(ax, proc_list, time_arr)

        perfect_label = entry == firstindex(nxny_arr) ? "Perfect scaling" : nothing
        lines!(ax, proc_list, perfect_arr, linestyle = :dash, color = :gray,label = perfect_label)
    end

    axislegend(ax)
    save("plots/strong_scaling_layer.png",fig)
end

############################
# Efficiency plot
############################

function plot_efficiency_nodes()
    nxny_arr = [30,60]
    num_layers = 2
    num_reps = 10
    proc_list = [1,2,4,8,16,32,64]
    efficiency_combined = []
    for nx in nxny_arr
        t1 = get_t1_time(nx, nx, num_layers)
        efficiency = get_efficiency_data(proc_list, t1, num_reps, nx, nx, num_layers)
        push!(efficiency_combined,efficiency)
    end 

    fig = Figure()
    ax = Axis(
        fig[1,1]; 
        title = "Efficiency - nodes", 
        xlabel = "Number of processors", 
        ylabel = "Efficiency", 
        xticks = proc_list,
        yticks = 0:0.2:1.2,
        xscale = log2,
        limits = (nothing, (0, 1.2))
    )

    perfect_eff = [1 for _ in proc_list]
    lines!(ax,proc_list, perfect_eff, color = :gray, label = "Perfect Efficiency")

    for entry in eachindex(nxny_arr)
        nx = nxny_arr[entry]
        num_el = nx*nx
        lines!(ax, proc_list, efficiency_combined[entry], label = "Efficiency, num_elements =$num_el")
        scatter!(ax, proc_list, efficiency_combined[entry])
    end

    axislegend(ax; position = :lb )
    save("plots/efficiency_nodes.png",fig)
end 

function plot_efficiency_layers()
    nx = ny = 50
    num_layers = [2,4,8]
    num_reps = 10
    proc_list = [1,2,4,8,16,32,64]
    efficiency_combined = []
    for layer in num_layers
        t1 = get_t1_time(nx, nx, layer)
        efficiency = get_efficiency_data(proc_list, t1, num_reps, nx, ny, layer)
        push!(efficiency_combined,efficiency)
    end 

    fig = Figure()
    ax = Axis(
        fig[1,1]; 
        title = "Efficiency - nodes", 
        xlabel = "Number of processors", 
        ylabel = "Efficiency", 
        xticks = proc_list,
        yticks = 0:0.2:1.4,
        xscale = log2,
        limits = (nothing, (0, 1.4))
    )

    perfect_eff = [1 for _ in proc_list]
    lines!(ax,proc_list, perfect_eff, color = :gray, label = "Perfect Efficiency")

    for (entry, layer) in enumerate(num_layers)
        lines!(ax, proc_list, efficiency_combined[entry], label = "Efficiency, $layer layers")
        scatter!(ax, proc_list, efficiency_combined[entry])
    end

    axislegend(ax; position = :lb )
    save("plots/efficiency_layers.png",fig)
end

############################
# Computation/Communication plot 
############################

procs_arr = [1,2,4,8,16,32,64] 
comm_arr_full = []
comp_arr_full = []
num_reps = 10

for proc in procs_arr
    data_eff = JSON.parsefile("data/exp_raw/mpi_exp/results_nx30ny30np$(proc)num_layers2.json")
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

save("plots/communication_vs_computation_30.png",fig)

# plot all 

plot_strong_scaling_nodes()
plot_strong_scaling_layers()
plot_efficiency_nodes()
plot_efficiency_layers()