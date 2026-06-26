using DrWatson
using GLMakie

@quickactivate "neural_field"

include(srcdir("data_processing.jl"))

function comm_box_plot(num_proc)
    num_reps = 50
    nx=ny=19
    layers_list = [16,32,64]
    comm_arr = []
    comp_arr = []

    for layer in layers_list
        comm, comp, barr = get_max_comm_comp_time_over_repetitions(num_proc,num_reps,nx,ny,layer)
        push!(comm_arr, comm.+barr)
        push!(comp_arr, comp)

    end 

    fig = Figure()

    ax = Axis(
        fig[1,1], 
        xticks = (1:3, ["16 layers", "32 layers", "64 layers"]), 
        limits = (nothing, (0,0.005)),
    )

    for i in eachindex(comm_arr)
        boxplot!(ax, fill(i, length(comm_arr[i])), comm_arr[i])
    end

    save("plots/boxplot_num_proc$(num_proc)barr_and_comm.png",fig)
end 

function plot_layer_comm_comp_barr_time(num_proc)
    comm_arr = []
    comp_arr = []
    barrier_arr = []
    layers_list = [2,4,8,16,32,64,128,256,512]
    nx=ny=19
    num_reps = 50

    for layer in layers_list
        comm, comp, barr = get_max_comm_comp_time_over_repetitions(num_proc,num_reps,nx,ny,layer)
        push!(comm_arr, comm)
        push!(comp_arr, comp)
        push!(barrier_arr, barr)
    end 

    fig = Figure()

    ax = Axis(
        fig[1,1], 
        xlabel = "Number of layers", 
        ylabel = "Wall Clock Time (s)",
        xticks= layers_list,
        xscale = log2,
        yscale = log2)

    median_comm = median.(comm_arr)
    max_comm = maximum.(comm_arr)
    min_comm = minimum.(comm_arr)

    median_barr = median.(barrier_arr)
    max_barr = maximum.(barrier_arr)
    min_barr = minimum.(barrier_arr)

    median_comp = median.(comp_arr)

    band!(ax, layers_list, min_comm, max_comm, color=(:blue,0.18))
    band!(ax, layers_list, min_barr, max_barr, color=(:yellow,0.18))

    lines!(ax, layers_list, median_comm, label = "Median All reduce time")
    lines!(ax, layers_list, median_barr, label = "Median Barrier time")
    lines!(ax, layers_list, median_comp, label = "Median Computation time")
    axislegend(ax, position =:lt)

    save("plots/comm_comp_barr_time_nx$(nx)_proc$(num_proc).png",fig)
end     

function plot_nodes_comp_comm_barr_time(num_proc)
    comm_arr = []
    comp_arr = []
    barrier_arr = []
    num_layer = 2
    nx_arr = [10,20,40,80]
    num_reps = 50

    for nx in nx_arr
        comm, comp, barr = get_max_comm_comp_time_over_repetitions(num_proc,num_reps,nx,nx,num_layer)
        push!(comm_arr, comm)
        push!(comp_arr, comp)
        push!(barrier_arr, barr)
    end 

    fig = Figure()

    ax = Axis(
        fig[1,1], 
        xlabel = "nx=ny", 
        ylabel = "Wall Clock Time (s)",
        xticks= nx_arr,
        xscale = log2,
        yscale = log2)

    median_comm = median.(comm_arr)
    max_comm = maximum.(comm_arr)
    min_comm = minimum.(comm_arr)

    median_barr = median.(barrier_arr)
    max_barr = maximum.(barrier_arr)
    min_barr = minimum.(barrier_arr)

    median_comp = median.(comp_arr)

    band!(ax, nx_arr, min_comm, max_comm, color=(:blue,0.18))
    band!(ax, nx_arr, min_barr, max_barr, color=(:yellow,0.18))

    lines!(ax, nx_arr, median_comm, label = "Median All reduce time")
    lines!(ax, nx_arr, median_barr, label = "Median Barrier time")
    lines!(ax, nx_arr, median_comp, label = "Median Computation time")
    axislegend(ax, position =:lt)

    save("plots/node_comm_comp_barr_time_layer$(num_layer)_proc$(num_proc).png",fig)
end

comm_box_plot(32)

plot_layer_comm_comp_barr_time(16)
plot_layer_comm_comp_barr_time(32)
plot_layer_comm_comp_barr_time(64)

plot_nodes_comp_comm_barr_time(16)
plot_nodes_comp_comm_barr_time(32)
plot_nodes_comp_comm_barr_time(64)