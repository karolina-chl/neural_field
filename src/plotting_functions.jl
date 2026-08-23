using DrWatson
using CairoMakie
using JSON
using Statistics

@quickactivate "neural_field"

include("parallel_delayed.jl")
include("data_processing.jl")

function result_file(nx, ny, np, num_layers)
    return datadir("exp_raw","mpi_exp","results_nx$(nx)ny$(ny)np$(np)num_layers$(num_layers).json",
    )
end

function get_available_strong_scaling_data(proc_list, num_reps, nx, ny, num_layers)
    available_procs = Int[]
    time_arr = Float64[]

    for proc in proc_list
        file = result_file(nx, ny, proc, num_layers)

        if isfile(file)
            t = get_strong_scaling_data([proc], num_reps, nx, ny, num_layers)[1]

            push!(available_procs, proc)
            push!(time_arr, t)
        else
            @warn "Missing file, skipping" file
        end
    end

    return available_procs, time_arr
end

function get_available_efficiency_data(proc_list, num_reps, nx, ny, num_layers)
    available_procs = Int[]

    for proc in proc_list
        file = result_file(nx, ny, proc, num_layers)

        if isfile(file)
            push!(available_procs, proc)
        else
            @warn "Missing file, skipping" file
        end
    end

    if isempty(available_procs)
        return Int[], Float64[], nothing
    end

    if isfile(result_file(nx, ny, 1, num_layers))
        p_ref = 1
        t_ref = get_t1_time(nx, ny, num_layers)
    else
        p_ref = available_procs[1]
        t_ref = get_strong_scaling_data([p_ref], num_reps, nx, ny, num_layers)[1]
    end

    efficiency_arr = Float64.(
        get_efficiency_data(
            available_procs,
            p_ref * t_ref,
            num_reps,
            nx,
            ny,
            num_layers,
        )
    )

    return available_procs, efficiency_arr, p_ref
end

function plot_strong_scaling_nodes(
    nxny_arr,
    proc_list,
    num_layers,
    num_reps;
    palette=[colorant"#C77A2F",colorant"#2F836C",colorant"#63C7C9", colorant"#A63E33"],
    )
    fig = Figure(fontsize = 18)

    ax = Axis(
        fig[1, 1];
        xlabel = "Number of processes",
        ylabel = "Wall clock time (s)",
        xticks = proc_list,
        xscale = log2,
        yscale = log2,
    )

    perfect_label_added = false

    for (i,(nx, ny)) in enumerate(nxny_arr)
        available_procs, time_arr = get_available_strong_scaling_data(
            proc_list,
            num_reps,
            nx,
            ny,
            num_layers,
        )

        if isempty(available_procs)
            @warn "No data available" nx ny
            continue
        end

        p0 = available_procs[1]
        t0 = time_arr[1]

        perfect_arr = t0 .* p0 ./ available_procs

        perfect_label = perfect_label_added ? nothing : "Perfect scaling"

        N = (nx+1)*(ny+1)
        M = N+num_layers*N^2
        M_bln = round(M/1000000000, digits=2)

        lines!(
            ax,
            available_procs,
            perfect_arr;
            linestyle = :dash,
            color = :gray,
            label = perfect_label,
        )

        perfect_label_added = true

        lines!(
            ax,
            available_procs,
            time_arr;
            label = "N= $N, $(M_bln) bln unknowns",
            color = palette[i],
        )

        scatter!(
            ax,
            available_procs,
            time_arr,
            markersize=13,
            color = :transparent,
            strokecolor = palette[i],
            strokewidth = 2,
        )
    end

    axislegend(ax, position = :lb)

    mkpath("plots")
    save("plots/large_strong_scaling_nodes_num_layers$(num_layers)_reps$(num_reps).pdf", fig)
end

function plot_strong_scaling_layers(
    nx,
    ny,
    proc_list,
    num_layers_arr,
    num_reps;
    palette=[colorant"#C77A2F",colorant"#2F836C",colorant"#63C7C9", colorant"#A63E33"],
    )

    fig = Figure(fontsize=18)

    ax = Axis(
        fig[1, 1];
        xlabel = "Number of MPI processes",
        ylabel = "Wall clock time (s)",
        xticks = proc_list,
        xscale = log2,
        yscale = log2,
    )


    
    perfect_label_added = false

    for (i,num_layers) in enumerate(num_layers_arr)
        available_procs, time_arr = get_available_strong_scaling_data(
            proc_list,
            num_reps,
            nx,
            ny,
            num_layers,
        )

        if isempty(available_procs)
            @warn "No data available" nx ny num_layers
            continue
        end

        p0 = available_procs[1]
        t0 = time_arr[1]

        perfect_arr = t0 .* p0 ./ available_procs

        perfect_label = perfect_label_added ? nothing : "Perfect scaling"

        N = (nx + 1) * (ny + 1)
        M = N + num_layers * N^2
        M_bln = round(M / 1_000_000_000, digits = 2)

        lines!(
            ax,
            available_procs,
            perfect_arr;
            linestyle = :dash,
            color = :gray,
            label = perfect_label,
        )

        perfect_label_added = true

        lines!(
            ax,
            available_procs,
            time_arr;
            label = "L=$(num_layers), $(M_bln) bln unknowns",
            color = palette[i]
        )

        scatter!(
            ax,
            available_procs,
            time_arr,
            markersize=13,
            color = :transparent,
            strokecolor = palette[i],
            strokewidth = 2,
        )
    end

    axislegend(ax, position = :lb)

    mkpath("plots")
    save("plots/large_strong_scaling_layers_nx$(nx)_ny$(ny)_reps$(num_reps).pdf", fig)
end

function plot_efficiency_layers(
    nx,
    ny,
    proc_list,
    num_layers_arr,
    num_reps;
    palette=[colorant"#C77A2F",colorant"#2F836C",colorant"#63C7C9", colorant"#A63E33"]
    )
    fig = Figure(fontsize=18)

    ax = Axis(
        fig[1, 1];
        xlabel = "Number of MPI processes",
        ylabel = "Parallel efficiency",
        xticks = proc_list,
        xscale = log2,
        limits = (nothing, (0, 1.4))
    )

    reference_label_added = false

    for (i,num_layers) in enumerate(num_layers_arr)
        available_procs, efficiency_arr, p_ref = get_available_efficiency_data(
            proc_list,
            num_reps,
            nx,
            ny,
            num_layers,
        )

        if isempty(available_procs)
            @warn "No data available" nx ny num_layers
            continue
        end

        reference_label = reference_label_added ? nothing : "Perfect efficiency"

        lines!(
            ax,
            available_procs,
            ones(length(available_procs));
            linestyle = :dash,
            color = :gray,
            label = reference_label,
        )

        reference_label_added = true

        N = (nx + 1) * (ny + 1)
        M = N + num_layers * N^2
        M_bln = round(M / 1_000_000_000, digits = 2)

        label = "L=$(num_layers), $(M_bln) bln unknowns"

        lines!(
            ax,
            available_procs,
            efficiency_arr;
            label = label,
            color = palette[i]
        )

        scatter!(
            ax,
            available_procs,
            efficiency_arr,
            markersize=13,
            color = :transparent,
            strokecolor = palette[i],
            strokewidth = 2,

        )
    end

    axislegend(ax, position = :lb)

    mkpath("plots")
    save("plots/large_efficiency_layers_nx$(nx)_ny$(ny)_reps$(num_reps).pdf", fig)
end

function plot_efficiency_nodes(
    nxny_arr, 
    proc_list, 
    num_layers, 
    num_reps;
    palette=[colorant"#C77A2F",colorant"#2F836C",colorant"#63C7C9", colorant"#A63E33"]
    )
    fig = Figure(fontsize = 18)

    ax = Axis(
        fig[1, 1];
        xlabel = "Number of MPI processes",
        ylabel = "Parallel efficiency",
        xticks = proc_list,
        xscale = log2,
        limits = (nothing, (0, 1.4))
    )

    hlines!(
        ax,
        [1.0];
        linestyle = :dash,
        color = :gray,
        label = "Perfect efficiency",
    )

    for (i,(nx, ny)) in enumerate(nxny_arr)
        available_procs, efficiency_arr, p_ref =
            get_available_efficiency_data(
                proc_list,
                num_reps,
                nx,
                ny,
                num_layers,
            )

        if isempty(available_procs)
            @warn "No data available" nx ny num_layers
            continue
        end

        N = (nx + 1) * (ny + 1)
        M = N + num_layers * N^2
        M_bln = round(M / 1_000_000_000, digits = 2)

        if p_ref == 1
            label = "$N nodes, $(M_bln) bln unknowns"
        else
            label = "$N nodes, $(M_bln) bln unknowns"
        end

        lines!(
            ax,
            available_procs,
            efficiency_arr;
            label = label,
            color = palette[i]
        )

        scatter!(
            ax,
            available_procs,
            efficiency_arr,
            markersize=13,
            color = :transparent,
            strokecolor = palette[i],
            strokewidth = 2,
        )
    end

    axislegend(ax, position = :lb)

    mkpath("plots")

    save("plots/large_efficiency_nodes_num_layers$(num_layers)_reps$(num_reps).pdf",fig)

    return fig
end

function plot_comunication_vs_computation(
    nx,
    ny,
    proc_list,
    num_layers,
    num_reps;
    palette = [colorant"#9EA8B0", colorant"#D7B667"]
    )

    comm_arr_full, comp_arr_full = get_comm_comp_data(proc_list,num_reps,nx,ny,num_layers)

    total = comm_arr_full .+ comp_arr_full
    ratio_comp = comp_arr_full./total
    total_1 = ones(length(total))

    fig = Figure(size = (800,350),fontsize=25)

    x = 1:length(proc_list)

    ax = Axis(
        fig[1, 1],
        xlabel = "MPI processes",
        ylabel = "Fraction",
        xticks = (x, string.(proc_list)),
    )


    barplot!(
        ax,
        x,
        ratio_comp,
        color = palette[1],
        label = "Computation"
    )

    barplot!(
        ax,
        x,
        total_1,
        fillto = ratio_comp,
        color = palette[2],
        label = "Communication"
    )

    axislegend(ax, position = :rb)
    save("plots/communication_vs_computation_nx$(nx)_ny$(ny)_num_layers$(num_layers).pdf",fig)
end

function plot_communication_vs_computation_legend()
    fig = Figure(size = (300, 40))

    Legend(
        fig[1, 1],
        [
            PolyElement(color = colorant"#9EA8B0"),
            PolyElement(color = colorant"#D7B667")
        ],
        ["Computation", "Communication"],
        orientation = :horizontal
    )

    save("plots/communication_vs_computation_legend.pdf", fig)
end


function plot_partial_strong_scaling(
    nx,
    ny,
    proc_list,
    num_layers;
    palette = [colorant"#58636D",colorant"#7FA58F",colorant"#A6ADB0",colorant"#D7B667",colorant"#C87955"],
    )
    fig = Figure(fontsize=18)

    components_labels= [
        "Preparation",
        "rhsIW",
        "Mapping u",
        "Communication",
        "Layers update",
    ]

    ax = Axis(
        fig[1, 1];
        xlabel = "Number of MPI processes",
        ylabel = "Wall clock time (s)",
        xticks = proc_list,
        xscale = log2,
        yscale = log2,
        limits=(nothing,(2^(-8),2^6))
    )

    for time_component in [2,4,5]
        time_arr = get_partial_strong_scaling_data(
            time_component,
            proc_list,
            nx,
            ny,
            num_layers
        )

        lines!(
            ax,
            proc_list,
            time_arr,
            label = components_labels[time_component],
            color = palette[time_component]
        )

        scatter!(ax, proc_list, time_arr, color = palette[time_component])
    end     

    axislegend(ax, position = :rt)

    save("plots/partial_strong_scaling_nx$(nx)_ny$(ny)_num_layers$(num_layers).pdf", fig)
end  

function plot_memory_comparison_per_rank(nxny_array, proc, num_layers)

    n = length(nxny_array)

    node_layers = 2
    layer_nx, layer_ny = (100,100)

    memory_nodes = [
        get_max_memory_across_processors(
            nxny_array[i][1],
            nxny_array[i][2],
            proc,
            node_layers,
            6,
        )
        for i in 1:n
    ]

    memory_layers = [
        get_max_memory_across_processors(
            layer_nx,
            layer_ny,
            proc,
            num_layers[i],
            6,
        )
        for i in 1:n
    ]

    N = [(nx + 1) * (ny + 1) for (nx, ny) in nxny_array]

    M = N .+ node_layers .* N.^2
    M_bln = M/1000000000
    memory_layers_GB = memory_layers/1000000000000
    memory_nodes_GB = memory_nodes/1000000000000
    x = 1:n

    fig = Figure(fontsize=18)

    palette = [colorant"#4682B4",colorant"#008080"]

    ax = Axis(
        fig[1, 1],
        xlabel = "Number of unknowns M (bln)",
        ylabel = "Max memory input memory GB",
        xticks = (x, round.(M_bln;digits=2))
    )

    barplot!(
        ax,
        x .- 0.2,
        memory_nodes_GB,
        width = 0.4,
        color = palette[1],
        label = "Nodes experiment"
    )

    barplot!(
        ax,
        x .+ 0.2,
        memory_layers_GB,
        width = 0.4,
        color = palette[2],
        label = "Layers experiment"
    )
    axislegend(ax, position = :lt)
    save("plots/total_mem_per_rank_M$(M_bln).pdf",fig)
end

function plot_memory_comparison(nxny_array, proc, num_layers)

    n = length(nxny_array)

    node_layers = 2
    layer_nx, layer_ny = (100,100)

    memory_nodes = [
        get_memory_data(
            nxny_array[i][1],
            nxny_array[i][2],
            proc,
            node_layers,
            6,
        )
        for i in 1:n
    ]

    memory_layers = [
        get_memory_data(
            layer_nx,
            layer_ny,
            proc,
            num_layers[i],
            6,
        )
        for i in 1:n
    ]
    summed_memory_layers_TB = sum.(memory_layers)/1000000000000
    summed_memory_nodes_TB = sum.(memory_nodes)/1000000000000

    N = [(nx + 1) * (ny + 1) for (nx, ny) in nxny_array]

    M = N .+ node_layers .* N.^2
    M_bln = M/1000000000
    x = 1:n

    fig = Figure(size = (900,500), fontsize=25)

    palette = [colorant"#4682B4",colorant"#008080"]

    ax = Axis(
        fig[1, 1],
        xlabel = "Number of unknowns M",
        ylabel = "Total input memory (TB)",
        xticks = (x, string.(round.(M_bln;digits=2)," bln"))
    )

    barplot!(
        ax,
        x .- 0.2,
        summed_memory_nodes_TB,
        width = 0.4,
        color = palette[1],
        label = "Nodes experiment"
    )

    barplot!(
        ax,
        x .+ 0.2,
        summed_memory_layers_TB,
        width = 0.4,
        color = palette[2],
        label = "Layers experiment"
    )
    axislegend(ax, position = :lt)
    save("plots/total_mem_per_rank_M$(M_bln).pdf",fig)
end

function plot_first_and_last(layer, filename)
    data = load("data/exp_raw/$filename")
    u_data = data["u"]
    t = data["t"]
    u_last = u_data[end] 
    u_initial = u_data[1]

    num_el = length(u_data[1])
    mesh = range(-30,30,length = num_el)
    fig = Figure(size=(800,400),fontsize=18)
    ax = Axis(
        fig[1,1],
        xlabel = "x", 
        ylabel = "u(x,t)",
        xlabelsize = 25,
        ylabelsize = 25,
        xticks = [-30,-15,0,15,30]
        )
    lines!(ax,mesh,u_last, color =:orange, linewidth = 2, label = "Final time, t=$(round(t[end]; digits = 2))")
    lines!(ax,mesh, u_initial, color=:blue, linewidth = 2, label = "Initial time, t=$(t[1])")
    axislegend(ax, position = :rt)
    save("plots/first_last_layer$layer.pdf", fig)
end 

function common_colorrange(filename)
    data = load("data/exp_raw/$filename")

    u_data = data["u"][end]
    z_data = data["z"][end]

    u_min = minimum(minimum.(u_data))
    u_max = maximum(maximum.(u_data))

    z_min = minimum(minimum.(Iterators.flatten(z_data)))
    z_max = maximum(maximum.(Iterators.flatten(z_data)))

    return (min(u_min, z_min), max(u_max, z_max))
end

function plot_u_evolution_heatmap(layer, filename)
    data = load("data/exp_raw/$filename")
    u_data = data["u"]
    u_heatmap = reduce(hcat,u_data)
    t_heatmap = data["t"]
    num_el = length(u_data[1])
    mesh = range(-30,30,length = num_el)
    colors = common_colorrange(filename)

    fig = Figure(size = (400,800), fontsize=18)
    ax = Axis(fig[1,1], xlabel = "x", ylabel = "t", 
    xticks = [-30,0,30], 
    yticks = [round(t_heatmap[1]; digits=1),round(t_heatmap[end]; digits=2)])
    hm = heatmap!(ax,mesh, t_heatmap,u_heatmap,colorrange = colors)
    Colorbar(fig[1, 2],hm,label = "u(x,t)")
    save("plots/u_over_time_$layer.pdf", fig)
end

function visualize_z(num_layers, filename)

    data = load("data/exp_raw/$filename")
    z_data = data["z"]
    t_data = data["t"]

    z_final = z_data[end]

    colors = common_colorrange(filename)

    fig = Figure(size=(800,400), fontsize = 18)

    heatmaps = []

    for layer in 1:length(z_final)

        ax = Axis(
            fig[1,layer],
        )

        hm = heatmap!(
            ax,
            permutedims(z_final[layer]);
            colorrange = colors
        )

        push!(heatmaps, hm)
    end

    # Colorbar(fig[1, length(z_final) + 1],heatmaps[1])

    save("plots/z_layers_$(num_layers)_final.pdf",fig)
end