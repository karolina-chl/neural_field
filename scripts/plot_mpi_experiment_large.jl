using DrWatson
using GLMakie
using JSON
using Statistics

@quickactivate "neural_field"

include(srcdir("parallel_delayed.jl"))
include(srcdir("data_processing.jl"))

function result_file(nx, ny, np, num_layers)
    return datadir(
        "exp_raw",
        "mpi_exp",
        "results_nx$(nx)ny$(ny)np$(np)num_layers$(num_layers).json",
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

    # Prefer true strong-scaling efficiency relative to one processor.
    if isfile(result_file(nx, ny, 1, num_layers))
        p_ref = 1
        t_ref = get_t1_time(nx, ny, num_layers)
    else
        # Otherwise compute relative efficiency from the first available proc count.
        p_ref = available_procs[1]
        t_ref = get_strong_scaling_data([p_ref], num_reps, nx, ny, num_layers)[1]
    end

    # Reuse your existing data_processing function.
    # Passing p_ref * t_ref gives:
    #   eff = (p_ref * t_ref) / (p * t_p)
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

function plot_strong_scaling_nodes(nxny_arr,proc_list,num_layers,num_reps)
    fig = Figure()

    ax = Axis(
        fig[1, 1];
        xlabel = "Number of processors",
        ylabel = "Wall clock time (s)",
        xticks = proc_list,
        xscale = log2,
        yscale = log2,
    )

    perfect_label_added = false

    for (nx, ny) in nxny_arr
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

        # Perfect scaling starts from the first available processor count.
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
            label = "$(M_bln) bln unknowns",
        )

        scatter!(
            ax,
            available_procs,
            time_arr,
        )
    end

    axislegend(ax, position = :lb)

    mkpath("plots")
    save("plots/large_strong_scaling_nodes.png", fig)
end

function plot_strong_scaling_layers(nx,ny,proc_list,num_layers_arr,num_reps)

    fig = Figure(size=(1200,600))

    ax = Axis(
        fig[1, 1];
        xlabel = "Number of processors",
        ylabel = "Wall clock time (s)",
        xticks = proc_list,
        xscale = log2,
        yscale = log2,
    )

    perfect_label_added = false

    for num_layers in num_layers_arr
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

        # Perfect scaling starts from the first available processor count.
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
        )

        scatter!(
            ax,
            available_procs,
            time_arr,
        )
    end

    Legend(fig[1,2],ax)

    mkpath("plots")
    save("plots/large_strong_scaling_layers_nx$(nx)_ny$(ny)_reps$(num_reps).png", fig)
end

function plot_efficiency_layers(nx,ny,proc_list,num_layers_arr,num_reps)
    fig = Figure(size=(1200,600))

    ax = Axis(
        fig[1, 1];
        xlabel = "Number of processors",
        ylabel = "Parallel efficiency",
        xticks = proc_list,
        xscale = log2,
        limits = (nothing, (0, 1.6))
    )

    reference_label_added = false

    for num_layers in num_layers_arr
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
        )

        scatter!(
            ax,
            available_procs,
            efficiency_arr,
        )
    end

    #axislegend(ax, position = :lt)
    Legend(fig[1,2],ax)

    mkpath("plots")
    save("plots/large_efficiency_layers_nx$(nx)_ny$(ny)_reps$(num_reps).png", fig)
end

# Usage: plotting 
nxny_arr = [(160,160),(320,160),(320,320),(640,320),(640,640)]
num_reps = 10
num_layers = 2
proc_list = [128,256,512,1024,2048]
nx=ny=130
num_layers_arr = [16,32,64,128,256,512]

#plot_strong_scaling_nodes(nxny_arr,proc_list,num_layers,num_reps)

plot_strong_scaling_layers(nx,ny,proc_list,num_layers_arr,num_reps)
plot_efficiency_layers(nx,ny,proc_list,num_layers_arr,num_reps)