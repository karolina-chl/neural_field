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

function plot_strong_scaling_nodes()
    proc_list = [1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096]

    nxny_arr = [
        (160,160),
        (320,160),
        (320,320),
        (640,320),
        (640,640)
    ]

    num_reps = 10
    num_layers = 2

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

function plot_strong_scaling_layers()
    proc_list = [1,2,4,8,16,32,64,128,256,512]

    # Choose the fixed mesh size for the layer-scaling experiment.
    nx = 70
    ny = 70

    # Choose the layer counts that you actually ran.
    num_layers_arr = [16,32,64,128,256,512,1024,2048]

    num_reps = 10

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

    axislegend(ax, position = :lb)

    mkpath("plots")
    save("plots/large_strong_scaling_layers.png", fig)

    return fig
end

# plot_strong_scaling_nodes()
plot_strong_scaling_layers()