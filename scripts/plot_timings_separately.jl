using DrWatson
using GLMakie

include(srcdir("data_processing.jl"))

@quickactivate "neural_field"

function plot_partial_strong_scaling_layers(time_component,proc_list,nx,ny,num_layers_arr)
    time_arr = []
    for layer in num_layers_arr
        data = get_partial_strong_scaling_data(
            time_component,
            proc_list,
            nx,
            ny,
            layer
        )

        push!(time_arr, data)
    end

    fig = Figure(size=(1400,800))

    ax = Axis(
        fig[1, 1];
        xlabel = "Number of processors",
        ylabel = "Wall clock time (s)",
        xticks = proc_list,
        xscale = log2,
        yscale = log2,
        xlabelsize = 30,
        ylabelsize = 30,
        xticklabelsize = 30,
        yticklabelsize = 30
    )

    for entry in eachindex(num_layers_arr)
        layer = num_layers_arr[entry]
        layer_times = time_arr[entry]

        perfect_arr = layer_times[1] .* proc_list[1] ./ proc_list

        lines!(
            ax,
            proc_list,
            layer_times,
            label = "$(layer) layers"
        )

        scatter!(ax, proc_list, layer_times)

        perfect_label =
            entry == firstindex(num_layers_arr) ? "Perfect scaling" : nothing

        lines!(
            ax,
            proc_list,
            perfect_arr,
            linestyle = :dash,
            color = :gray,
            label = perfect_label
        )
    end

    Legend(fig[1,2],ax; labelsize = 30)

    save("plots/partial_strong_scaling_layer_$(time_component)time_component.png", fig)
end

proc_list = [1,2,4,8,16,32,64]
nx = ny = 60
num_layers_arr = [16,32,64,128,256]
time_component = 1

plot_partial_strong_scaling_layers(time_component,proc_list,nx,ny,num_layers_arr)