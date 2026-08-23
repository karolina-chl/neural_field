using DrWatson
using JLD2
using LinearAlgebra
using CairoMakie
import GalerkinToolkit as GT

@quickactivate "neural_field"

include(srcdir("equations.jl"))
include(srcdir("mesh.jl"))

# plot norm 

function plot_norm_increasing_layers(layers_arr, layer_reff)
    data = load("data/exp_raw/solution_layers_$layer_reff.jld2")
    u_data_reff = data["u"]
    u_norm_reff = LinearAlgebra.norm(u_data_reff, Inf)
    norm_arr = []
    diff_layer_arr = [] # temp

    for layer in layers_arr
        data_layer = load("data/exp_raw/solution_layers_$layer.jld2")
        u_data_layer = data_layer["u"]
        diff_layer = u_data_reff.-u_data_layer
        norm_layer = LinearAlgebra.norm(diff_layer, Inf)
        norm_normalize = norm_layer/u_norm_reff
        push!(norm_arr, norm_normalize)
        push!(diff_layer_arr, diff_layer) # temp
    end

    fig=Figure(size = (800,400), fontsize=18)
    ax = Axis(
        fig[1,1],
        xlabel = "Number of layers",
        ylabel = "Relative infinity-norm error",
        xticks = layers_arr,
        yscale = log2,
        )

    lines!(ax,layers_arr,norm_arr, color= :gray, linewidth = 2)
    scatter!(ax,layers_arr,norm_arr, color= :gray, markersize = 10)

    save("plots/norm_increasing_layers_inf_norm.pdf",fig)

    return diff_layer_arr # temp
end     

# plot first and last and heatmap

function plot_first_and_last(layer)
    data = load("data/exp_raw/solution_layers_$layer.jld2")
    u_data = data["u"]
    t = data["t"]
    last = length(t)
    u_last = u_data[last]
    u_initial = u_data[1]

    mesh = range(-30,30,length = 501)
    fig = Figure(fontsize=18)
    ax = Axis(
        fig[1,1],
        xlabel = "Mesh", 
        ylabel = "U value"
        )
    lines!(ax,mesh,u_last, color =:orange, label = "Final U")
    lines!(ax,mesh, u_initial, color=:blue, label = "Initial U")
    axislegend(ax, position = :rt)
    save("plots/first_last_layer$layer.pdf", fig)
end 

function plot_u_evolution_heatmap(layer)
    data = load("data/exp_raw/solution_layers_$layer.jld2")
    u_data = data["u"]
    u_heatmap = reduce(hcat,u_data)
    t_heatmap = data["t"]
    mesh = range(-30,30,length = 501)

    heatmap = Figure()
    ax = Axis(heatmap[1,1], xlabel = "U", ylabel = "Time")
    heatmap!(ax,mesh, t_heatmap,u_heatmap)
    save("plots/u_over_time_$layer.png", heatmap)
end  

## Example 
# layers_arr = [2,4,6,8,10,12,14,16,18]
# diff_layer_arr = plot_norm_increasing_layers(layers_arr,20)

# fig = Figure()
# ax = Axis(fig[1,1])
# data = hcat(diff_layer_arr[11]...)
# hm = heatmap!(data)
# Colorbar(fig[1,2],hm)
# fig

plot_first_and_last(7)
plot_u_evolution_heatmap(7)


#plot_u_evolution_heatmap("equilibrium_model_layers_20.jld2")

