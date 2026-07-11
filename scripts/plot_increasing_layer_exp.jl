using DrWatson
using JLD2
using LinearAlgebra
using GLMakie
import GalerkinToolkit as GT

@quickactivate "neural_field"

include(srcdir("equations.jl"))
include(srcdir("mesh.jl"))

# plot norm 

function plot_norm_increasing_layers(layers_arr, layer_reff)
    data = load("data/exp_raw/solution_layers_$layer_reff.jld2")
    u_data_reff = data["u"]
    norm_arr = []

    for layer in layers_arr
        data_layer = load("data/exp_raw/solution_layers_$layer.jld2")
        u_data_layer = data_layer["u"]
        diff_layer = u_data_reff.-u_data_layer
        norm_layer = LinearAlgebra.norm(diff_layer)
        push!(norm_arr, norm_layer)
    end

    fig=Figure()
    ax = Axis(
        fig[1,1],
        xlabel = "Number of layers",
        ylabel = "Norm",
        xticks = layers_arr
        )

    lines!(ax,layers_arr,norm_arr, color= :purple)
    scatter!(ax,layers_arr,norm_arr, color= :purple)

    save("plots/layer_behaviour/norm_increasing_layers.png",fig)
end     

# plot first and last and heatmap

function plot_first_and_last(layer)
    data = load("data/exp_raw/solution_layers_$layer.jld2")
    u_data = data["u"]
    u_last = u_data[20]
    u_initial = u_data[1]

    mesh = range(-30,30,length = 501)
    fig = Figure()
    ax = Axis(
        fig[1,1],
        xlabel = "Mesh", 
        ylabel = "U value"
        )
    lines!(ax,mesh,u_last, color =:yellow)
    lines!(ax,mesh, u_initial, color=:blue)
    save("plots/layer_behaviour/first_last_layer$layer.png", fig)
end 

function plot_u_evolution_heatmap(layer)
    data = load("data/exp_raw/solution_layers_$layer.jld2")
    u_data = data["u"]
    u_heatmap = reduce(hcat,u_data)
    t_heatmap = data["t"]

    heatmap = Figure()
    ax = Axis(heatmap[1,1], xlabel = "U", ylabel = "Time")
    heatmap!(ax,mesh, t_heatmap,u_heatmap)
    save("plots/layer_behaviour/u_over_time_$layer.png", heatmap)
end  

## Example 
layers_arr = [2,3,4,5,6,7,8,9,10,11,12,13,14]
plot_norm_increasing_layers(layers_arr,14)