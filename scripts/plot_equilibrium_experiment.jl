using DrWatson
using CairoMakie
using JSON
using Statistics

@quickactivate "neural_field"

include(srcdir("plotting_functions.jl"))

# Plotting 
# plot_first_and_last(2, "equilibrium_tol_1.0e-7_solution_layers_2.jld2")
# plot_u_evolution_heatmap(2, "equilibrium_tol_1.0e-7_solution_layers_2.jld2")
# visualize_z(2, "equilibrium_tol_1.0e-7_solution_layers_2.jld2")

# Infinity norm 
reference_data = load("data/exp_raw/equilibrium_tol_1.0e-7_solution_layers_2.jld2")
u_reference = reference_data["u"][end]

for tol in ["0.01","0.001", "0.0001", "1.0e-5", "1.0e-6"]
    data = load("data/exp_raw/equilibrium_tol_$(tol)_solution_layers_2.jld2")
    u_last = data["u"][end]
    inf_norm = norm(u_reference .- u_last, Inf)
    println("Tolerance = $tol, infinity norm = $inf_norm")
end