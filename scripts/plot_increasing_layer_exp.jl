using JLD2

# visualize few 

data = load("data/exp_raw/solution_layers_2.jld2")
println(data["u"])