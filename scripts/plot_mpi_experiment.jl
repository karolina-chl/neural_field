using DrWatson
using CairoMakie
using JSON
using Statistics

@quickactivate "neural_field"

include(srcdir("plotting_functions.jl"))


#nxny_arr = [(100,100),(200,200),(300,300)]
nxny_arr = [(400,400),(500,500),(600,600)]
num_reps = 50
num_layers = 2
#proc_list = [1,2,4,8,16,32,64]
proc_list = [128,256,512,1024,2048,4096]
nx=ny=100
num_layers_arr = [497,1211,2508]

plot_strong_scaling_nodes(nxny_arr,proc_list,num_layers,num_reps)
plot_efficiency_nodes(nxny_arr,proc_list,num_layers,num_reps)

plot_strong_scaling_layers(nx,ny,proc_list,num_layers_arr,num_reps)
plot_efficiency_layers(nx,ny,proc_list,num_layers_arr,num_reps)

plot_comunication_vs_computation(100,100,proc_list,497,50)
plot_comunication_vs_computation(400,400,proc_list,2,50)

plot_partial_strong_scaling(100,100,proc_list,497)
plot_partial_strong_scaling(400,400,proc_list,2)

plot_memory_comparison(
    [(100,100),(200,200),(300,300),(400,400),(500,500),(600,600)],
    1024, 
    [2,31,158,497,1211,2508]
)