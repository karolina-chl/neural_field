using DrWatson
using CairoMakie
using JSON
using Statistics

@quickactivate "neural_field"

include(srcdir("plotting_functions.jl"))

# Strong scaling and efficiency - large scenarios 
plot_strong_scaling_nodes([(400,400),(500,500),(600,600)],[128,256,512,1024,2048,4096],2,50)
plot_efficiency_nodes([(400,400),(500,500),(600,600)],[128,256,512,1024,2048,4096],2,50)

plot_strong_scaling_layers(100,100,[128,256,512,1024,2048,4096],[497,1211,2508],50)
plot_efficiency_layers(100,100,[128,256,512,1024,2048,4096],[497,1211,2508],50)

# Strong scaling and efficiency - small scenarios 
plot_strong_scaling_nodes([(100,100),(200,200),(300,300)],[1,2,4,8,16,32,64,128,256,512,1024],2,50)
plot_efficiency_nodes([(100,100),(200,200),(300,300)],[1,2,4,8,16,32,64,128,256,512,1024],2,50)

plot_strong_scaling_layers(100,100,[1,2,4,8,16,32,64,128,256,512,1024],[2,31,158],50)
plot_efficiency_layers(100,100,[1,2,4,8,16,32,64,128,256,512,1024],[2,31,158],50)

# Communication and Computation ratios 

plot_comunication_vs_computation(100,100,[1,2,4,8,16,32,64,128,256,512,1024],2,50)

plot_comunication_vs_computation(100,100,[1,2,4,8,16,32,64,128,256,512,1024],31,50)
plot_comunication_vs_computation(200,200,[2,4,8,16,32,64,128,256,512,1024],2,50)

plot_comunication_vs_computation(100,100,[32,64,128,256,512,1024],158,50)
plot_comunication_vs_computation(300,300,[32,64,128,256,512,1024],2,50)

plot_comunication_vs_computation(100,100,[128,256,512,1024,2048,4096],497,50)
plot_comunication_vs_computation(400,400,[128,256,512,1024,2048,4096],2,50)

plot_comunication_vs_computation(100,100,[512,1024,2048,4096],1211,50)
plot_comunication_vs_computation(500,500,[512,1024,2048,4096],2,50)

plot_comunication_vs_computation(100,100,[1024,2048,4096],2508,50)
plot_comunication_vs_computation(600,600,[1024,2048,4096],2,50)

plot_communication_vs_computation_legend()

# Communication and Computation ratios
plot_partial_strong_scaling(100,100,[1,2,4,8,16,32,64,128,256,512,1024],2)

plot_partial_strong_scaling(100,100,[1,2,4,8,16,32,64,128,256,512,1024],31)
plot_partial_strong_scaling(200,200,[2,4,8,16,32,64,128,256,512,1024],2)

plot_partial_strong_scaling(100,100,[32,64,128,256,512,1024],158)
plot_partial_strong_scaling(300,300,[32,64,128,256,512,1024],2)

plot_partial_strong_scaling(100,100,[128,256,512,1024,2048,4096],497)
plot_partial_strong_scaling(400,400,[128,256,512,1024,2048,4096],2)

plot_partial_strong_scaling(100,100,[512,1024,2048,4096],1211)
plot_partial_strong_scaling(500,500,[512,1024,2048,4096],2)

plot_partial_strong_scaling(100,100,[1024,2048,4096],2508)
plot_partial_strong_scaling(600,600,[1024,2048,4096],2)

# Memory comparison 
plot_memory_comparison(
    [(100,100),(200,200),(300,300),(400,400),(500,500),(600,600)],
    1024, 
    [2,31,158,497,1211,2508]
)