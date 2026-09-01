using DrWatson
import GalerkinToolkit as GT
import PartitionedArrays as PA
import JSON
using LinearAlgebra
using GLMakie
using RecursiveArrayTools: ArrayPartition
using Statistics

@quickactivate "neural_field"
include(srcdir("parallel_delayed.jl"))

np = 1

# # generate the data 
num_elements = [16,32,64,128,256,512,1024,2048]
mesh_sizes = [
    (4, 4),
    (4, 8),
    (8, 8),
    (8, 16),
    (16, 16),
    (16, 32),
    (32, 32),
    (32, 64)
]

# for (nx,ny) in mesh_sizes
#     main_debug(nx,ny,1,2)
# end 


#extract the data - memory
mem_array = Float64[] 

for (idx, value) in enumerate(mesh_sizes)
    nx, ny = value 
    time_data = JSON.parsefile("data/exp_raw/parallel_time_10x/results_nx$(nx)ny$(ny)np$(np).json")
    push!(mem_array, time_data[1]["mem"][1][1])
end

# plotting 

mem_MB = round.(mem_array./1000000, digits = 2)

fig = Figure()
ax = Axis(fig[1,1];
    xlabel = "Number of elements",
    ylabel = "Memory of setup (MB)", 
    xscale = log2, 
    yscale = log2,
    xticks = (num_elements),
    yticks = (mem_MB),
    title = "Memory of setup"
)
lines!(ax, num_elements, mem_MB, color = :purple)
scatter!(ax, num_elements, mem_MB, color = :purple)
save("plots/memory_scaling.png",fig)

# extract the data - time 
time_arr = Float64[]
for value in mesh_sizes[2:8]
    nx, ny = value
    exp_data = JSON.parsefile("data/exp_raw/parallel_time_10x/results_nx$(nx)ny$(ny)np$(np).json")
    rhs_time = exp_data[1]["rhs"]
    summed = sum.(rhs_time)

    mean_rhs_time = mean(summed)
    push!(time_arr, mean_rhs_time)
end     

# plotting - times

time_arr_round = round.(time_arr, digits = 4)

fig = Figure()
ax = Axis(fig[1,1];
    xlabel = "Number of elements",
    ylabel = "Time (s)", 
    xscale = log2, 
    yscale = log2, 
    xticks = (num_elements[2:8]),
    yticks = time_arr_round,
    title = "Mean time of 1x rhs execution"
)
lines!(ax, num_elements[2:8], time_arr_round, color = :purple)
scatter!(ax, num_elements[2:8], time_arr_round, color = :purple)
save("plots/time_scaling.png",fig)
