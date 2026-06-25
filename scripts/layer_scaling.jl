using DrWatson

num_reps = 50
nx=ny=19
layers_list = [2,32,512]
comm_arr = []
comp_arr = []

num_proc = 8

for layer in layers_list
    comm, comp = get_max_comm_comp_time_over_repetitions(num_proc,num_reps,nx,ny,layer)
    push!(comm_arr, comm)
    push!(comp_arr, comp)

    println(median(comm))
end 

fig = Figure()

ax = Axis(
    fig[1,1], 
    xticks = (1:3, ["2 layers", "32 layers", "512 layers"]), 
    limits = (nothing, (0,0.005)),
)

for i in eachindex(comm_arr)
    boxplot!(ax, fill(i, length(comm_arr[i])), comm_arr[i])
end

save("plots/layers_boxplot_num_proc$(num_proc)with_barrier.png",fig)