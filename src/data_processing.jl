

function get_max_rhs_time_over_repetitions(num_proc, num_reps, nx, ny, num_layers)
    """ 
    For each MPI rank, the individual RHS timing components are summed for every
    repetition. Then, for each repetition, the maximum summed RHS time over all
    MPI ranks is computed. The function finally returns an array of max times. 
    """

    data = JSON.parsefile("data/exp_raw/mpi_exp/results_nx$(nx)ny$(ny)np$(num_proc)num_layers$(num_layers).json")
    rhs_times = [entry["rhs"] for entry in data]
    # get summed times (all rhs parts)
    summed_arr = []
    for num in 1:num_proc
        summed = sum.(rhs_times[num])
        push!(summed_arr, summed)
    end 
    max_arr = []
    for num in 1:num_reps
        max_time = maximum(summed_arr[rank][num] for rank in 1:num_proc)
        push!(max_arr, max_time)
    end 
    return max_arr
end  

function get_strong_scaling_data(proc_list, num_reps, nx, ny, num_layers)
    best_time_arr = []
    for proc in proc_list
        max_arr = get_max_rhs_time_over_repetitions(proc, num_reps, nx, ny, num_layers)
        push!(best_time_arr, minimum(max_arr))
    end     
    return best_time_arr
end     

proc_list = [2,4,8,16,32,64]
nx = 30
ny = 30 
num_reps = 10
num_layers = 2
best_time_arr = get_strong_scaling_data(proc_list, num_reps, nx, ny, num_layers)
print(best_time_arr)