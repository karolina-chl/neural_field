using JSON
using Statistics

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
    slowest_worker = []
    for num in 1:num_reps
        max_time = maximum(summed_arr[rank][num] for rank in 1:num_proc)
        push!(slowest_worker, max_time)
    end 
    return slowest_worker
end  

function get_strong_scaling_data(proc_list, num_reps, nx, ny, num_layers)
    best_time_arr = []
    for proc in proc_list
        max_arr = get_max_rhs_time_over_repetitions(proc, num_reps, nx, ny, num_layers)
        push!(best_time_arr, mean(max_arr))
    end     
    return best_time_arr
end 

function get_efficiency_data(proc_list, t1, num_reps, nx, ny, num_layers)
    efficiency = []
    for proc in proc_list
        max_arr = get_max_rhs_time_over_repetitions(proc, num_reps, nx, ny, num_layers)
        tp = mean(max_arr)
        eff = t1/(proc*tp)
        push!(efficiency, eff)
    end 
    return efficiency     
end

function get_max_comm_comp_time_over_repetitions(num_proc,num_reps,nx,ny,num_layers)
    communication_arr = Float64[]
    computation_arr = Float64[]
    data = JSON.parsefile("data/exp_raw/mpi_exp/results_nx$(nx)ny$(ny)np$(num_proc)num_layers$(num_layers).json")
    rhs_times = [entry["rhs"] for entry in data]
    
    summed_arr = [sum.(rhs_times[rank]) for rank in 1:num_proc]

    for rep in 1:num_reps
        slowest_idx = argmax(summed_arr[rank][rep] for rank in 1:num_proc)
        communication = rhs_times[slowest_idx][rep][4]
        computation = summed_arr[slowest_idx][rep]-communication
        push!(communication_arr, communication)
        push!(computation_arr, computation)
    end 
    return communication_arr, computation_arr
end  

function get_comm_comp_data(proc_list,num_reps,nx,ny,num_layers)
    communication_arr_full = Float64[]
    computation_arr_full = Float64[]
    for proc in proc_list
        communication_arr, computation_arr = get_max_comm_comp_time_over_repetitions(proc,num_reps,nx,ny,num_layers)
        push!(communication_arr_full, mean(communication_arr))
        push!(computation_arr_full, mean(computation_arr))
    end   
    return communication_arr_full, computation_arr_full
end     

function get_t1_time(nx, ny, num_layers)
    data_t1 = JSON.parsefile("data/exp_raw/mpi_exp/results_nx$(nx)ny$(ny)np1num_layers$(num_layers).json")
    t1 = minimum(sum.(data_t1[1]["rhs"]))
    return t1       
end 