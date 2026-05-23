using JSON
using Plots

function read_rhs_times(json_file)
    data = JSON.parsefile(json_file)

    # Your JSON is a Vector with one entry per debug rank.
    # Each entry has "rhs" as a Vector of RHS runs.
    #
    # Example:
    # data[1]["rhs"] == [[step1, step2, ...]]
    rhs_by_rank = [rank_data["rhs"] for rank_data in data]

    # For now, take the first rank.
    # In your debug output all ranks are identical anyway.
    rhs_runs = rhs_by_rank[1]

    # If there are multiple RHS repetitions, average them step-by-step.
    n_runs = length(rhs_runs)
    n_steps = length(rhs_runs[1])

    rhs_mean = [
        sum(rhs_runs[r][step] for r in 1:n_runs) / n_runs
        for step in 1:n_steps
    ]

    return rhs_mean
end

function plot_rhs_times(json_file, output_file)
    rhs_times = read_rhs_times(json_file)

    step_labels = ["step $i" for i in eachindex(rhs_times)]

    bar(
        step_labels,
        rhs_times,
        xlabel = "RHS step",
        ylabel = "Time spent in each step [s]",
        title = "Time spent in RHS steps",
        legend = false,
        xrotation = 45
    )

    savefig(output_file)

    return rhs_times
end

json_file = "debug.json"
output_file = "debug_rhs_times.png"

rhs_times = plot_rhs_times(json_file, output_file)

println("RHS step times:")
for i in eachindex(rhs_times)
    println("step ", i, ": ", rhs_times[i], " s")
end

println("Saved plot to: ", output_file)