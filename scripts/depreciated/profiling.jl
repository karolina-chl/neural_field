using DrWatson
using BenchmarkTools
using Profile
using ProfileView
using RecursiveArrayTools: ArrayPartition
import GalerkinToolkit as GT

@quickactivate "neural_field"

include(srcdir("FEM.jl"))
include(srcdir("equations.jl"))
include(srcdir("main_function.jl"))
include(srcdir("utils.jl"))

### Experiment: Profiling the code

duz, uz, p, t, workspace = setup_rhs_delayed_problem(;
    L = 10,
    num_el = 100,
    num_layer = 2,
    new_fun = true
) 

@time run_rhs_delayed(duz, uz, p, t, workspace; n = 1, new_fun = false)

Profile.clear()
ProfileView.closeall()

ProfileView.@profview begin
    run_rhs_delayed(duz, uz, p, t, workspace; n = 10, new_fun = false)
end