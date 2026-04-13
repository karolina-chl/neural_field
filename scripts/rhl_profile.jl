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

function setup_rhs_delayed_problem(;
    L = 10,
    num_el = 100,
    num_layer = 2,
)
    interpolation_degree = 1
    integration_degree = 1

    # Build mesh 
    mesh = one_d_mesh(L, num_el)
    Ω = GT.interior(mesh)
    V = GT.lagrange_space(Ω, interpolation_degree)
    node_x = GT.node_coordinates(V)
    dΩ = GT.quadrature(Ω, integration_degree)

    # Build workspace 
    W = synaptic_matrix(V, dΩ, w)
    workspace = (; W, V, dΩ, f)

    # Initial conditions
    node_u = φ.(node_x)
    all_z = z_initial(num_layer, node_x)
    uz = ArrayPartition(node_u, all_z...)
    duz = similar(uz)
    t = 0.0
    p = 0

    return duz, uz, p, t, workspace
end

function run_rhs_delayed(duz, uz, p, t, workspace; n = 10)
    for _ in 1:n
        rhs_delayed!(duz, uz, p, t; workspace)
    end
end


# Setup
duz, uz, p, t, workspace = setup_rhs_delayed_problem(
    L = 10,
    num_el = 100,
    num_layer = 2,
)

# run once before 
run_rhs_delayed(duz, uz, p, t, workspace; n = 1)

#profile 
Profile.clear()
ProfileView.closeall()

ProfileView.@profview begin
    run_rhs_delayed(duz, uz, p, t, workspace; n = 10)
end