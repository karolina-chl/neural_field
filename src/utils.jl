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
include(srcdir("parallel_delayed.jl"))


########################################
# Comparing GT vs non-GT implementation 
########################################


function setup_rhs_delayed_problem(;
    L = 10,
    num_el = 100,
    num_layer = 2,
    new_fun = true
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
    n_nodes = length(node_x)

    #new
    V_faces= GT.each_face(V,dΩ;tabulate=(GT.value,))
    I_face = transpose(V_faces.accessor.reference_space_face.workspace.values[1])
    face_to_dofs = GT.face_dofs(V)

    #new vectores
    u_face_nodes = zeros(size(I_face,2))
    u_face_points = zeros(size(I_face,1))
    workspace = (;W,V,dΩ,f,node_x, I_face, face_to_dofs, u_face_nodes, u_face_points)

    # Initial conditions
    node_u = φ.(node_x)
    all_z = z_initial(num_layer, node_x)
    uz = ArrayPartition(node_u, all_z...)
    duz = similar(uz)
    t = 0.0
    p = 0

    return duz, uz, p, t, workspace
end

function run_rhs_delayed(duz, uz, p, t, workspace; n = 10, new_fun = true)
    if new_fun == false
        for _ in 1:n
            rhs_delayed!(duz, uz, p, t; workspace)
        end
    else
        for _ in 1:n
            rhs_delayed_corrected!(duz, uz, p, t; workspace)
        end 
    end       
end

########################################
# Testing parallel implementation
########################################

function materialize_debug_pvector(u)
    parts = PA.local_values(u)
    return reduce(vcat, parts.items)
end

function materialize_debug_layer(z_layer)
    return reduce(hcat, z_layer.items)
end

function materialize_debug_duz(duz)
    du = materialize_debug_pvector(duz.x[1])

    dz_layers = [
        materialize_debug_layer(duz.x[i])
        for i in 2:length(duz.x)
    ]

    return ArrayPartition(du, dz_layers...)
end

function setup_for_parallel_test(nx, ny, num_layer)
    mesh = GT.cartesian_mesh((0,1,0,1),(nx,ny))
    Ω = GT.interior(mesh)
    dΩ = GT.quadrature(Ω,2)
    V = GT.lagrange_space(Ω,2)
    node_x = GT.node_coordinates(V)

    # Build workspace 
    W = synaptic_matrix(V, dΩ, w)
    #n_nodes = length(node_x)

    #new
    V_faces= GT.each_face(V,dΩ;tabulate=(GT.value,))
    I_face = transpose(V_faces.accessor.reference_space_face.workspace.values[1])
    face_to_dofs = GT.face_dofs(V)

    #new vectores
    u_face_nodes = zeros(size(I_face,2))
    u_face_points = zeros(size(I_face,1))
    workspace = (;W,V,dΩ,f,node_x, I_face, face_to_dofs, u_face_nodes, u_face_points)

    # Initial conditions
    node_u = φ.(node_x)
    all_z = z_initial(num_layer, node_x)
    uz = ArrayPartition(node_u, all_z...)
    duz = similar(uz)
    t = 0.0
    p = 0
    return duz, uz, p, t, workspace
end

function test_parallel_implementation(duz, uz, p, t, workspace;nx,ny,np,num_layers, parallel = false)
    if parallel == false
        rhs_delayed_corrected!(duz, uz, p, t; workspace)
    else
        main_debug(nx,ny,np,num_layers)
    end       
end