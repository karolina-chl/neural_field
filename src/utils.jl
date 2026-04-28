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