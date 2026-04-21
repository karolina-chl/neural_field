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

    if new_fun == true
        I_qp = build_qp_interpolation_matrix(V, dΩ, n_nodes)
        n_qp = size(I_qp, 1)
        z_qp_buf = Vector{Float64}(undef, n_qp)
        workspace = (; W, V, dΩ, f, node_x, I_qp, z_qp_buf)
    else 
        workspace = (;W,V,dΩ,f,node_x, I_face, face_to_dofs, u_face_nodes, u_face_points)
    end         
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
            #rhs_delayed_new!(duz, uz, p, t; workspace)
            rhs_delayed_GT_test!(duz, uz, p, t; workspace)
        end 
    end       
end



########################################################################################
### Experiment 1: Check the implementations and compare 
# duz_1, uz_1, p, t, workspace_1 = setup_rhs_delayed_problem(;
#     L = 10,
#     num_el = 100,
#     num_layer = 2,
#     new_fun = true
# )

# run_rhs_delayed(duz_1, uz_1, p, t, workspace_1; n=1, new_fun = true)

# duz_2, uz_2, p, t, workspace_2 = setup_rhs_delayed_problem(;
#     L = 10,
#     num_el = 100,
#     num_layer = 2,
#     new_fun = false
# )

# run_rhs_delayed(duz_2, uz_2, p, t, workspace_2; n=1, new_fun = false)

# @show isapprox(duz_1.x[1], duz_2.x[1]; rtol=1e-12, atol=1e-12)

# for k in 2:length(duz_1.x)
#     @show k isapprox(duz_1.x[k], duz_2.x[k]; rtol=1e-12, atol=1e-12)
# end

# @show uz_1 == uz_2

###############################################################################################
### Experiment 2: Profiling the code

duz, uz, p, t, workspace = setup_rhs_delayed_problem(;
    L = 10,
    num_el = 100,
    num_layer = 2,
    new_fun = false
) 

@time run_rhs_delayed(duz, uz, p, t, workspace; n = 1, new_fun = false)

Profile.clear()
ProfileView.closeall()

ProfileView.@profview begin
    run_rhs_delayed(duz, uz, p, t, workspace; n = 10, new_fun = false)
end





















##### Discarted 
# function check_allocation(duz, uz, p, t, workspace)
#     node_node_last_z = uz.x[3]

#     node_z = view(node_node_last_z, :, 1) # put i the workspace
#     disc_Ω = GT.discrete_field(workspace.V, node_z)
#     z_faces = GT.each_face(disc_Ω, workspace.dΩ; tabulate = (GT.value,))

#     node_z = view(node_node_last_z, :, 1)
#     z_faces = GT.replace_free_values(z_faces,node_z)
  
# end 

# function one_layer_update(duz, uz, workspace)
#     num_layer = length(uz.x) - 1
#     node_u = uz.x[1]
#     node_node_dz1 = duz.x[2]
#     node_node_z1 = uz.x[2]
#     n_nodes = length(node_u)

#     for i in 1:n_nodes
#         for j in 1:n_nodes
#             node_i = workspace.node_x[i]
#             node_j = workspace.node_x[j]
#             node_node_dz1[j,i] = α(node_i, node_j, num_layer) * (node_u[j] - node_node_z1[j,i])
#         end
#     end
#     return nothing
# end

# function one_rhs_u_update!(node_du, node_u, node_node_last_z, workspace, i)
#     node_z1 = view(node_node_last_z, :, 1)
#     disc_Ω = GT.discrete_field(workspace.V, node_z1)
#     z_faces = GT.each_face(disc_Ω, workspace.dΩ; tabulate = (GT.value,))

#     node_zi = view(node_node_last_z, :, i)
#     z_faces = GT.replace_free_values(z_faces, node_zi)

#     qp = 0
#     node_du[i] = 0.0
#     for z_face in z_faces
#         for z_point in GT.each_point(z_face)
#             qp += 1
#             point_node_z = GT.field(GT.value, z_point)
#             fz = workspace.f(point_node_z)
#             node_du[i] += workspace.W[i, qp] * fz
#         end
#     end
#     node_du[i] -= node_u[i]
#     return nothing
# end

#@btime one_layer_update!(duz, uz, workspace)
#@btime check_allocation(duz, uz, p, t, workspace)


# num_layer = length(uz.x) - 1
# node_u = uz.x[1]
# node_du = duz.x[1]
# node_node_last_z = uz.x[num_layer + 1]
#@btime one_rhs_u_update!($node_du, $node_u, $node_node_last_z, $workspace, 1)