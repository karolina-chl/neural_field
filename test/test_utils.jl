using DrWatson
using BenchmarkTools
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

function main_debug_test(nx,ny,np,num_layers)
    PA.with_debug() do backend
        return main_test(backend,np,nx,ny,num_layers)
    end
end

function main_test(backend,np,nx,ny,num_layers; save = false, file_name = title)
    ranks = backend(1:np)

    # Setup
    elap_setup = zeros(1)
    
    elap_setup[1] = @elapsed p_setup = map(ranks) do rank
        prepare_setup_on_rank(rank, np, nx, ny, num_layers)
    end

    mem = Base.summarysize(p_setup)

    # Initial condition
    ngn = PA.getany(map(setup->setup.ngn,p_setup))
    gn_partition = PA.uniform_partition(ranks, np, ngn)
    p_ln_u0  = map(setup_ln_u0, p_setup)
    
    u0 = PA.PVector(p_ln_u0, gn_partition)
    u = similar(u0)
    du = similar(u0)
    
    z_layers = [map(setup_ln_z0, p_setup, p_ln_u0) for _ in 1:num_layers]
    dz_layers = [map(z -> zero(z), z_layers[layer]) for layer in 1:num_layers]
    
    uz = ArrayPartition(u, z_layers...)
    duz = ArrayPartition(du, dz_layers...)

    nr = 1
    r_elap_rhs = [zeros(8) for _ in 1:nr]
    for r in 1:nr

        # reset the repetition
        elap_rhs = r_elap_rhs[r]
        copy!(u,u0)

        rhs!(duz,uz,p_setup,elap_rhs)
    end

    elap = Dict{Symbol,Vector{Vector{Float64}}}()
    elap[:setup] = [elap_setup]
    elap[:rhs] = r_elap_rhs
    elap[:mem] = [[mem]]
    p_elap_main = PA.gather(map(_->elap,ranks))
    if save == true
        file_title = file_name(nx, ny, np, num_layers)
        PA.map_main(p_elap_main) do p_elap
            open("data/exp_raw/mpi_exp/$file_title.json", "w") do io
                JSON.print(io, p_elap)
            end
        end
        println("Experiment finished.Results saved as $file_title")
    end
    return duz     
end

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
    V = GT.lagrange_space(Ω,1)
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
        main_debug_test(nx,ny,np,num_layers)
    end       
end