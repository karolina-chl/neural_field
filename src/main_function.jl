using DrWatson
import GalerkinToolkit as GT
#using Makie
using GLMakie
import DifferentialEquations as DE
import ProgressMeter as PM
import JSON
using LinearAlgebra
using SparseArrays
using RecursiveArrayTools



function main(G_params, S_params, params)

    ### Create a mesh 2D
    # mesh_size = G_params.Ω.Ω_args.mesh_size
    # R = G_params.Ω.Ω_args.R
    # build_Ω = G_params.Ω.build_Ω

    # mesh = GT.with_gmsh(gmsh -> build_Ω(gmsh, mesh_size, R))
    # Ω = GT.interior(mesh)

    L = G_params.Ω.Ω_args.L
    num_el = G_params.Ω.Ω_args.num_el
    
    mesh = one_d_mesh(L, num_el)
    Ω = GT.interior(mesh)

    ### Finite element interpolation
    interpolation_degree = G_params.other.interpolation_degree
    V = GT.lagrange_space(Ω,interpolation_degree)
    node_x = GT.node_coordinates(V)

    ### Initial conditions
    initialize_u = G_params.state_initialization.initialize_u
    initialize_z = G_params.state_initialization.initialize_z
    num_layer = G_params.state_initialization.num_layer

    node_u = initialize_u.(node_x)
    all_z = initialize_z(num_layer, node_x)
    
    ### Numerical integration 
    integration_degree = G_params.other.integration_degree
    dΩ = GT.quadrature(Ω,integration_degree)

    ### Synaptic matrix
    synaptic_builder = G_params.synaptic_matrix.synaptic_builder
    w = G_params.synaptic_matrix.w

    W = synaptic_builder(V,dΩ,w) 
    
    #ODE right-hand-side
    f = G_params.firing_function.f

    #new
    V_faces= GT.each_face(V,dΩ;tabulate=(GT.value,))
    I_face = transpose(V_faces.accessor.reference_space_face.workspace.values[1])
    face_to_dofs = GT.face_dofs(V)

    #new vectores
    u_face_nodes = zeros(size(I_face,2))
    u_face_points = zeros(size(I_face,1))
    workspace = (;W,V,dΩ,f,node_x, I_face, face_to_dofs, u_face_nodes, u_face_points)
    
    #ODE solution
    T = S_params.simulation_time
    tspan = (0.0, float(T))
    uz_array = ArrayPartition(node_u, all_z...) # three dots to treat matrix separate
    ode = DE.ODEProblem(uz_array,tspan) do args...
        rhs_delayed!(args...;workspace)
    end
    
    save_dt = S_params.save_time_step

    if isnothing(save_dt)
        ode_solution = DE.solve(ode,DE.Tsit5())
    else
        ode_solution = DE.solve(ode,DE.Tsit5();saveat=save_dt)
    end

    #save the data
    state_uz = ode_solution.u
    timesteps = ode_solution.t
    file_name = params.datafile_name

    ode_solution_extracted = Dict(
        "t" => timesteps, 
        "u" => [state_uz[t].x[1] for t in 1:length(timesteps)],
        "z" => [state_uz[t].x[2:num_layer+1] for t in 1:length(timesteps)] # do I need this?
    )

    if params.save_data
        open("data/exp_raw/$file_name", "w") do io
            JSON.json(io, ode_solution_extracted) 
        end
    end 

    println("Solution saved to JSON")
end