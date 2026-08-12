using DrWatson
import GalerkinToolkit as GT
import DifferentialEquations as DE
import JSON
using JLD2
using LinearAlgebra
using SparseArrays
using RecursiveArrayTools

function mesh_setup_sequencial(mesh, interpolation_degree, integration_degree)
    Ω = GT.interior(mesh)
    V = GT.lagrange_space(Ω,interpolation_degree)
    node_x = GT.node_coordinates(V)
    dΩ = GT.quadrature(Ω,integration_degree)
    V_faces= GT.each_face(V,dΩ;tabulate=(GT.value,))
    I_face = transpose(V_faces.accessor.reference_space_face.workspace.values[1])
    face_to_dofs = GT.face_dofs(V)
    u_face_nodes = zeros(size(I_face,2))
    u_face_points = zeros(size(I_face,1))
    workspace = (;V,dΩ,node_x, I_face, face_to_dofs, u_face_nodes, u_face_points)
    return workspace
end

function save_solution_json(ode_solution,num_layer,file_name; save_layers_data = false)
    state_uz = ode_solution.u
    timesteps = ode_solution.t
    
    if save_layers_data
        ode_solution_extracted = Dict(
            "t" => timesteps, 
            "u" => [state_uz[t].x[1] for t in 1:length(timesteps)],
            "z" => [state_uz[t].x[2:num_layer+1] for t in 1:length(timesteps)]
        )
    else 
        ode_solution_extracted = Dict(
            "t" => timesteps, 
            "u" => [state_uz[t].x[1] for t in 1:length(timesteps)]
        ) 
    end     

    open("data/exp_raw/$file_name", "w") do io
        JSON.json(io, ode_solution_extracted) 
    end
end  

function save_solution_jld(ode_solution,num_layer,file_name; save_layers_data = false)
    state_uz = ode_solution.u
    timesteps = ode_solution.t
    
    u_data = [state_uz[t].x[1] for t in 1:length(timesteps)]

    if save_layers_data
        z_data = [state_uz[t].x[2:num_layer+1] for t in 1:length(timesteps)]
        jldsave("data/exp_raw/$file_name"; t = timesteps, u = u_data, z = z_data)
    else 
        jldsave("data/exp_raw/$file_name"; t = timesteps, u = u_data)
    end     
end

function solver(G_params, S_params, params)

    ### Mesh setup
    build_mesh = G_params.Ω.build_Ω
    mesh_args = G_params.Ω.Ω_args
    interpolation_degree = G_params.other.interpolation_degree
    integration_degree = G_params.other.integration_degree
    arg1,arg2 = mesh_args
    
    mesh = build_mesh(arg1,arg2)
    workspace = mesh_setup_sequencial(mesh, interpolation_degree, integration_degree)

    ### Initial conditions
    initialize_u = G_params.state_initialization.initialize_u
    initialize_z = G_params.state_initialization.initialize_z
    num_layer = G_params.state_initialization.num_layer

    node_u = initialize_u.(workspace.node_x)
    all_z = initialize_z(num_layer, workspace.node_x)


    ### Synaptic matrix and firing function 
    synaptic_builder = G_params.synaptic_matrix.synaptic_builder
    w = G_params.synaptic_matrix.w

    W = synaptic_builder(workspace.V,workspace.dΩ,w) 
    f = G_params.firing_function.f
    
    workspace = (;workspace...,W,f)

    # #ODE solution
    T = S_params.simulation_time
    uz_array = ArrayPartition(node_u, all_z...) # three dots to treat matrix separately
    ode = DE.ODEProblem(uz_array,[0,T]) do args...
        rhs_delayed_corrected!(args...;workspace)
    end
    
    # determine additional arguments for the solver 
    solve_args = (;)

    save_dt = S_params.save_time_step
    if save_dt !== nothing 
        solve_args = merge(solve_args,(;saveat=save_dt))
    end 

    if_equlibrium = S_params.stop_at_equilibrium
    if if_equlibrium
        tol = S_params.equilibrium_tolerance
        termination_condition = DE.TerminateSteadyState(tol, tol)
        solve_args = merge(solve_args, (;callback = termination_condition))
        println("Equilibrium stoping conditions tolerance set to $tol")
    end 
    
    
    #solve the problem 
    ode_solution = DE.solve(ode,DE.Tsit5(); abstol = 1e-10, reltol = 1e-10, solve_args...)

    println("Solver finished. Solver returned: $(ode_solution.retcode)")
    println("Simulated until t= $(ode_solution.t[end])")

    #save the data
    if params.save_data
        file_name = params.datafile_name
        if if_equlibrium
            file_name = params.datafile_name
            file_name = "equilibrium_tol_$(tol)_" * file_name
        end     
        save_layers_data = params.save_layers_data
        save_solution_jld(ode_solution,num_layer,file_name;save_layers_data)
        println("Solution saved in data folder as $file_name")
    end 
end