using DrWatson
import GalerkinToolkit as GT
import GLMakie as Makie
import DifferentialEquations as DE
import ProgressMeter as PM
import JSON
using LinearAlgebra
using SparseArrays
using RecursiveArrayTools


function main(G_params, S_params, params)

    ### Create a mesh
    mesh_size = G_params.Ω.Ω_args.mesh_size
    R = G_params.Ω.Ω_args.R
    build_Ω = G_params.Ω.build_Ω

    axis = (aspect = Makie.DataAspect(),)
    colormap=:viridis
    mesh = GT.with_gmsh(gmsh -> build_Ω(gmsh, mesh_size, R))
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
    dim_u = length(node_u)
    all_z = initialize_z(num_layer, dim_u)

    ### Numerical integration 
    integration_degree = G_params.other.integration_degree
    dΩ = GT.quadrature(Ω,integration_degree)

    ### Synaptic matrix
    synaptic_builder = G_params.synaptic_matrix.synaptic_builder
    w = G_params.synaptic_matrix.w

    W = synaptic_builder(V,dΩ,w) 
    
    #ODE right-hand-side
    f = G_params.firing_function.f

    n_nodes = length(node_u)
    n_points = size(W,2) 
    node_wfz = similar(node_x, Float64) 
    point_node_fz = zeros(n_points, n_nodes) # use similar like before?
    workspace = (;W,V,dΩ,node_wfz,point_node_fz, f)

    #ODE solution
    T = S_params.simulation_time
    uz_array = ArrayPartition(node_u, all_z...) # three dots to treat matrix separate
    ode = DE.ODEProblem(uz_array,[0,T]) do args...
        rhs_delayed!(args...;workspace)
    end

    ode_solution = DE.solve(ode,DE.Tsit5();saveat=0.01) #saves every 0.01 time units

    #save the data
    if params.save_data
        open("data/solution.json", "w") do io
            JSON.json(io, ode_solution)
        end
    end 

    #plot the solution and save as mp3
    if params.post_processing.movie
        color = GT.discrete_field(V,node_u)
        fig = Makie.Figure()
        _,sc = GT.makie_surfaces(fig[1,1],Ω;color,axis,refinement=3,colormap)
        fn = "solution.mp4"
        integrator = DE.init(ode, DE.Tsit5())
        prog = PM.ProgressThresh(0.0)
        Makie.record(fig,fn,DE.tuples(integrator);framerate=10) do (uz,t)
            node_u = uz.x[1] 
            sc.color = GT.discrete_field(V,node_u)
            PM.update!(prog,T-t)
        end
    end
end