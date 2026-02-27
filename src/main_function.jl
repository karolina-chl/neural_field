using DrWatson
import GalerkinToolkit as GT
import GLMakie as Makie
import DifferentialEquations as DE
import ProgressMeter as PM
using LinearAlgebra
using SparseArrays
using RecursiveArrayTools

function main()

    ### Create a mesh
    mesh_size = 7
    R = 30
    axis = (aspect = Makie.DataAspect(),)
    colormap=:viridis
    mesh = GT.with_gmsh(gmsh -> circle_mesh(gmsh, mesh_size, R))
    Ω = GT.interior(mesh)

    ### Plot a mesh - optional
    plot_circle_mesh(Ω)

    ### Finite element interpolation
    interpolation_degree = 1
    V = GT.lagrange_space(Ω,interpolation_degree)
    node_x = GT.node_coordinates(V)

    ### Initial conditions 
    node_u = φ.(node_x)
    dim_u = length(node_u)
    num_layer = 4
    all_z = z_initial(num_layer, dim_u)

    ### Numerical integration 
    integration_degree = 2*interpolation_degree
    dΩ = GT.quadrature(Ω,integration_degree)
    face_lpoint_x = GT.sample(x->x,dΩ)
    point_x = face_lpoint_x.data
    npoints = length(point_x)

    ### Synaptic matrix
    W = synaptic_matrix(V,dΩ)

    #Finite element function accessors # is there something to be change here?
    u = GT.discrete_field(V,node_u)
    u_faces = GT.each_face(u,dΩ;tabulate=(GT.value,))


    #ODE right-hand-side 
    n_nodes = length(node_u)
    n_points = size(W,2) 
    node_wfz = similar(node_x, Float64) 
    point_node_fz = zeros(n_points, n_nodes) #use similar like before?
    workspace = (;W,V,dΩ,node_wfz,point_node_fz, f)

    #ODE solution
    T = 10 # Use 400 for nicer results
    uz_array = ArrayPartition(node_u, all_z...) #three dots to treat matrix separate
    ode = DE.ODEProblem(uz_array,[0,T]) do args...
        rhs_delayed!(args...;workspace)
    end

    #plot the solution and save as mp3
    color = u
    fig = Makie.Figure()
    ax,sc = GT.makie_surfaces(fig[1,1],Ω;color,axis,refinement=3,colormap)
    fn = "solution.mp4"
    integrator = DE.init(ode, DE.Tsit5())
    prog = PM.ProgressThresh(0.0)
    Makie.record(fig,fn,DE.tuples(integrator);framerate=10) do (node_u,t)
        node_u = uz_array.x[1] 
        sc.color = GT.discrete_field(V,node_u)
        PM.update!(prog,T-t)
    end

end