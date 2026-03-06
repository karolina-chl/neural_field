using DrWatson
import GalerkinToolkit as GT

@quickactivate "neural_field"

include(srcdir("FEM.jl"))
include(srcdir("equations.jl"))
include(srcdir("main_function.jl"))

make_G_params(mesh_size) = (;
    Ω = (
        build_Ω = circle_mesh, 
        Ω_args = (mesh_size = mesh_size, R = 30)
        ),
    firing_function = (;
        f = f
        ),
    forcing_function = (;
        g = nothing
        ),
    state_initialization = (
        initialize_u = φ,
        initialize_z = z_initial,
        num_layer = 1
        ),
    synaptic_matrix = (
        w = w, 
        synaptic_builder = synaptic_matrix, 
        ),
    other = (
        interpolation_degree = 1, 
        integration_degree = 1
        )
    )    

S_params = (
    simulation_time = 10, 
    solution_time_step = 5,
    save_time_step = [10] # if you want to save only last timestep T, just insert [T]
)

make_params(datafile_name) = ( 
    post_processing = (
        movie = false, 
        movie_timestep = 1
    ),
    save_data = true,
    datafile_name = datafile_name
)    

mesh_size_testset = (7,6,5)

for size in mesh_size_testset
    G_params = make_G_params(size)
    params = make_params("data_$size")
    main(G_params, S_params, params)
end 

### Create a mesh
G_params = make_G_params(6) # must be adjutsed each time 
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

### Numerical integration 
integration_degree = G_params.other.integration_degree
dΩ = GT.quadrature(Ω,integration_degree)

# read the file 
data = JSON.parsefile("data/exp_raw/data_6") # must be adjusted 
u_final_state = Float64.(data.u[1])

uh = GT.solution_field(V,u_final_state)
integral = GT.∫(x -> abs(uh(x)),dΩ) |> sum 
