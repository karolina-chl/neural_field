import JLD2
using DrWatson

@quickactivate "neural_field"

include(srcdir("FEM.jl"))
include(srcdir("mesh.jl"))
include(srcdir("equations.jl"))
include(srcdir("main_function.jl"))

layers = parse(Int,ARGS[1])

make_G_params(layers) = (;
    Ω = (
        build_Ω = create_cartesian_mesh,
        Ω_args = ((-30,30),(300,))
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
        num_layer = layers
        ),
    synaptic_matrix = (
        w = w, 
        synaptic_builder = synaptic_matrix, 
        ),
    other = (
        interpolation_degree = 1, 
        integration_degree = 2
        )
)

S_params = (
    simulation_time = 300, 
    solution_time_step = nothing,
    save_time_step = 1,
    stop_at_equilibrium = true,
    equilibrium_tolerance = 1e-7 # only plays a role when stop_at_equilibrium set to true 
)

params = (
    post_processing = (
        movie = false, 
        movie_timestep = 1
    ),
    save_data = true, 
    save_layers_data = true, 
    datafile_name = "solution_layers_$layers.jld2"
)

println("Executing the code with $layers layers")
G_params = make_G_params(layers)
time = @elapsed solver(G_params, S_params, params)

println("It took $time seconds to finish")




