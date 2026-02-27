using DrWatson

@quickactivate "neural_field"

include(srcdir("FEM.jl"))
include(srcdir("equations.jl"))
include(srcdir("main_function.jl"))

G_params = (;
    Ω = (
        build_Ω = circle_mesh, 
        Ω_args = (mesh_size = 7, R = 30)
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
        num_layer = 4
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
    simulation_time = 400, 
    time_step = 1 # not used yet
)

params = ( # not used yet
    post_processing = (
        movie = true, 
        movie_timestep = 1
    ),
    save_data = false
)

main(G_params, S_params, params)
