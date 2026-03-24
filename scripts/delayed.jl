using DrWatson

@quickactivate "neural_field"

include(srcdir("FEM.jl"))
include(srcdir("equations.jl"))
include(srcdir("main_function.jl"))

G_params = (;
    Ω = (
        build_Ω = one_d_mesh, 
        Ω_args = (L = 15, num_el = 100)
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
        num_layer = 2,
        delay_function = compute_delay_matrix
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
    simulation_time = 50, 
    solution_time_step = nothing,
    save_time_step = 1
)

L = G_params.Ω.Ω_args.L
num_el = G_params.Ω.Ω_args.num_el
simulation_time = S_params.simulation_time

params = (
    post_processing = (
        movie = false, 
        movie_timestep = 1
    ),
    save_data = true, 
    datafile_name = "1D_data_newinitial_2403_L$(L)_num_el$(num_el)T$(simulation_time)_danielez"
)


print("Executing the code")
main(G_params, S_params, params)

