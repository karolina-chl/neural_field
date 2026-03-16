using DrWatson
import GalerkinToolkit as GT

@quickactivate "neural_field"

include(srcdir("FEM.jl"))
include(srcdir("equations.jl"))
include(srcdir("main_function.jl"))

G_params = (;
    Ω = (
        build_Ω = circle_mesh, 
        Ω_args = (mesh_size = 3, R = 30)
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

make_S_params(solution_time_step) = (
    simulation_time = 100, 
    solution_time_step = solution_time_step,
    save_time_step = [100] # if you want to save only last timestep T, just insert [T]
)

make_params(datafile_name) = ( 
    post_processing = (
        movie = false, 
        movie_timestep = 1
    ),
    save_data = true,
    datafile_name = datafile_name
)

### Create a mesh
mesh_size = G_params.Ω.Ω_args.mesh_size
R = G_params.Ω.Ω_args.R
build_Ω = G_params.Ω.build_Ω

mesh = GT.with_gmsh(gmsh -> build_Ω(gmsh, mesh_size, R))
Ω = GT.interior(mesh)

### Finite element interpolation
interpolation_degree = G_params.other.interpolation_degree
V = GT.lagrange_space(Ω,interpolation_degree)

### Numerical integration 
integration_degree = G_params.other.integration_degree
dΩ = GT.quadrature(Ω,integration_degree)

timesteps = (0.1, 0.05, 0.01)
results = Dict()

for t in timesteps
    params = make_params("data_time_$t")
    main(G_params, S_params, params)

    V, dΩ = integral_mesh(G_params)
    data = JSON.parsefile("data/exp_raw/data_time_$t")
    u_final_state = Float64.(data.u[1])

    uh = GT.solution_field(V,u_final_state)
    integral = GT.∫(x -> abs(uh(x)),dΩ) |> sum 
    results["timestep: $t"] = integral

    open("data/exp_raw/final-results_time_$t", "w") do io
        JSON.json(io, results)
    end
end