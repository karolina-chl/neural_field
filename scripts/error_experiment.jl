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
    solution_time_step = 1,
    save_time_step = [400] # if you want to save only last timestep T, just insert [T]
)

make_params(datafile_name) = ( 
    post_processing = (
        movie = false, 
        movie_timestep = 1
    ),
    save_data = true,
    datafile_name = datafile_name
)

function integral_mesh(G_params)
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
    return V, dΩ
end

mesh_size_testset = (10,9,8,7,6,5,4,3,2,1,1/2,1/4)
results = Dict()

for size in mesh_size_testset
    G_params = make_G_params(size)
    params = make_params("data_$size")
    main(G_params, S_params, params)

    V, dΩ = integral_mesh(G_params)
    data = JSON.parsefile("data/exp_raw/data_$size")
    u_final_state = Float64.(data.u[1])

    uh = GT.solution_field(V,u_final_state)
    integral = GT.∫(x -> abs(uh(x)),dΩ) |> sum 
    results["mesh_size: $size"] = integral
end 


# save the file with final results 
    open("data/exp_raw/final-results", "w") do io
        JSON.json(io, results; allownan=true)
    end


