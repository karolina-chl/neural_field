using JSON
using GLMakie
using LinearAlgebra

function visualize_first_and_last(datafile, L)
    """Visualizes U at the first timestep and at the last timestep"""
    data = JSON.parsefile(datafile)
  
    z = data.u[1]
    y = data.u[51]

    x = range(-L,L; length = length(z))
    fig = Figure(size = (700, 450))

    ax = Axis(
        fig[1,1],
        yscale = log10,
        xlabel = "Mesh",
        ylabel = "U value", 
        title = "U - final and initial state"
    )
    lines!(ax, x, z, label = "Initial time", linewidth = 2)
    lines!(ax, x, y, label = "Final time", linewidth = 2)
    Legend(fig[1,2],ax)
    save("first_and_last.png", fig)
end


function visualize_U_history(datafile)
    data = JSON.parsefile(datafile)
    u = data["u"]             
    u_matrix = hcat(u...) 

    fig = Figure()
    ax = Axis(fig[1, 1], xlabel = "U mesh", ylabel = "Time")
    hm = GLMakie.heatmap!(ax, u_matrix; colorscale = log10, colorrange = (1e-10, 5))
    Colorbar(fig[1, 2], hm)
    save("U_history.png", fig)
end

function visualize_z(datafile)
    data = JSON.parsefile(datafile)
    z = data.z[51]
    z_matrix = hcat(z[2]...)

    fig = Figure()
    ax = Axis(fig[1,1])
    hm = GLMakie.heatmap!(ax, z_matrix)
    Colorbar(fig[1,2], hm)
    save("z_final.png",fig)
end 

#visualize_first_and_last("data/exp_raw/1D_data_newinitial_2403_L15_num_el100T50_danielez", 15)
visualize_U_history("data/exp_raw/1D_data_newinitial_2403_L15_num_el100T50_danielez")
#visualize_z("data/exp_raw/1D_data_newinitial_2403_L15_num_el100T50_danielez")


