using JSON
using GLMakie
using LinearAlgebra

include(srcdir("FEM.jl"))

function visualize_first_and_last(datafile, L)
    """Visualizes U at the first timestep and at the last timestep"""
    data = JSON.parsefile(datafile)
  
    z = data.u[1]
    y = data.u[51]

    x = range(-L,L; length = length(z))
    fig = Figure(size = (700, 450))

    ax = Axis(
        fig[1,1],
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
    hm = GLMakie.heatmap!(ax, u_matrix)
    Colorbar(fig[1, 2], hm)
    save("U_history.png", fig)
end

function visualize_z(datafile)
    data = JSON.parsefile(datafile)
    z = data.z[51]
    z_matrix_1 = hcat(z[1]...)
    z_matrix_2 = hcat(z[2]...)

    fig_1 = Figure()
    ax_1 = Axis(fig_1[1,1])
    hm_1 = GLMakie.heatmap!(ax_1, z_matrix_1)
    Colorbar(fig_1[1,2], hm_1)
    save("z1_final.png",fig_1)

    fig_2 = Figure()
    ax_2 = Axis(fig_2[1,1])
    hm_2 = GLMakie.heatmap!(ax_2, z_matrix_2)
    Colorbar(fig_2[1,2], hm_2)
    save("z2_final.png",fig_2)

end 

#visualize_first_and_last("data/exp_raw/Test_test_test_23.04", 30)
#visualize_U_history("data/exp_raw/Test_test_test_23.04")
#visualize_z("data/exp_raw/Test_test_test_23.04")





