using JSON
using GLMakie

function visualize_first_and_last(datafile, L)
    data = JSON.parsefile(datafile)
    #print(data.u)
    z = data.u[1]
    y = data.u[2]

    x = range(-L,L; length = length(z))
    fig = Figure(size = (700, 450))

    ax = Axis(
        fig[1,1],
        xlabel = "x",
        ylabel = "value", 
        title = "U - final and initial state"
    )
    lines!(ax, x, z, label = "Initial time", linewidth = 2)
    lines!(ax, x, y, label = "Final time", linewidth = 2)
    Legend(fig[1,2],ax)
    fig
end     

visualize_first_and_last("data/exp_raw/1D_data",15)