using DrWatson
import GalerkinToolkit as GT
using CairoMakie
using LaTeXStrings

@quickactivate "neural_field"

include(srcdir("main_function.jl"))
include(srcdir("mesh.jl"))
include(srcdir("equations.jl"))

# Spatial discretization
mesh = create_cartesian_mesh([-30, 30], (100,))
Ω = GT.interior(mesh)
V = GT.lagrange_space(Ω, 1)

node_x = GT.node_coordinates(V)
x = range(-30, 30, length = 101)

# Initial conditions
initial_u = φ.(node_x)
initial_z = one_layer(node_x)

# Figure
fig = Figure(size = (1300, 520))

# Initial u
ax_u = Axis(
    fig[1, 1],
    xlabel = L"x",
    ylabel = L"u(x, 0)",
    xticks = [-30, -15, 0, 15, 30],
    xlabelsize = 30,
    ylabelsize = 30,
    xticklabelsize = 25,
    yticklabelsize = 25,
    xgridvisible = false,
    ygridvisible = false,
)

lines!(
    ax_u,
    x,
    initial_u,
    color = :darkblue,
    linewidth = 4,
)

xlims!(ax_u, -30, 30)

# Initial z
ax_z = Axis(
    fig[1, 2],
    xlabel = L"x",
    ylabel = L"y",
    xticks = [-30, -15, 0, 15, 30],
    yticks = [-30, -15, 0, 15, 30],
    xlabelsize = 30,
    ylabelsize = 30,
    xticklabelsize = 25,
    yticklabelsize = 25,
    xgridvisible = false,
    ygridvisible = false,
)

hm_z = heatmap!(
    ax_z,
    x,
    x,
    initial_z,
    colormap = :viridis,
)

Colorbar(
    fig[1, 3],
    hm_z,
    label = L"z(x, y, 0)",
    labelsize = 30,
    ticklabelsize = 25,
)

colsize!(fig.layout, 2, Aspect(1, 1))
colgap!(fig.layout, 1, 40)
colgap!(fig.layout, 2, 10)

save("plots/initial_conditions.png", fig, px_per_unit = 3)
save("plots/initial_conditions.pdf", fig)