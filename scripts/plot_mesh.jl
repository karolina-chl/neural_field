using DrWatson
import GalerkinToolkit as GT
using GLMakie

include(srcdir("mesh.jl"))

function create_mesh(nx, ny, interpolation_degree, integration_degree)
    mesh = GT.cartesian_mesh((0, 1, 0, 1), (nx, ny))
    Ω = GT.interior(mesh)
    dΩ = GT.quadrature(Ω, integration_degree)
    V = GT.lagrange_space(Ω, interpolation_degree)

    V_faces = GT.each_face(V, dΩ, tabulate = (GT.value, GT.gradient))
    fq_fn_I = transpose(V_faces.accessor.reference_space_face.workspace.values[1])

    gf_fn_gn = GT.face_nodes(V)
    gn_x = GT.node_coordinates(V)

    nfq = size(fq_fn_I, 1)
    ngf = length(gf_fn_gn)
    ngq = nfq * ngf
    ngn = length(gn_x)

    println("interpolation degree = $interpolation_degree, integration degree = $integration_degree")
    println("Number of global faces is $ngf")
    println("Number of global nodes is $ngn")
    println("Number of global quadrature points is $ngq")

    return V, Ω, dΩ
end

function quadrature_coordinates(V, dΩ)
    V_faces = GT.each_face(V, dΩ, tabulate = (GT.value, GT.gradient))
    fq_fn_I = transpose(V_faces.accessor.reference_space_face.workspace.values[1])

    gf_fn_gn = GT.face_nodes(V)
    gn_x = GT.node_coordinates(V)

    q_x = eltype(gn_x)[]

    for face in eachindex(gf_fn_gn)
        face_nodes = gf_fn_gn[face]

        for q in axes(fq_fn_I, 1)
            xq = zero(gn_x[1])

            for local_node in eachindex(face_nodes)
                global_node = face_nodes[local_node]
                xq += fq_fn_I[q, local_node] * gn_x[global_node]
            end

            push!(q_x, xq)
        end
    end

    return q_x
end

function plot_mesh(V,Ω,qn_x)
    fig = Figure()
    ax = Axis(fig[1, 1], aspect = Makie.DataAspect())
    node_x = GT.node_coordinates(V)

    GT.makie_edges!(
        ax,
        Ω;
        color = :gray,
        linewidth = 0.5,
        label = "Element edges",
    )
    
    scatter!(
        ax,
        node_x;
        color = :steelblue,
        markersize = 10,
        label = "Nodes",
    )

    scatter!(
        ax,
        qn_x;
        color = :darkorange,
        markersize = 11,
        marker = :xcross,
        label = "Quadrature points",
    )

    Legend(fig[1,2], ax)
    save("plots/mesh.png",fig)

    println("Plotted mesh saved in the folder plots as mesh.png")
end 

function plot_2d_mesh(Ω)
    """Plot 2D mesh"""
    axis = (
        aspect = Makie.DataAspect(),
        xticks=[0,1], 
        yticks=[0,1],
        xticklabelsize = 30,
        yticklabelsize = 30)
    shading = Makie.NoShading
    fig = GT.makie_surfaces(Ω;color=:white,axis,shading)
    GT.makie_edges!(Ω;color=:purple,linewidth = 2)
    save("plots/domain.png", fig;px_per_unit = 6)
    save("plots/domain.pdf", fig;px_per_unit = 6)
    println("Figure saved in the plots folder as 'domain.png'")
end 


###### Create a mesh and save plotted mesh in the plots folder 

V, Ω, dΩ = create_mesh(5, 5, 1, 1)
qn_x = quadrature_coordinates(V, dΩ)
plot_mesh(V, Ω, qn_x)

mesh = create_cartesian_mesh((0,1,0,1), (3,3))
plot_2d_mesh(mesh)