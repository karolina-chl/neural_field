import GalerkinToolkit as GT
using GLMakie

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

function cartesian_mesh_edges(nx, ny)
    xs = Float64[]
    ys = Float64[]

    xgrid = range(0, 1, length = nx + 1)
    ygrid = range(0, 1, length = ny + 1)

    for x in xgrid
        push!(xs, x)
        push!(ys, 0.0)
        push!(xs, x)
        push!(ys, 1.0)
        push!(xs, NaN)
        push!(ys, NaN)
    end

    for y in ygrid
        push!(xs, 0.0)
        push!(ys, y)
        push!(xs, 1.0)
        push!(ys, y)
        push!(xs, NaN)
        push!(ys, NaN)
    end

    return xs, ys
end

function plot_mesh_panel!(ax, nx, ny, interpolation_degree, integration_degree)
    V, Ω, dΩ = create_mesh(nx, ny, interpolation_degree, integration_degree)

    node_x = GT.node_coordinates(V)
    quad_x = quadrature_coordinates(V, dΩ)
    edge_x, edge_y = cartesian_mesh_edges(nx, ny)

    lines!(
        ax,
        edge_x,
        edge_y,
        color = :gray,
        linewidth = 0.8,
        label = "Mesh edges"
    )

    scatter!(
        ax,
        node_x,
        color = :steelblue,
        markersize = 10,
        label = "Nodes"
    )

    scatter!(
        ax,
        quad_x,
        color = :darkorange,
        markersize = 11,
        marker = :xcross,
        label = "Quadrature points"
    )

    return ax
end

function plot_two_meshes_with_shared_legend(nx, ny)
    fig = Figure(size = (1200, 500))

    ax1 = Axis(
        fig[1, 1],
        #title = "V = GT.lagrange_space(Ω, 1), dΩ = GT.quadrature(Ω, 1)",
        aspect = DataAspect()
    )

    ax2 = Axis(
        fig[1, 3],
        #title = "V = GT.lagrange_space(Ω, 2), dΩ = GT.quadrature(Ω, 2)",
        aspect = DataAspect()
    )

    plot_mesh_panel!(ax1, nx, ny, 1, 1)
    plot_mesh_panel!(ax2, nx, ny, 2, 2)

    Legend(fig[1, 2], ax1)

    colsize!(fig.layout, 1, Relative(0.42))
    colsize!(fig.layout, 2, Relative(0.16))
    colsize!(fig.layout, 3, Relative(0.42))

    fig
end

fig = plot_two_meshes_with_shared_legend(5, 5)
save("plots/mesh.png",fig)