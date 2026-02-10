"""
Functions needed for FEM method 
"""

import GalerkinToolkit as GT
import GLMakie as Makie
import Gmsh

include("equations.jl")

function circle_mesh(gmsh, mesh_size, R)
    dim = 2
    gmsh.option.setNumber("General.Verbosity", 2)
    circle_tag = gmsh.model.occ.add_circle(0,0,0,R)
    circle_curve_tag = gmsh.model.occ.add_curve_loop([circle_tag])
    circle_surf_tag = gmsh.model.occ.add_plane_surface([circle_curve_tag])
    gmsh.model.occ.synchronize()
    gmsh.model.model.add_physical_group(dim,[circle_surf_tag],-1,"cortex")
    gmsh.option.setNumber("Mesh.MeshSizeMax",mesh_size)
    gmsh.model.mesh.generate(dim)
    GT.mesh_from_gmsh(gmsh)
end

function plot_circle_mesh(Ω)
    """Plot circle mesh"""
    axis = (aspect = Makie.DataAspect(),)
    shading = Makie.NoShading
    fig = GT.makie_surfaces(Ω;color=:pink,axis,shading)
    GT.makie_edges!(Ω;color=:blue)
    Makie.save("domain.png", fig)
end 

function synaptic_matrix(V,dΩ)
    node_x = GT.node_coordinates(V)
    # Count number of nz in matrix
    nnz = 0
    for dΩ_face in GT.each_face(dΩ)
        for dΩ_point in GT.each_point(dΩ_face)
            y = GT.coordinate(dΩ_point)
            for (node,x) in enumerate(node_x)
                wxy = w(x,y)
                if wxy != 0
                    nnz += 1
                end
            end
        end
    end
    # Allocate COO arrays
    WI = zeros(Int32,nnz)
    WJ = zeros(Int32,nnz)
    WV = zeros(Float64,nnz)
    # Fill COO arrays
    nnz = 0
    point = 0
    for dΩ_face in GT.each_face(dΩ)
        for dΩ_point in GT.each_point(dΩ_face)
            point += 1
            y =  GT.coordinate(dΩ_point)
            dy = GT.weight(dΩ_point)
            for (node,x) in enumerate(node_x)
                wxy = w(x,y)
                if wxy != 0
                    nnz += 1                    
                    WI[nnz] = node
                    WJ[nnz] = point
                    WV[nnz] = wxy*dy
                end
            end
        end
        
    end
    npoints = point
    nnodes = length(node_x)
    # Compress into CSC format
    W = sparse(WI,WJ,WV,nnodes,npoints)
    return W
end


function rhs!(node_du,node_u,p,t;workspace)
    (;W,V,dΩ,node_wfu,point_fu) = workspace
    u = GT.discrete_field(V,node_u)
    u_faces = GT.each_face(u,dΩ;tabulate=(GT.value,))
    point = 0
    for u_face in u_faces
        for u_point in GT.each_point(u_face)
            point += 1
            u = GT.field(GT.value,u_point)
            point_fu[point] = f(u)
        end
    end
    mul!(node_wfu,W,point_fu)
    node_du .= node_wfu .- node_u
end

function rhs_delayed!(
    (node_du, node_node_dz),  # L = 1
    (node_u, node_node_z),
    p,
    t;
    workspace::NamedTuple,
)
    (; W, shape_space, dΩ, node_wfz, point_fz, f) = workspace
    n_nodes = length(node_u)
    n_points = size(W,2)
    # Interpolation of the z_L from nodes to quadratures
    for i in 1:n_nodes
        node_z = view(node_node_z, :, i)
        disc_Ω = GT.discrete_field(shape_space, node_z)
        z_faces = GT.each_face(disc_Ω, dΩ; tabulate = (GT.value,))
        p = 0  # μ is going to denote the quadrature points in the field (in the documentation it replaces ũ)
        for z_face in z_faces
            for z_point in GT.each_point(z_face)
                p += 1
                z = GT.field(GT.value, z_point)
                point_node_fz[p, i] = f(z) #consider removing the matrix point_node_fz
            end
        end
        node_du[i] = 0
        for p in 1:n_points
            node_du[i] +=  W[i,p]*point_node_fz[p, i]
        end 
        node_du[i] -= node_u[i] 
    end
    ## Q part is missing llop over j and i and both loops are lenths 

    for i in 1:n_nodes
        for j in 1:n_nodes
            node_node_dz[j,i] = alpha_ji *(node_u[j] - node_node_z[j,i]) #need to compute alpha (separate function)
        end     
    end    
end
