"""
Functions needed for FEM method 
"""

import GalerkinToolkit as GT
#import GLMakie as Makie
using RecursiveArrayTools
import Gmsh

include("equations.jl")

function one_d_mesh(L, num_el) 
    return GT.cartesian_mesh((-L, L),(num_el,))
end    

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

function synaptic_matrix(V,dΩ,w)
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
    duz,
    uz,
    p,
    t;
    workspace::NamedTuple,
)
    (; W, V, dΩ, f) = workspace

    num_layer = length(uz.x) - 1

    node_u  = uz.x[1] # modifies in place 
    node_du = duz.x[1]  

    n_nodes = length(node_u)
    node_node_last_z = uz.x[num_layer + 1]

    # Interpolation of the z_L from nodes to quadratures
    for i in 1:n_nodes
        node_z = view(node_node_last_z, :, i)
        disc_Ω = GT.discrete_field(V, node_z)
        z_faces = GT.each_face(disc_Ω, dΩ; tabulate = (GT.value,))
        qp = 0  
        node_du[i] = 0
        for z_face in z_faces
            for z_point in GT.each_point(z_face)
                qp += 1
                point_node_z = GT.field(GT.value, z_point) #this is an element from point_node z matrix
                fz = f(point_node_z)
                node_du[i] +=  W[i,qp]*fz
            end
        end
        node_du[i] -= node_u[i]
    end

    node_x = GT.node_coordinates(V)
    node_node_dz1 = duz.x[2]
    node_node_z1 = uz.x[2]

    #for the first layer 
    for i in 1:n_nodes #iterating column by column 
        for j in 1:n_nodes
            node_i = node_x[i]
            node_j = node_x[j]
            node_node_dz1[j,i] = α(node_i, node_j, num_layer)*(node_u[j] - node_node_z1[j,i]) 
        end
    end
    
    # for all other layers
    for layer in 2:num_layer
        node_node_dzn = duz.x[layer+1] #indeks switched by 1 because of the u
        node_node_z1 = uz.x[layer]
        node_node_z2 = uz.x[layer+1]
        for i in 1:n_nodes #iterating column by column 
            for j in 1:n_nodes
                node_i = node_x[i]
                node_j = node_x[j]
                node_node_dzn[j,i] = α(node_i, node_j, num_layer) *(node_node_z1[j,i] - node_node_z2[j,i]) 
            end     
        end
    end  
end


