import GalerkinToolkit as GT
import Gmsh: gmsh
import GLMakie as Makie

function create_cartesian_mesh(domain, num_elements_per_dir)
    mesh = GT.cartesian_mesh(domain,num_elements_per_dir)
    return mesh
end

function circle_mesh(gmsh, face_size, R)
    dim = 2
    gmsh.option.setNumber("General.Verbosity", 2)
    circle_tag = gmsh.model.occ.add_circle(0,0,0,R)
    circle_curve_tag = gmsh.model.occ.add_curve_loop([circle_tag])
    circle_surf_tag = gmsh.model.occ.add_plane_surface([circle_curve_tag])
    gmsh.model.occ.synchronize()
    gmsh.model.model.add_physical_group(dim,[circle_surf_tag],-1,"cortex")
    gmsh.option.setNumber("Mesh.MeshSizeMax",face_size)
    gmsh.model.mesh.generate(dim)
    GT.mesh_from_gmsh(gmsh)
end

function create_circle_mesh(face_size, R) 
    mesh = GT.with_gmsh(gmsh -> circle_mesh(gmsh, face_size, R)) 
    return mesh 
end

function plot_2d_mesh(Ω)
    """Plot 2D mesh"""
    axis = (aspect = Makie.DataAspect(),)
    shading = Makie.NoShading
    fig = GT.makie_surfaces(Ω;color=:pink,axis,shading)
    GT.makie_edges!(Ω;color=:blue)
    Makie.save("plots/domain.png", fig)
    println("Figure saved in the plots folder as 'domain.png'")
end 

