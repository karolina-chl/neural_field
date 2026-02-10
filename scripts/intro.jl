using DrWatson
import GalerkinToolkit as GT
import GLMakie as Makie
import DifferentialEquations as DE
import ProgressMeter as PM
using LinearAlgebra
using SparseArrays

@quickactivate "neural_field"

# Here you may include files from the source directory
include(srcdir("FEM.jl"))
include(srcdir("equations.jl"))

### Create a mesh
mesh_size = 7
R = 30
axis = (aspect = Makie.DataAspect(),)
colormap=:viridis
mesh = GT.with_gmsh(gmsh -> circle_mesh(gmsh, mesh_size, R))
Ω = GT.interior(mesh)
plot_circle_mesh(Ω)

### Finite element interpolation
interpolation_degree = 1
V = GT.lagrange_space(Ω,interpolation_degree)
node_x = GT.node_coordinates(V)

#initial condition
node_u = φ.(node_x)

#numerical integration 
integration_degree = 2*interpolation_degree
dΩ = GT.quadrature(Ω,integration_degree)
face_lpoint_x = GT.sample(x->x,dΩ)
point_x = face_lpoint_x.data
npoints = length(point_x)

#Synaptic matrix 

W = synaptic_matrix(V,dΩ)

#Finite element function accessors 
u = GT.discrete_field(V,node_u)
u_faces = GT.each_face(u,dΩ;tabulate=(GT.value,))
face = 4
u_face = u_faces[face]
u_points = GT.each_point(u_face)
lpoint = 3
u_point = u_points[lpoint]
GT.field(GT.value,u_point)

#ODE right-hand-side 
node_wfu = similar(node_x,Float64)
point_fu = similar(node_wfu,npoints)
workspace = (;W,V,dΩ,node_wfu,point_fu)

#ODE solution
T = 400 # Use 400 for nicer results
ode = DE.ODEProblem(node_u,[0,T]) do args...
    rhs!(args...;workspace)
end

#plot the solution and save as mp3
color = u
fig = Makie.Figure()
ax,sc = GT.makie_surfaces(fig[1,1],Ω;color,axis,refinement=3,colormap)
fn = "solution.mp4"
integrator = DE.init(ode, DE.Tsit5())
prog = PM.ProgressThresh(0.0)
Makie.record(fig,fn,DE.tuples(integrator);framerate=10) do (node_u,t)
    sc.color = GT.discrete_field(V,node_u)
    PM.update!(prog,T-t)
end