using DrWatson
import GalerkinToolkit as GT
import GLMakie as Makie
import DifferentialEquations as DE
import ProgressMeter as PM
using LinearAlgebra
using SparseArrays

@quickactivate "neural_field"

include(srcdir("FEM.jl"))
include(srcdir("equations.jl"))

### Create a mesh
mesh_size = 7
R = 30
axis = (aspect = Makie.DataAspect(),)
colormap=:viridis
mesh = GT.with_gmsh(gmsh -> circle_mesh(gmsh, mesh_size, R))
Ω = GT.interior(mesh)

### Plot a mesh - optional
plot_circle_mesh(Ω)

### Finite element interpolation
interpolation_degree = 1
V = GT.lagrange_space(Ω,interpolation_degree)
node_x = GT.node_coordinates(V)

### Initial conditions 
node_u = φ.(node_x)
#node_z = 

### Numerical integration 
integration_degree = 2*interpolation_degree
dΩ = GT.quadrature(Ω,integration_degree)
face_lpoint_x = GT.sample(x->x,dΩ)
point_x = face_lpoint_x.data
npoints = length(point_x)

### Synaptic matrix
W = synaptic_matrix(V,dΩ)

#Finite element function accessors # is there something to be change here?
u = GT.discrete_field(V,node_u)
u_faces = GT.each_face(u,dΩ;tabulate=(GT.value,))
face = 4
u_face = u_faces[face]
u_points = GT.each_point(u_face)
lpoint = 3
u_point = u_points[lpoint]
GT.field(GT.value,u_point)

#ODE right-hand-side 
n_nodes = length(node_u)
n_points = size(W,2)
point_fz = zeros(n_points, n_nodes)
node_wfz = similar(node_x,Float64) # will this actually be the same size? 
workspace = (;W,V,dΩ,node_wfu,point_fz, f)

#ODE solution
T = 400 # Use 400 for nicer results
ode = DE.ODEProblem(node_u,[0,T]) do args...
    rhs!(args...;workspace)
end


