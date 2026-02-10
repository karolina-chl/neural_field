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

### Initial conditions L=1
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





