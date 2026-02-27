using DrWatson
# import GalerkinToolkit as GT
# import GLMakie as Makie
# import DifferentialEquations as DE
# import ProgressMeter as PM
# using LinearAlgebra
# using SparseArrays
# using RecursiveArrayTools

@quickactivate "neural_field"

include(srcdir("FEM.jl"))
include(srcdir("equations.jl"))
include(srcdir("main_function.jl"))

main()
