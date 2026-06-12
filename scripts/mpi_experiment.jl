using MPI
using DrWatson

@quickactivate "neural_field"

include(srcdir("parallel_delayed.jl"))

nx = parse(Int, ARGS[1])
ny = parse(Int, ARGS[2])
num_layers = parse(Int, ARGS[3])

main_mpi(nx, ny, num_layers)