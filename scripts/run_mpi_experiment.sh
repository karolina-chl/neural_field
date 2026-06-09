#!/bin/bash

export PATH="/var/scratch/tln716/julia/julia-1.12.4/bin:$PATH"
export JULIA_DEPOT_PATH="/var/scratch/tln716/julia_depot"

NP=64
NX=30
NY=30
NUM_LAYERS=2

mpiexecjl -np "$NP" julia --project=. scripts/MPI_experiment.jl "$NX" "$NY" "$NUM_LAYERS"