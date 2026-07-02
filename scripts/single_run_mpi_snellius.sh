#!/bin/bash

NP=$SLURM_NTASKS
NX=$1
NY=$2
NUM_LAYERS=$3

module load 2023
module load juliaup/1.14.5-GCCcore-12.3.0
export PATH="$HOME/.julia/bin:$PATH"

echo "Running the job with NP=$NP, NX=$NX, NY=$NY, NUM_LAYERS=$NUM_LAYERS"

srun -n $SLURM_NTASKS hostname | sort | uniq -c

mpiexecjl -n $NP julia --project=. scripts/mpi_experiment.jl $NX $NY $NUM_LAYERS