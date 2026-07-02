#!/bin/bash

export PATH="/var/scratch/tln716/julia/julia-1.12.4/bin:$PATH"
export JULIA_DEPOT_PATH="/var/scratch/tln716/julia_depot"

NP=$SLURM_NTASKS
NX=$1
NY=$2
NUM_LAYERS=$3

echo "Running the job with NP=$NP, NX=$NX, NY=$NY, NUM_LAYERS=$NUM_LAYERS"
echo "Slurm tasks per node $SLURM_NTASKS_PER_NODE"

srun -n $SLURM_NTASKS hostname | sort | uniq -c

env -u SLURM_NTASKS_PER_NODE \
    mpiexecjl -n $NP julia --project=. scripts/mpi_experiment.jl $NX $NY $NUM_LAYERS