#!/bin/bash

NP=$SLURM_NTASKS
NX=$1
NY=$2
NUM_LAYERS=$3

source snellius/modules.sh

echo "Running the job with NP=$NP, NX=$NX, NY=$NY, NUM_LAYERS=$NUM_LAYERS"

srun -n $SLURM_NTASKS hostname | sort | uniq -c

srun --ntasks=$NP \
    --cpu-bind=cores \
    --distribution=block:cyclic \
    julia --project=. scripts/mpi_experiment.jl $NX $NY $NUM_LAYERS