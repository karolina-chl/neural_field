#!/bin/bash

NP=$SLURM_NTASKS
NX=$1
NY=$2
NUM_LAYERS=$3

source snellius/modules.sh

echo "Running the job with NP=$NP, NX=$NX, NY=$NY, NUM_LAYERS=$NUM_LAYERS"

srun -n $SLURM_NTASKS hostname | sort | uniq -c

SECONDS=0

srun \
    --ntasks=$NP \
    --cpu-bind=cores \
    --distribution=block:cyclic \
        julia \
            --threads=1 \
            -O3 \
            --check-bounds=no \
            --project=. \
            scripts/mpi_experiment.jl $NX $NY $NUM_LAYERS

echo "Total runtime: $SECONDS"