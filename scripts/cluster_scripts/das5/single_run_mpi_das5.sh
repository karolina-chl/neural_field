#!/bin/bash 

NP=$SLURM_NTASKS
NX=$1
NY=$2
NUM_LAYERS=$3

source das5/modules.sh

echo "Running the job with NP=$NP, NX=$NX, NY=$NY, NUM_LAYERS=$NUM_LAYERS"
echo "Slurm tasks per node $SLURM_NTASKS_PER_NODE"

srun -n $SLURM_NTASKS hostname | sort | uniq -c

mpirun -np $NP \
    julia --project=. scripts/mpi_experiment.jl "$NX" "$NY" "$NUM_LAYERS"