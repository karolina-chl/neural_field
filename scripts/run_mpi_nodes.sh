#!/bin/bash

export PATH="/var/scratch/tln716/julia/julia-1.12.4/bin:$PATH"
export JULIA_DEPOT_PATH="/var/scratch/tln716/julia_depot"

power_arr=(0 1 2 3 4 5)
proc_arr=(2 4 8 16 32 64 128)

for i in ${power_arr[@]}
do
    for proc in ${proc_arr[@]}
    do
        NP=$proc
        NX=$((30*2**i))
        NY=$((30*2**i))
        NUM_LAYERS=2

        sbatch --ntasks=$NP run_one_mpi.sh $NX $NY $NUM_LAYERS
    done    
done  