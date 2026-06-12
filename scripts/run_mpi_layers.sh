#!/bin/bash


export PATH="/var/scratch/tln716/julia/julia-1.12.4/bin:$PATH"
export JULIA_DEPOT_PATH="/var/scratch/tln716/julia_depot"

layer_values_arr=(2 3 4 5 6 7 8 9 10)
proc_arr=(2 4 8 16 32 64 128)

for layer in ${layer_values_arr[@]}
do
    for proc in ${proc_arr[@]}
    do
        NP=$proc
        NX=50
        NY=50
        NUM_LAYERS=$layer

        sbatch --ntasks=$NP run_one_mpi.sh $NX $NY $NUM_LAYERS
    done    
done  