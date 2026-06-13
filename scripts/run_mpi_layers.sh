#!/bin/bash


export PATH="/var/scratch/tln716/julia/julia-1.12.4/bin:$PATH"
export JULIA_DEPOT_PATH="/var/scratch/tln716/julia_depot"
mkdir -p data/exp_raw/mpi_exp

layer_values_arr=(2 3 4 5 6 7 8)
proc_arr=(1 2 4 8 16 32 64)

for layer in ${layer_values_arr[@]}
do
    for proc in ${proc_arr[@]}
    do
        NP=$proc
        NX=50
        NY=50
        NUM_LAYERS=$layer

        sbatch --ntasks=$NP scripts/single_run_mpi.sh $NX $NY $NUM_LAYERS
    done    
done  