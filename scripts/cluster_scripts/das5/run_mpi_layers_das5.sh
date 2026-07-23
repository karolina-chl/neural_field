#!/bin/bash

export PATH="/var/scratch/tln716/julia/julia-1.12.4/bin:$PATH"
export JULIA_DEPOT_PATH="/var/scratch/tln716/julia_depot"
mkdir -p data/exp_raw/mpi_exp

layer_values_arr=(780)
proc_arr=(8)

for layer in ${layer_values_arr[@]}
do
    for proc in ${proc_arr[@]}
    do
        NP=$proc
        NX=50
        NY=50
        NUM_LAYERS=$layer

        sbatch \
            --ntasks=$NP \
            --ntasks-per-node=4 \
            scripts/single_run_mpi_das5.sh $NX $NY $NUM_LAYERS
    done    
done  