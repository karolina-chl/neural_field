#!/bin/bash

export PATH="/var/scratch/tln716/julia/julia-1.12.4/bin:$PATH"
export JULIA_DEPOT_PATH="/var/scratch/tln716/julia_depot"

mkdir -p data/exp_raw/mpi_exp

nx_arr=(160)
ny_arr=(80)
proc_arr=(1 2 4 8 16 32 64 128 256)

for idx in ${!nx_arr[@]}
do
    for proc in ${proc_arr[@]}
    do
        NP=$proc
        NX=${nx_arr[$idx]}
        NY=${ny_arr[$idx]}
        NUM_LAYERS=2

        sbatch \
            --ntasks=$NP \
            --ntasks-per-node=8 \
            scripts/single_run_mpi.sh $NX $NY $NUM_LAYERS
    done    
done  