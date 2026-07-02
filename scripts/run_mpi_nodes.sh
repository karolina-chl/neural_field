#!/bin/bash

export PATH="/var/scratch/tln716/julia/julia-1.12.4/bin:$PATH"
export JULIA_DEPOT_PATH="/var/scratch/tln716/julia_depot"

mkdir -p data/exp_raw/mpi_exp

nx_arr=(60)
ny_arr=(60)
proc_arr=(1)

for idx in ${!nx_arr[@]}
do
    for proc in ${proc_arr[@]}
    do
        NP=$proc
        NX=${nx_arr[$idx]}
        NY=${ny_arr[$idx]}
        NUM_LAYERS=512

        sbatch \
            --ntasks=$NP \
            --ntasks-per-node=4 \
            scripts/single_run_mpi.sh $NX $NY $NUM_LAYERS
    done    
done  