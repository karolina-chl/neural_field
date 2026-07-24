#!/bin/bash

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
            scripts/cluster_scripts/das5/single_run_mpi_das5.sh $NX $NY $NUM_LAYERS
    done    
done  