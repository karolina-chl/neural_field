#!/bin/bash

mkdir -p data/exp_raw/mpi_exp

nx_arr=(200)
ny_arr=(200)
proc_arr=(1 2 4)

for idx in ${!nx_arr[@]}
do
    for proc in ${proc_arr[@]}
    do
        NP=$proc
        NX=${nx_arr[$idx]}
        NY=${ny_arr[$idx]}
        NUM_LAYERS=2

        sbatch \
            --partition=fat_rome \
            --ntasks=$NP \
            --time=03:00:00 \
            --exclusive \
            scripts/cluster_scripts/snellius/single_run_mpi_snellius.sh $NX $NY $NUM_LAYERS
    done    
done  