#!/bin/bash

mkdir -p data/exp_raw/mpi_exp

nx_arr=(160)
ny_arr=(160)
proc_arr=(2048 4096 8192)

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
            --time=00:15:00 \
            scripts/single_run_mpi.sh $NX $NY $NUM_LAYERS
    done    
done  