#!/bin/bash

mkdir -p data/exp_raw/mpi_exp

layer_values_arr=(2)
proc_arr=(1 2 4 8 16 32 64 128 256 512 1024)

for layer in ${layer_values_arr[@]}
do
    for proc in ${proc_arr[@]}
    do
        NP=$proc
        NX=100
        NY=100
        NUM_LAYERS=$layer

        sbatch \
            --partition=fat_rome \
            --ntasks=$NP \
            --time=02:00:00 \
            --exclusive \
            scripts/cluster_scripts/snellius/single_run_mpi_snellius.sh $NX $NY $NUM_LAYERS
    done    
done  