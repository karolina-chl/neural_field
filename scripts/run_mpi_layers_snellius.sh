#!/bin/bash

mkdir -p data/exp_raw/mpi_exp

layer_values_arr=(2 32 512)
proc_arr=(1 2 4 8 16 32 64)

for layer in ${layer_values_arr[@]}
do
    for proc in ${proc_arr[@]}
    do
        NP=$proc
        NX=19
        NY=19
        NUM_LAYERS=$layer

        sbatch \
            --partition=fat_rome \
            --ntasks=$NP \
            --time=00:15:00 \
            scripts/single_run_mpi_snellius.sh $NX $NY $NUM_LAYERS
    done    
done  