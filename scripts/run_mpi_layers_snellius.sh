#!/bin/bash

mkdir -p data/exp_raw/mpi_exp

layer_values_arr=(32)
proc_arr=(128)

for layer in ${layer_values_arr[@]}
do
    for proc in ${proc_arr[@]}
    do
        NP=$proc
        NX=60
        NY=60
        NUM_LAYERS=$layer

        sbatch \
            --partition=fat_rome \
            --ntasks=$NP \
            --time=00:15:00 \
            --exclusive \
            scripts/single_run_mpi_snellius.sh $NX $NY $NUM_LAYERS
    done    
done  