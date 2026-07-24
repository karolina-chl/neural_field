#!/bin/bash

mkdir -p data/exp_raw/mpi_exp

layer_values_arr=(2 8 31 70 158)
proc_arr=(1 2 4 8 16 32 64)

for layer in ${layer_values_arr[@]}
do
    for proc in ${proc_arr[@]}
    do
        NP=$proc
        NX=100
        NY=100
        NUM_LAYERS=$layer

        sbatch \
            --ntasks=$NP \
            --ntasks-per-node=2 \
            scripts/cluster_scripts/das5/single_run_mpi_das5.sh $NX $NY $NUM_LAYERS
    done    
done  