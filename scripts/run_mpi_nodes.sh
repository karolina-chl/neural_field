#!/bin/bash

export PATH="/var/scratch/tln716/julia/julia-1.12.4/bin:$PATH"
export JULIA_DEPOT_PATH="/var/scratch/tln716/julia_depot"

mkdir -p data/exp_raw/mpi_exp

nx_arr=(19 39 79)
proc_arr=(1 2 4 8 16 32 64)

for i in ${nx_arr[@]}
do
    for proc in ${proc_arr[@]}
    do
        NP=$proc
        NX=$i
        NY=$i
        NUM_LAYERS=2

        sbatch --ntasks=$NP scripts/single_run_mpi.sh $NX $NY $NUM_LAYERS
    done    
done  