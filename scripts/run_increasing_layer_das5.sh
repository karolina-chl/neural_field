#!/bin/bash

export PATH="/var/scratch/tln716/julia/julia-1.12.4/bin:$PATH"
export JULIA_DEPOT_PATH="/var/scratch/tln716/julia_depot"
mkdir -p data/exp_raw

layer_values_arr=(2)

for layer in ${layer_values_arr[@]}
do
    julia --project=. scripts/exp_increasing_layer_number.jl $layer 
done