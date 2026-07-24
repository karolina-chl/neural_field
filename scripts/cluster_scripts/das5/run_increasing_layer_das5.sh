#!/bin/bash

source das5/modules.sh
mkdir -p data/exp_raw

layer_values_arr=(20)

for layer in ${layer_values_arr[@]}
do
    julia --project=. scripts/exp_increasing_layer_number.jl $layer 
done