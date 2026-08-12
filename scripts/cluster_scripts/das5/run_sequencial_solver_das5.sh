#!/bin/bash

source das5/modules.sh
mkdir -p data/exp_raw

layer_values_arr=(2)

for layer in ${layer_values_arr[@]}
do
    julia --project=. scripts/run_sequencial_solver.jl $layer 
done