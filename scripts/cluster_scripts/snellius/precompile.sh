#!/bin/bash
#SBATCH --partition=rome
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=00:30:00

source snellius/modules.sh

julia \
    --threads=1 \
    -O3 \
    --check-bounds=no \
    --project=. \
    -e 'using Pkg; Pkg.instantiate(); Pkg.precompile()'