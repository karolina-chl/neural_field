#!/bin/bash
#SBATCH --partition=rome
#SBATCH --time=00:30:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --output=snellius-configure.o
#SBATCH --error=snellius-configure.e
set -e
source snellius/modules.sh
julia --project=. -e 'using Pkg; Pkg.instantiate()'
julia --project=. -e 'using MPIPreferences;use_system_binary(;force=true)'
julia --project=.  -e 'using Pkg; Pkg.precompile()'
julia --project=. -O3 --check-bounds=no -e 'using Pkg; Pkg.precompile()'