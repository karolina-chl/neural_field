#!/bin/bash
#SBATCH --time=00:30:00
#SBATCH --nodes=1
#SBATCH --output=das5-configure.o
#SBATCH --error=das5-configure.e
set -e
source das5/modules.sh
julia --project=. -e 'using Pkg; Pkg.instantiate()'
julia --project=. -e 'using MPIPreferences;use_system_binary(;force=true)'
julia --project=.  -e 'using Pkg; Pkg.precompile()'
julia --project=. -O3 --check-bounds=no -e 'using Pkg; Pkg.precompile()'