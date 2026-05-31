using DrWatson
import GalerkinToolkit as GT
import PartitionedArrays as PA
using MPI
import JSON
using LinearAlgebra
using RecursiveArrayTools: ArrayPartition

MPI.Init()

include(srcdir("equations.jl"))
include(srcdir("MPI/functions.jl"))

#Input 
nx = 2
ny = 2
num_layers = 2

comm = MPI.COMM_WORLD
nranks = MPI.Comm_size(comm)
rank = MPI.Comm_rank(comm)
root = 0

# setup in main and bradcast  
setup_root = Ref(0) # consider setting this as an actual size 
ln_gn_all = Ref(0)
if rank == root
    setup_root = do_setup(nx, ny, num_layers)
    (;gn_x) = setup_root
    ngn = length(gn_x) 
    ln_gn_all = partition_nodes(nranks,ngn) 
end
setup_root = MPI.bcast(setup_root,comm;root) # should I use Bcast! ?
ln_gn_all = MPI.bcast(ln_gn_all,comm;root) # maybe they don't all need full mapping 

ln_gn = ln_gn_all[rank+1]
print("My rank is $rank and my ln_gn is $ln_gn")

MPI.Finalize()
