using DrWatson
import GalerkinToolkit as GT
import PartitionedArrays as PA
import MPI; MPI.Init()
import JSON
using LinearAlgebra
using RecursiveArrayTools: ArrayPartition

include(srcdir("equations.jl"))

function prepare_setup_on_rank(rank, np, nx, ny, num_layers)
    setup = mesh_setup(nx, ny)

    (; gn_x) = setup
    ngn = length(gn_x)

    gn_partition_sequential = PA.uniform_partition(1:np, np, ngn)

    nids = gn_partition_sequential[rank]
    ln_gn = PA.local_to_global(nids)

    return restrict_setup(setup, ln_gn, num_layers)
end

function mesh_setup(nx, ny)
    mesh = GT.cartesian_mesh((0,1,0,1),(nx,ny))
    Ω = GT.interior(mesh)
    dΩ = GT.quadrature(Ω,2)
    V = GT.lagrange_space(Ω,1)
    V_faces= GT.each_face(V,dΩ;tabulate=(GT.value,GT.gradient))
    fq_fn_I = transpose(V_faces.accessor.reference_space_face.workspace.values[1])
    fq_fn_dI = transpose(V_faces.accessor.reference_space_face.workspace.gradients[1])
    fq_refdy = GT.weights(GT.reference_quadratures(dΩ)[1])
    gf_fn_gn = GT.face_nodes(V)
    gn_x = GT.node_coordinates(V)
    (;fq_fn_I,fq_fn_dI,fq_refdy,gf_fn_gn,gn_x)
end   

function restrict_setup(setup, ln_gn, num_layers)
    (;fq_fn_I, fq_fn_dI, fq_refdy, gf_fn_gn, gn_x) = setup

    ln_gn = collect(ln_gn)

    ln_x = gn_x[ln_gn]

    ngf = length(gf_fn_gn)
    ngn = length(gn_x)

    nfq, nfn = size(fq_fn_I)
    nln = length(ln_x)

    # Full global u vector, needed later for dz_1.
    gn_u = zeros(ngn)
    fn_ln_zl = zeros(nfn, nln)
    fq_ln_zl_tilde = zeros(nfq, nln)

    return (;
        fq_fn_I,
        fq_fn_dI,
        fq_refdy,
        gf_fn_gn,
        gn_x,
        ln_gn,
        ln_x,
        ngf,
        ngn,
        gn_u,
        fn_ln_zl,
        fq_ln_zl_tilde,
        num_layers
    )
end

function setup_ln_u0(setup)
    (;ln_x) = setup
    return φ.(ln_x)
end

function setup_ln_z0(setup, ln_u0)
    (; gn_x) = setup
    gn_0 = φ.(gn_x)
    return gn_0 * transpose(ln_u0)
end


function all_reduce!(p_gq_u::AbstractArray)
    gq_u_1 = p_gq_u[1]
    for p in 2:length(p_gq_u)
        gq_u = p_gq_u[p]
        for gq in 1:length(gq_u)
            gq_u_1[gq] += gq_u[gq]
        end
    end
    for p in 2:length(p_gq_u)
        gq_u = p_gq_u[p]
        for gq in 1:length(gq_u)
            gq_u[gq] = gq_u_1[gq]
        end
    end
end

function all_reduce!(p_gq_u::PA.DebugArray)
    all_reduce!(p_gq_u.items)
end

function all_reduce!(p_gn_u::PA.MPIArray)
    buf = p_gn_u.item
    op = +
    comm = p_gn_u.comm
    MPI.Allreduce!(buf,op,comm)
    p_gn_u
end

function rhs_IW!(setup, gn_ln_zl,ln_du, ln_u)
    (;fq_fn_I, fq_fn_dI, fq_refdy, gf_fn_gn, gn_x, ln_gn, ln_x, ngf, gn_u, fq_ln_zl_tilde, fn_ln_zl) = setup

    # I - preparation 
    nfq, nfn = size(fq_fn_I)
    nln = size(gn_ln_zl, 2)

    # W - preparation
    fill!(ln_du,0)
    fill!(gn_u,0)

    Ty = eltype(gn_x) # type of coordinates stored in gn_x 
    zy = zero(Ty) # zero coordinate vectore
    zJt = zero(zy*transpose(zy)) # zero matrix-like object for the Jacobian transpose 
    fn_y = zeros(Ty,nfn) # physical coordinates of the nodes of one face 
    fq_y = zeros(Ty,nfq) # physical coordinates of the quadrature points of one face
    fq_dy = zeros(nfq) # physical quadrature weights for one face. 

    for gf in 1:ngf # loop over all faces
        # z_tilde calculation
        fn_gn = gf_fn_gn[gf]
        for fn in 1:nfn
            gn = fn_gn[fn]
            fn_y[fn] = gn_x[gn]
            for ln in 1:nln
                fn_ln_zl[fn,ln] = gn_ln_zl[gn, ln]
            end     
        end     
        mul!(fq_ln_zl_tilde, fq_fn_I, fn_ln_zl)

        for fq in 1:nfq
            # preparation 
            y = zy
            Jt = zJt
            for fn in 1:nfn
                y += fq_fn_I[fq,fn]*fn_y[fn]
                Jt += fn_y[fn]*fq_fn_dI[fq,fn]
            end
            fq_y[fq] = y
            refdy = fq_refdy[fq]
            fq_dy[fq] = abs(det(Jt))*refdy

            # w*f(z_tilde) calculation
            for ln in 1:nln
                x = ln_x[ln]
                dy = fq_dy[fq] 
                wnq = w(x,y)*dy
                zqn = fq_ln_zl_tilde[fq, ln]  
                ln_du[ln] += wnq * f(zqn)
            end    
        end
    end

    for ln in 1:nln
        ln_du[ln] -= ln_u[ln]
        # writing global u values to the gn_u (for all reduce to work properly)
        gn = ln_gn[ln]
        gn_u[gn] = ln_u[ln]
    end    
end 

function rhs_Zall!(setup, z_all, dz_all, gn_u)
    (;gn_x, ln_x, num_layers) = setup
    ngn, nln = size(z_all[1])

    # layer 1
    for i in 1:nln 
        y = ln_x[i]
        for j in 1:ngn
            x = gn_x[j]
            dz_all[1][j,i] = α(x, y, num_layers)*(gn_u[j] - z_all[1][j,i]) 
        end
    end

    # layer 2:all
    for layer in 2:num_layers
        for i in 1:nln  
            y = ln_x[i]
            for j in 1:ngn
                x = gn_x[j]
                dz_all[layer][j,i] = α(x, y, num_layers) *(z_all[layer-1][j,i] - z_all[layer][j,i]) 
            end     
        end
    end     
end     


function rhs!(duz,uz, p_setup, elap_rhs)
    elap_rhs[1] = @elapsed begin 
        num_layers = PA.getany(map(setup->setup.num_layers, p_setup))
        zl = uz.x[num_layers + 1]
        z_all = uz.x[2:num_layers + 1]
        du = duz.x[1]
        dz_all = duz.x[2:num_layers + 1]
        u = uz.x[1]

        p_ln_du = PA.local_values(du)
        p_ln_u = PA.local_values(u)
    end  
    elap_rhs[2] = @elapsed foreach(rhs_IW!, p_setup, zl, p_ln_du, p_ln_u)
    elap_rhs[3] = @elapsed p_gn_u = map(setup-> setup.gn_u,p_setup)
    elap_rhs[4] = @elapsed all_reduce!(p_gn_u)
    elap_rhs[5] = @elapsed begin 
        p_z_all = map(tuple, z_all...)
        p_dz_all = map(tuple, dz_all...)
        foreach(rhs_Zall!,p_setup, p_z_all, p_dz_all, p_gn_u)
    end
end

function main(backend,np,nx,ny,num_layers; save = true, file_name = title)
    ranks = backend(1:np) # creates a partitioned representation of the ranks /process IDs

    # Setup
    elap_setup = zeros(1)
    
    elap_setup[1] = @elapsed p_setup = map(ranks) do rank
        prepare_setup_on_rank(rank, np, nx, ny, num_layers)
    end

    mem = Base.summarysize(p_setup)

    # Initial condition
    ngn = PA.getany(map(setup->setup.ngn,p_setup))
    gn_partition = PA.uniform_partition(ranks, np, ngn)
    p_ln_u0  = map(setup_ln_u0, p_setup)
    
    u0 = PA.PVector(p_ln_u0, gn_partition)
    u = similar(u0)
    du = similar(u0)
    
    z_layers = [map(setup_ln_z0, p_setup, p_ln_u0) for _ in 1:num_layers]
    dz_layers = [map(z -> zero(z), z_layers[layer]) for layer in 1:num_layers]
    
    uz = ArrayPartition(u, z_layers...)
    duz = ArrayPartition(du, dz_layers...)

    nr = 1
    r_elap_rhs = [zeros(5) for _ in 1:nr]
    for r in 1:nr
        # reset the repetition
        elap_rhs = r_elap_rhs[r]
        copy!(u,u0)
        rhs!(duz,uz,p_setup,elap_rhs)
    end

    if save == true
        elap = Dict{Symbol,Vector{Vector{Float64}}}()
        elap[:setup] = [elap_setup]
        elap[:rhs] = r_elap_rhs
        elap[:mem] = [[mem]]
        p_elap_main = PA.gather(map(_->elap,ranks))
        file_title = file_name(nx, ny, np, num_layers)
        PA.map_main(p_elap_main) do p_elap
            open("data/exp_raw/mpi_exp/$file_title.json", "w") do io
                JSON.print(io, p_elap)
            end
        end

        rank = MPI.Comm_rank(MPI.COMM_WORLD)
        comm = MPI.COMM_WORLD

        local_peak_gb = Sys.maxrss() / 1000000000
        maximum_peak_gb = MPI.Allreduce(local_peak_gb, MPI.MAX, comm)

        if rank == 0
            println("Experiment finished. Results saved as $file_title")
            println("Maximum peak memory across all ranks was $maximum_peak_gb GB")
        end

    end     
end

function title(nx,ny,np,num_layers)
   "results_nx$(nx)ny$(ny)np$(np)num_layers$num_layers"
end

function main_mpi(nx,ny, num_layers)
    PA.with_mpi() do backend
        comm = MPI.COMM_WORLD
        np = MPI.Comm_size(comm)
        main(backend,np,3,3,num_layers;save = false, file_name = title)
        MPI.Barrier(comm)
        main(backend,np,nx,ny,num_layers;save = true, file_name = title)
        MPI.Barrier(comm) 
    end
end

function main_debug(nx,ny,np,num_layers)
    PA.with_debug() do backend
        main(backend,np,nx,ny,num_layers;save = true, file_name = title)
    end
end