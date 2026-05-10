module Demo

import GalerkinToolkit as GT
import PartitionedArrays as PA
import MPI; MPI.Init()
import JSON
using LinearAlgebra

function prepare_setups_on_main(np,nx,ny)
    setup = do_setup(nx,ny)
    (;gn_x) = setup
    ngn = length(gn_x)
    gn_partition_sequential = PA.uniform_partition(1:np,np,ngn)
    p_fem_setup = map(gn_partition_sequential) do nids
        ln_gn = PA.local_to_global(nids)
        gn_ln = PA.global_to_local(nids)
        restrict_setup(setup,ln_gn,gn_ln)
    end
end

function do_setup(nx,ny)
    mesh = GT.cartesian_mesh((0,1,0,1),(nx,ny))
    Ω = GT.interior(mesh)
    dΩ = GT.quadrature(Ω,2)
    V = GT.lagrange_space(Ω,2)
    V_faces= GT.each_face(V,dΩ;tabulate=(GT.value,GT.gradient))
    fq_fn_I = transpose(V_faces.accessor.reference_space_face.workspace.values[1])
    fq_fn_dI = transpose(V_faces.accessor.reference_space_face.workspace.gradients[1])
    fq_refdy = GT.weights(GT.reference_quadratures(dΩ)[1])
    gf_fn_gn = GT.face_nodes(V)
    gn_x = GT.node_coordinates(V)
    nfq = size(fq_fn_I,1)
    ngf = length(gf_fn_gn)
    ngq = nfq*ngf
    ngn = length(gn_x)
    (;fq_fn_I,fq_fn_dI,fq_refdy,gf_fn_gn,gn_x)
end

function restrict_setup(setup,ln_gn,gn_ln)
    (;fq_fn_I,fq_fn_dI,fq_refdy,gf_fn_gn,gn_x) = setup
    ln_x = gn_x[ln_gn] # Do not use view
    ngf = length(gf_fn_gn)
    ngn =length(gn_x)
    gn_in_part = fill(false,ngn) #!!!
    gn_in_part[ln_gn] .= true
    gf_in_part = fill(false,ngf) # !!!
    for gf in 1:ngf #!!!
        fn_gn= gf_fn_gn[gf]
        in_part = false
        for gn in fn_gn
            in_part = in_part || gn_in_part[gn]
        end
        gf_in_part[gf]= in_part
    end
    lf_gf = findall(gf_in_part) # !!!
    lf_fn_gn = PA.jagged_array(gf_fn_gn[lf_gf])
    data = lf_fn_gn.data
    for k in 1:length(data)
        gn = data[k]
        ln = gn_ln[gn]
        data[k] = ln
    end
    ptrs = lf_fn_gn.ptrs
    lf_fn_ln =PA.jagged_array(data,ptrs)
    nfq = size(fq_fn_I,1)
    ngq = ngf * nfq
    gq_u = zeros(ngq)# !!!
    (;fq_fn_I,fq_fn_dI,fq_refdy,lf_fn_ln,ln_x,ngf,ngn,gq_u,lf_gf,gf_fn_gn,gn_x)
end

function setup_ln_u0(setup)
    (;ln_x) = setup
    ln_u0 = sum.(ln_x) # use the right function
    ln_u0
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

function all_reduce!(p_gq_u::PA.MPIArray)
    buf = p_gq_u.item
    op = +
    comm = p_gq_u.comm
    MPI.Allreduce!(buf,op,comm)
    p_gq_u
end

function rhs_I!(setup, ln_u)
    (;fq_fn_I,lf_fn_ln,gq_u,lf_gf)= setup
    nfq,nfn = size(fq_fn_I)
    fq_u = zeros(nfq)
    fn_u = zeros(nfn)
    nlf = length(lf_fn_ln)
    for lf in 1:nlf
        fn_ln = lf_fn_ln[lf]
        for fn in 1:nfn
            ln = fn_ln[fn]
            if ln != 0
                fn_u[fn]  = ln_u[ln]
            else
                fn_u[fn] = 0
            end
        end
        mul!(fq_u,fq_fn_I,fn_u)
        gf = lf_gf[lf]
        for fq in 1:nfq
            gq = (gf-1)*nfq + fq
            gq_u[gq] = fq_u[fq]
        end
    end
end

function rhs_W!(setup, ln_du)
    (;fq_fn_I,fq_fn_dI,fq_refdy,lf_fn_ln,ln_x,ngf,ngn,gq_u,lf_gf,gf_fn_gn,gn_x) = setup

    fill!(ln_du,0)

    nln = length(ln_x)
    nfq,nfn = size(fq_fn_I)
    ngq = ngf * nfq
    Ty = eltype(gn_x)
    zy = zero(Ty)
    zJt = zero(zy*transpose(zy))
    fn_y = zeros(Ty,nfn)
    fq_y = zeros(Ty,nfq)
    fq_dy = zeros(nfq)
    for gf in 1:ngf
        fn_gn = gf_fn_gn[gf]
        for fn in 1:nfn
            gn = fn_gn[fn]
            fn_y[fn] = gn_x[gn]
        end
        for fq in 1:nfq
            y = zy
            for fn in 1:nfn
                y += fq_fn_I[fq,fn]*fn_y[fn]
            end
            fq_y[fq] = y
        end
        for fq in 1:nfq
            Jt = zJt
            for fn in 1:nfn
                Jt += fn_y[fn]*fq_fn_dI[fq,fn]
            end
            refdy = fq_refdy[fq]
            fq_dy[fq] = abs(det(Jt))*refdy
        end
        for fq in 1:nfq
            y = fq_y[fq]
            dy = fq_dy[fq]
            gq = (gf-1)*nfq + fq
            u = gq_u[gq]
            for ln in 1:nln
                x = ln_x[ln]
                W = norm(y-x)*dy # Use the right function here
                ln_du[ln] += W*u # Use firing rate here
            end
        end
    end

end

# It is a good idea to measure the time for all lines
function rhs!(du,u,p_setup,elap)
    elap[1] = @elapsed foreach(setup->fill!(setup.gq_u,0),p_setup) # This is negligible as long as N/P >> 1
    elap[2] = @elapsed p_ln_u = PA.local_values(u)
    elap[3] = @elapsed foreach(rhs_I!,p_setup,p_ln_u)
    elap[4] = @elapsed p_gq_u = map(setup->setup.gq_u,p_setup)
    elap[5] = @elapsed all_reduce!(p_gq_u)
    elap[6] = @elapsed p_ln_du = PA.local_values(du)
    elap[7] = @elapsed foreach(rhs_W!,p_setup,p_ln_du)
end

function main(backend,np,nx,ny,title)
    ranks = backend(1:np)

    # Setup
    elap_setup = zeros(2)
    elap_setup[1] = @elapsed p_setup_on_main = PA.map_main(ranks) do _
                                 prepare_setups_on_main(np,nx,ny)
                             end
    elap_setup[2] = @elapsed p_setup = PA.scatter(p_setup_on_main)

    mem = Base.summarysize(p_setup)

    # Initial condition
    p_ln_u0 = map(setup_ln_u0,p_setup)
    ngn = PA.getany(map(setup->setup.ngn,p_setup))
    gn_partition = PA.uniform_partition(ranks,np,ngn)
    u0 = PA.PVector(p_ln_u0,gn_partition)

    u = similar(u0)
    du = similar(u0)
    nr = 10
    r_elap_rhs = [zeros(7) for _ in 1:nr]
    for r in 1:nr
        elap_rhs = r_elap_rhs[r]
        copy!(u,u0)
        rhs!(du,u,p_setup,elap_rhs)
    end

    elap = Dict{Symbol,Vector{Vector{Float64}}}()
    elap[:setup] = [elap_setup]
    elap[:rhs] = r_elap_rhs
    elap[:mem] = [[mem]]
    p_elap_main = PA.gather(map(_->elap,ranks))
    PA.map_main(p_elap_main) do p_elap
        JSON.json("$title.json",p_elap)
    end
    title
end

function title(nx,ny,np)
   "results_nx$(nx)ny$(ny)np$np"
end

function main_mpi(nx,ny)
    PA.with_mpi() do backend
        comm = MPI.COMM_WORLD
        np = MPI.Comm_size(comm)
        main(backend,np,3,3,"warmup")
        MPI.Barrier(comm)
        main(backend,np,nx,ny,title(nx,ny,np))
    end
end

function main_debug(nx,ny,np)
    PA.with_debug() do backend
        main(backend,np,nx,ny,"debug")
    end
end

end # module

Demo.main_debug(3, 3, 2)
