using DrWatson, Test
@quickactivate "neural_field"


include(srcdir("utils.jl"))
include(srcdir("parallel_delayed.jl"))


# Run test suite
println("Starting tests")
ti = time()

@testset "model correctnes" begin
    duz_1, uz_1, p, t, workspace_1 = setup_rhs_delayed_problem(;
    L = 10,
    num_el = 100,
    num_layer = 2,
    new_fun = true
    )

    duz_2, uz_2, p, t, workspace_2 = setup_rhs_delayed_problem(;
    L = 10,
    num_el = 100,
    num_layer = 2,
    new_fun = false
    )
    run_rhs_delayed(duz_1, uz_1, p, t, workspace_1; n=1, new_fun = true)
    run_rhs_delayed(duz_2, uz_2, p, t, workspace_2; n=1, new_fun = false)

    for i in eachindex(duz_1.x)
        @test duz_1.x[i] == duz_2.x[i]
    end
end

@testset "correctness of parallel implementation" begin
    
    duz_s, uz_s, p_s, t_s, workspace_s = setup_for_parallel_test(3,3,2) #nx, ny, num_layer
    test_parallel_implementation(duz_s, uz_s, p_s, t_s, workspace_s;nx=3,ny=3,np=3,num_layers=2, parallel = false)
    duz_p_debug = test_parallel_implementation(duz_s, uz_s, p_s, t_s, workspace_s;nx=3,ny=3,np=3,num_layers=2, parallel = true)

    duz_p = materialize_debug_duz(duz_p_debug)

    for i in eachindex(duz_s.x)
        println("==")
        @test duz_s.x[i] == duz_p.x[i]
        maxdiff = maximum(abs.(duz_s[i] .- duz_p[i]))
        println("  max abs diff:   ", maxdiff)
        @test isapprox(duz_s.x[i], duz_p.x[i]; atol = 1e-10, rtol = 1e-10)
    end

end

ti = time() - ti
println("\nTest took total time of:")
println(round(ti/60, digits = 3), " minutes")

