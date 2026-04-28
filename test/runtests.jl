using DrWatson, Test
@quickactivate "neural_field"


include(srcdir("utils.jl"))

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

ti = time() - ti
println("\nTest took total time of:")
println(round(ti/60, digits = 3), " minutes")
