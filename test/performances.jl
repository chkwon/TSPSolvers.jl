# using TSPSolvers

using Printf, Statistics 

include("../src/main.jl")

function perf_test(n::Int; rep=100, algorithms = supported_algorithms)

    algo = algorithms
    times = zeros(length(algo))
    costs = Matrix{Int}(undef, length(algo), rep)

    for i in 1:rep
        M = rand_int_dist_mtx(n)

        for k in eachindex(algo)
            times[k] += @elapsed tour, cost = solve_tsp(M; algorithm=algo[k])    
            costs[k, i] = cost
        end

    end

    best_costs = minimum(costs, dims=1)

    gaps = Matrix{Float64}(undef, length(algo), rep)
    
    for i in 1:rep
        for k in eachindex(algo)
            gaps[k, i] = (costs[k, i] - best_costs[i]) / best_costs[i] * 100
        end
    end



    println("(n = $n, repetition = $rep)")
    @printf("%20s \t %12s \t %10s \t %10s \n", "Algorithm", "Avg Time (s)", "Avg Cost", "Avg Gap (%)")  
    @printf("%20s \t %12s \t %10s \t %10s \n", "-"^20, "-"^12, "-"^10, "-"^10)  
    for k in eachindex(algo)
        @printf("%20s \t %12.5f \t %10d \t %10.3f\n", algo[k], times[k]/rep, mean(costs[k, :]), mean(gaps[k, :]))
    end

end

perf_test(10; rep=5)
perf_test(10; rep=5, algorithms=["TwoOpt", "LKH"])

perf_test(200; rep=10)
