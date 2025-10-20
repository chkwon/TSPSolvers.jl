using TSPSolvers, Test, LinearAlgebra

function test()
    M = [
        0  16   7   14
        16   0   3    5
        7   3   0   16
        14   5  16    0 
    ]

    # @assert issymmetric(M)

    for algo in TSPSolvers.supported_algorithms
        @code_warntype solve_tsp(M; algorithm=algo, firstcity=rand(1:4))
    end
end

@time test()

