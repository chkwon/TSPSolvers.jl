using TSPSolvers
using Test 

@testset "simple TSP" begin
    M = [
        0  16   7  14
        16   0   3   5
        7   3   0  16
        14   5  16   0 
    ]

    for algo in TSPSolvers.supported_algorithms
        tour, cost = solve_tsp(M; algorithm=algo)
        @test cost == 29
    end

end