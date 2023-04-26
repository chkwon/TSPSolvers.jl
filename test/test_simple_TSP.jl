using TSPSolvers
using Test 
using LinearAlgebra

@testset "simple TSP" begin
    M = [
         0  16   7   14
        16   0   3    5
         7   3   0   16
        14   5  16    0 
    ]

    @assert issymmetric(M)

    for algo in TSPSolvers.supported_algorithms
        @testset "$algo" begin
            tour, cost = solve_tsp(M; algorithm=algo, firstcity=rand(1:4))
            @test cost == 29
        end
    end

    @testset "TwoOpt Tests" begin
        algo = "TwoOpt"

        firstcity = 3
        tour, cost = solve_tsp(M; algorithm=algo, firstcity=firstcity)
        @test cost == 29        
        @test tour[1] == firstcity

        init_tour = Int[]
        tour, cost = solve_tsp(M; algorithm=algo, init_tour=init_tour)
        @test cost == 29        

        init_tour = [2, 3, 4, 1, 2]
        firstcity = 3
        @test_logs (:warn, "The first city of the initial tour is not $firstcity. It will be ignored.") tour, cost = solve_tsp(M; algorithm=algo, init_tour=init_tour, firstcity=firstcity)
        @test cost == 29      
        
        init_tour = [2, 3, 4, 1]
        tour, cost = solve_tsp(M; algorithm=algo, init_tour=init_tour, firstcity = 3)
        @test cost == 29              
    end

end