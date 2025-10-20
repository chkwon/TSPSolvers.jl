using TSPSolvers
using Test
using TSPLIB

@testset "Held-Karp lower bound" begin
    @testset "Basic functionality" begin
        instances = [
            ("simple4", [
                0  16   7  14;
               16   0   3   5;
                7   3   0  16;
               14   5  16   0
            ]),
            ("simple5", [
                0  12  10  19   8;
               12   0   3   7   2;
               10   3   0   6  20;
               19   7   6   0   4;
                8   2  20   4   0
            ]),
            ("simple6", [
                 0  11  17   6  15  12;
                11   0  13   9   7  10;
                17  13   0  21   8  14;
                 6   9  21   0  12   5;
                15   7   8  12   0  16;
                12  10  14   5  16   0
            ])
        ]

        for (name, dist) in instances
            @testset "$name" begin
                lb = lowerbound(dist)
                _, opt_cost = TSPSolvers.solve_tsp(dist; algorithm="Concorde")
                @test lb isa Int
                @test lb <= opt_cost
                @test lb > 0

                # Test that lower bound is reasonably tight (within 50% of optimal)
                gap_ratio = (opt_cost - lb) / opt_cost
                @test gap_ratio < 0.5
            end
        end
    end

    @testset "TSPLIB instances" begin
        tsp_path = joinpath(@__DIR__, "gr17.tsp")
        if isfile(tsp_path)
            @testset "gr17" begin
                tsp = TSPLIB.readTSP(tsp_path)
                dist = Int.(tsp.weights)
                lb = lowerbound(dist)
                _, opt_cost = TSPSolvers.solve_tsp(dist; algorithm="Concorde")
                @test lb isa Int
                @test lb <= opt_cost
                @test lb > 0
            end
        else
            @info "gr17.tsp not found, skipping TSPLIB test"
        end
    end

    @testset "Custom parameters" begin
        M = [
            0  16   7  14;
           16   0   3   5;
            7   3   0  16;
           14   5  16   0
        ]

        # Test with custom max_iter
        lb1 = lowerbound(M; max_iter=50)
        lb2 = lowerbound(M; max_iter=200)
        @test lb1 isa Int
        @test lb2 isa Int
        @test lb1 > 0
        @test lb2 > 0

        # Test with custom upper bound
        _, opt_cost = TSPSolvers.solve_tsp(M; algorithm="Concorde")
        lb3 = lowerbound(M; ub=opt_cost)
        @test lb3 isa Int
        @test lb3 <= opt_cost

        # Test with both parameters
        lb4 = lowerbound(M; max_iter=150, ub=opt_cost)
        @test lb4 isa Int
        @test lb4 <= opt_cost
    end

    @testset "2-opt upper bound integration" begin
        # Test that 2-opt is automatically used when ub is not provided
        M = [
            0  12  10  19   8;
           12   0   3   7   2;
           10   3   0   6  20;
           19   7   6   0   4;
            8   2  20   4   0
        ]

        lb = lowerbound(M)
        _, opt_cost = TSPSolvers.solve_tsp(M; algorithm="Concorde")
        @test lb isa Int
        @test lb <= opt_cost
        @test lb > 0
    end

    @testset "Edge cases" begin
        # Test with 3x3 matrix (minimum meaningful TSP)
        M3 = [
            0  10  15;
           10   0  20;
           15  20   0
        ]
        lb3 = lowerbound(M3)
        @test lb3 isa Int
        @test lb3 <= 45  # Optimal tour: 1-2-3-1 = 10+20+15 = 45
        @test lb3 > 0
    end
end
