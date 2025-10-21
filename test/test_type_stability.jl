using TSPSolvers
using Test

@testset "Type Stability Tests" begin
    M = [
        0  16   7   14
        16   0   3    5
        7   3   0   16
        14   5  16    0
    ]

    @testset "solve_tsp type stability" begin
        for algo in TSPSolvers.supported_algorithms
            @testset "$algo" begin
                # Test that solve_tsp returns a type-stable result
                result = @inferred solve_tsp(M; algorithm=algo, firstcity=1)
                @test result isa Tuple{Vector{Int}, Int}
            end
        end
    end

    @testset "lowerbound type stability" begin
        D = Matrix{Int}(M)

        # Test with different parameter combinations
        @testset "default parameters" begin
            lb = @inferred lowerbound(D)
            @test lb isa Int
        end

        @testset "with max_iter" begin
            lb = @inferred lowerbound(D; max_iter=100)
            @test lb isa Int
        end

        @testset "with ub" begin
            lb = @inferred lowerbound(D; ub=35)
            @test lb isa Int
        end

        @testset "with both parameters" begin
            lb = @inferred lowerbound(D; max_iter=100, ub=35)
            @test lb isa Int
        end
    end

    @testset "heldkarp_bound type stability" begin
        D = Matrix{Int}(M)
        n = size(D, 1)
        S = collect(2:n)
        depot = 1
        max_iter = 100
        ub = 35

        # Test the internal heldkarp_bound function (returns Float64 by design)
        lb = @inferred TSPSolvers.heldkarp_bound(D, S; depot=depot, max_iter=max_iter, ub=ub)
        @test lb isa Float64
    end
end

