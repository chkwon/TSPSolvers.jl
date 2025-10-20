using TSPSolvers
using Test
using TSPLIB

@testset "Held-Karp lower bound" begin
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

    tsp_path = joinpath(@__DIR__, "gr17.tsp")
    if isfile(tsp_path)
        tsp = TSPLIB.readTSP(tsp_path)
        push!(instances, ("gr17", Int.(tsp.weights)))
    end

    for (name, dist) in instances
        @testset "$name" begin
            lb = lowerbound(dist)
            _, opt_cost = TSPSolvers.solve_tsp(dist; algorithm="Concorde")
            @test lb isa Int
            @test lb <= opt_cost
        end
    end
end
