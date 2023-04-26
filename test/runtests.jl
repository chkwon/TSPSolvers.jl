using TSPSolvers
using Test

@testset "TSPSolvers.jl" begin
    # Write your tests here.
    include("test_shift_tour.jl")
    include("test_simple_TSP.jl")
end