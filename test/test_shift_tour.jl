using TSPSolvers
using Test 
using LinearAlgebra
using Random 

@testset "Test Shift Tour" begin
    for _ in 1:10
        n = rand(1:100)
        tour = shuffle(1:n)
        first_city = rand(1:n)
        TSPSolvers.shift_tour!(tour, first_city)
        @test tour[1] == first_city
    end
end