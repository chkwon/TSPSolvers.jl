import Concorde
import LKH
import TravelingSalesmanHeuristics
import TSPLIB
using Random

function dist_mat(x::Vector{T}, y::Vector{T}) where T <: Real
    @assert length(x) == length(y)
    n = length(x)
    M = Matrix{T}(undef, n, n)
    @inbounds for i in 1:n, j in 1:n
        M[i, j] = sqrt((x[i]-x[j])^2 + (y[i]-y[j])^2)
    end

    return M
end

function rand_int_dist_mtx(n::Int; digits=3)
    x = rand(n) .* 10 ^ digits
    y = rand(n) .* 10 ^ digits
    dist = round.(Int, dist_mat(x, y))
    return dist
end

function shift_tour!(tour::Vector{Int}, firstcity)
    if tour[1] != firstcity 
        idx = findfirst(x -> x==firstcity, tour)
        circshift!(tour, idx - 1)
    end
end

function rand_tour_for_heuristics(n::Int; firstcity=1)
    tour = shuffle(1:n)
    shift_tour!(tour, firstcity)
    push!(tour, firstcity) # TravelingSalesmanHeuristics requires the first node at the end of the tour again.
    return tour
end



supported_algorithms = [
    "Concorde", 
    "LKH", 
    "NearestNeighbor", 
    "FarthestInsertion", 
    "CheapestInsertion",
    "TwoOpt", 
    "SimulatedAnnealing"
]


function solve_tsp(dist_mtx::Matrix{Int}; algorithm="LKH", firstcity=1, kwargs...) 
    n = size(dist_mtx, 1)

    if algorithm == "Concorde"
        tour, cost = Concorde.solve_tsp(dist_mtx; kwargs...)
        shift_tour!(tour, firstcity)
        return tour, cost

    elseif algorithm == "LKH"
        tour, cost = LKH.solve_tsp(dist_mtx; kwargs...)
        shift_tour!(tour, firstcity)
        return tour, cost

    elseif algorithm == "NearestNeighbor"
        _tour, cost = TravelingSalesmanHeuristics.nearest_neighbor(dist_mtx; firstcity=firstcity, kwargs...)
        tour = _tour[1:end-1]
        return tour, cost

    elseif algorithm == "FarthestInsertion"
        _tour, cost = TravelingSalesmanHeuristics.farthest_insertion(dist_mtx; firstcity=firstcity, kwargs...)
        tour = _tour[1:end-1]
        return tour, cost

    elseif algorithm == "CheapestInsertion"
        _tour, cost = TravelingSalesmanHeuristics.cheapest_insertion(dist_mtx; firstcity=firstcity, kwargs...)
        tour = _tour[1:end-1]      
        return tour, cost  

    elseif algorithm == "TwoOpt"
        init_tour = rand_tour_for_heuristics(n; firstcity=firstcity)
        _tour, cost = TravelingSalesmanHeuristics.two_opt(dist_mtx, init_tour; kwargs...)
        tour = _tour[1:end-1]
        return tour, cost

    elseif algorithm == "SimulatedAnnealing"
        _tour, cost = TravelingSalesmanHeuristics.farthest_insertion(dist_mtx; firstcity=firstcity, kwargs...)
        tour = _tour[1:end-1]
        return tour, cost

    else
        error("Algorithm \"$(algorithm)\" is not supported. Choose from $supported_algorithms.")
    end

end



function solve_tsp(file::String; algorithm="LKH", firstcity=1, kwargs...) 

    if algorithm == "Concorde"
        tour, cost = Concorde.solve_tsp(file; kwargs...)
        shift_tour!(tour, firstcity)
        return tour, cost

    elseif algorithm == "LKH"
        tour, cost =  LKH.solve_tsp(file; kwargs...)
        shift_tour!(tour, firstcity)
        return tour, cost

    else
        tsp = TSPLIB.readTSP(file)
        M = Int.(tsp.weights)    
        return solve_tsp(M; algorithm=algorithm, firstcity=firstcity, kwargs...)
    end

end