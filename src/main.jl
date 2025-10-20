import Concorde
import LKH
import Hygese
import TravelingSalesmanHeuristics
import TSPLIB
using LinearAlgebra
using Random


include("utils.jl")
include("performances.jl")
include("heldkarp.jl")

const supported_algorithms = [
    "Concorde", 
    "LKH", 
    "HGS", 
    "NearestNeighbor", 
    "FarthestInsertion", 
    "CheapestInsertion",
    "TwoOpt", 
    "SimulatedAnnealing"
]


function solve_tsp(dist_mtx::Matrix{Int}; algorithm="LKH", firstcity=1, init_tour::Vector{Int}=Int[], kwargs...) 
    n = size(dist_mtx, 1)

    S = dist_mtx

    if algorithm == "Concorde"
        tour, cost = Concorde.solve_tsp(S; kwargs...)
        shift_tour!(tour, firstcity)
        return tour, round(Int, cost)

    elseif algorithm == "LKH"
        # LKH can handle ATSP, so it doesn't require atsp2tsp
        tour, cost = LKH.solve_tsp(dist_mtx; kwargs...)
        shift_tour!(tour, firstcity)
        return tour, round(Int, cost)

    elseif algorithm == "HGS"
        if haskey(kwargs, :nbIter)
            ap = Hygese.AlgorithmParameters(;kwargs...) 
        else
            nb = round(Int, 400 * n / 100)
            ap = Hygese.AlgorithmParameters(;nbIter=nb, kwargs...) 
        end
        result = Hygese.solve_tsp(dist_mtx, ap, verbose=false)
        tour = vcat(1, result.routes[1])
        shift_tour!(tour, firstcity)
        return tour, round(Int, result.cost)

    elseif algorithm == "NearestNeighbor"
        _tour, cost = TravelingSalesmanHeuristics.nearest_neighbor(S; firstcity=firstcity, kwargs...)
        tour = _tour[1:end-1]
        return tour, round(Int, cost)

    elseif algorithm == "FarthestInsertion"
        _tour, cost = TravelingSalesmanHeuristics.farthest_insertion(S; firstcity=firstcity, kwargs...)
        tour = _tour[1:end-1]
        return tour, round(Int, cost)
        
    elseif algorithm == "CheapestInsertion"
        _tour, cost = TravelingSalesmanHeuristics.cheapest_insertion(S; firstcity=firstcity, kwargs...)
        tour = _tour[1:end-1]      
        return tour, round(Int, cost)
        
    elseif algorithm == "TwoOpt"
        if isempty(init_tour)
            init_tour = rand_tour_for_heuristics(size(S, 1); firstcity=firstcity)
        end
        
        if length(init_tour) == n 
            init_tour = shift_tour!(init_tour, firstcity)
            push!(init_tour, firstcity)
        end

        @assert length(init_tour) == n + 1
        @assert init_tour[1] == init_tour[end]

        if init_tour[1] != firstcity 
            @warn("The first city of the initial tour is not $firstcity. It will be ignored.")
        end

        _tour, cost = TravelingSalesmanHeuristics.two_opt(S, init_tour; kwargs...)
        tour = _tour[1:end-1]
        return tour, round(Int, cost)
        
    elseif algorithm == "SimulatedAnnealing"
        _tour, cost = TravelingSalesmanHeuristics.simulated_annealing(dist_mtx; kwargs...)
        tour = _tour[1:end-1]
        shift_tour!(tour, firstcity)
        return tour, round(Int, cost)
        
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