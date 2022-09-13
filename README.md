# TSPSolvers


```julia
] add https://github.com/chkwon/TSPSolvers.jl
```


A common interface to [Concorde.jl](https://github.com/chkwon/Concorde.jl), [LKH.jl](https://github.com/chkwon/LKH.jl), [TravelingSalesmanHeuristics.jl](https://github.com/evanfields/TravelingSalesmanHeuristics.jl), and [Hygese.jl](https://github.com/chkwon/Hygese.jl).


# Example
```julia
using TSPSolvers
M = [
    0  16   7  14
   16   0   3   5
    7   3   0  16
   14   5  16   0 
]
tour, tour_len = TSPSolvers.solve_tsp(M) 
# ([1, 3, 2, 4], 29)
```
The distance matrix M must be symmetric and integer-valued. 
It can be asymmetric for the LKH solver only. 
Note that the output `tour` does not repeat the first city at the end of the tour.

Also supports the TSPLIB input:
```julia
tour, tour_len = solve_tsp("gr17.tsp")
```

You can specify the `algorithm` and `firstcity` keywords:
```julia
tour, tour_len = TSPSolvers.solve_tsp(M; algorithm="FarthestInsertion", firstcity=2)
# ([2, 4, 1, 3], 29)
```
Supported algorithms are:
```julia
supported_algorithms = [
    "Concorde",             # Concorde.jl
    "LKH",                  # LKH.jl
    "HGS",                  # Hygese.jl
    "NearestNeighbor",      # TravelingSalesmanHeuristics.jl
    "FarthestInsertion",    # TravelingSalesmanHeuristics.jl
    "CheapestInsertion",    # TravelingSalesmanHeuristics.jl
    "TwoOpt",               # TravelingSalesmanHeuristics.jl
    "SimulatedAnnealing"    # TravelingSalesmanHeuristics.jl
]
```

You may also pass the solver-specific keyword arguments. For example:
```julia
tour, tour_len = TSPSolvers.solve_tsp(M; algorithm="LKH", INITIAL_TOUR_ALGORITHM="GREEDY", RUNS=5, TIME_LIMIT=10.0)
tour, tour_len = TSPSolvers.solve_tsp(M; algorithm="FarthestInsertion", do2opt=false)
tour, tour_len = TSPSolvers.solve_tsp(M; algorithm="SimulatedAnnealing", firstcity=3, steps=10, num_starts=3)
tour, tour_len = TSPSolvers.solve_tsp(M; algorithm="HGS", nbIter=100)
```

By default, `nbIter` is set to `4 * n` for HGS.