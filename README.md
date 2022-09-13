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
tour, tour_len = TSPSolvers.solve_tsp(M) # ([1, 3, 2, 4], 29)
```
The distance matrix M must be symmetric and integer-valued. 
It can be asymmetric for the LKH solver only. 

Also supports the TSPLIB input:
```julia
tour, tour_len = solve_tsp("gr17.tsp")
```

```julia
supported_algorithms = [
    "Concorde", 
    "LKH", 
    "HGS", 
    "NearestNeighbor", 
    "FarthestInsertion", 
    "CheapestInsertion",
    "TwoOpt", 
    "SimulatedAnnealing"
]
```
