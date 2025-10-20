module TSPSolvers

# Write your package code here.



include("main.jl")

export solve_tsp, lowerbound



# for precompilation 
M = [
    0  16   7  14;
   16   0   3   5;
    7   3   0  16;
   14   5  16   0
]
lb = lowerbound(M)
tour, cost = solve_tsp(M; algorithm="Concorde")
@assert lb <= cost


end
