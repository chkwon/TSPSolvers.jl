using TSPSolvers
using InteractiveUtils

M = [
    0  16   7  14;
   16   0   3   5;
    7   3   0  16;
   14   5  16   0
]

D = Matrix{Int}(M)
n = size(D, 1)
S = collect(2:n)
depot = 1
max_iter = 100
ub = 35

println("Testing type stability of heldkarp_bound:")
println("=" ^ 60)

@code_warntype TSPSolvers.heldkarp_bound(D, S; depot=depot, max_iter=max_iter, ub=ub)
