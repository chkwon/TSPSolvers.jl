"""
    onetree_cost_and_degrees!(D::Matrix{Int}, S::Vector{Int}, depot::Int, pi::Vector{Float64}, deg::Vector{Int})

Compute the minimum 1-tree cost and node degrees for the Held-Karp relaxation.

A 1-tree is a spanning tree on the subset S with two edges connecting the depot to S.
Uses Prim's algorithm to compute the MST on S, then adds the two cheapest edges from depot to S.

The `!` indicates that this function modifies the `deg` vector in-place.

# Arguments
- `D::Matrix{Int}`: Distance matrix
- `S::Vector{Int}`: Subset of nodes (excluding depot)
- `depot::Int`: Depot node index
- `pi::Vector{Float64}`: Lagrangian multipliers for modified costs
- `deg::Vector{Int}`: Pre-allocated degree vector (modified in-place)

# Returns
- `total::Float64`: Total cost of the 1-tree with modified costs
"""
function onetree_cost_and_degrees!(D::Matrix{Int}, S::Vector{Int},
                                   depot::Int, pi::Vector{Float64}, deg::Vector{Int})
    # Modified cost: c̃_ij = D[i,j] + pi[i] + pi[j]

    nS = length(S)
    if nS == 0
        fill!(deg, 0)
        return 0.0
    end

    # Initialize degree array
    fill!(deg, 0)

    # 1) Compute MST on S using Prim's algorithm
    # Pre-allocate arrays for Prim's
    inT = falses(nS)
    key = Vector{Float64}(undef, nS)
    fill!(key, Inf)
    parent = Vector{Int}(undef, nS)
    fill!(parent, -1)
    key[1] = 0.0

    total = 0.0

    # Prim's algorithm with inlining
    @inbounds for iteration = 1:nS
        # Extract minimum key node not in tree
        u = 1
        best = Inf
        for i in 1:nS
            if !inT[i] && key[i] < best
                best = key[i]
                u = i
            end
        end

        inT[u] = true

        # Add edge to MST
        if parent[u] != -1
            iu = S[u]
            ip = S[parent[u]]
            # Inline modified cost calculation
            total += Float64(D[iu, ip]) + pi[iu] + pi[ip]
            deg[iu] += 1
            deg[ip] += 1
        end

        # Relax adjacent edges
        iu = S[u]
        for v in 1:nS
            if !inT[v]
                jv = S[v]
                # Inline modified cost calculation
                cij = Float64(D[iu, jv]) + pi[iu] + pi[jv]
                if cij < key[v]
                    key[v] = cij
                    parent[v] = u
                end
            end
        end
    end

    # 2) Find two cheapest edges from depot to S
    best1_node = -1
    best1_cost = Inf
    best2_node = -1
    best2_cost = Inf

    @inbounds for j in S
        # Inline modified cost calculation
        cij = Float64(D[depot, j]) + pi[depot] + pi[j]
        if cij < best1_cost
            best2_node = best1_node
            best2_cost = best1_cost
            best1_node = j
            best1_cost = cij
        elseif cij < best2_cost
            best2_node = j
            best2_cost = cij
        end
    end

    # Add depot edges to 1-tree
    if best1_node == -1 || best2_node == -1
        # Special case when S has size 1
        if best1_node != -1
            total += best1_cost
            deg[best1_node] += 1
            deg[depot] += 1
        end
    else
        total += best1_cost + best2_cost
        deg[best1_node] += 1
        deg[best2_node] += 1
        deg[depot] += 2
    end

    return total
end

"""
    heldkarp_bound(D::Matrix{Int}, S::Vector{Int};
                   depot::Int=1, max_iter::Int=100, ub::Int=typemax(Int))

Compute the Held-Karp lower bound using subgradient optimization.

The Held-Karp bound is obtained by Lagrangian relaxation of the degree constraints
in the TSP. The relaxation is solved using subgradient optimization with 1-tree computations.

# Arguments
- `D::Matrix{Int}`: Distance matrix
- `S::Vector{Int}`: Subset of nodes to visit (excluding depot)
- `depot::Int=1`: Depot node (default: 1)
- `max_iter::Int=100`: Maximum number of subgradient iterations
- `ub::Int=typemax(Int)`: Upper bound for Polyak step size

# Returns
- `bestLB::Float64`: Best lower bound found
"""
function heldkarp_bound(D::Matrix{Int}, S::Vector{Int};
                        depot::Int=1, max_iter::Int=100, ub::Int=typemax(Int))::Float64

    n = size(D, 1)
    pi = zeros(Float64, n)
    deg = zeros(Int, n)  # Pre-allocate degree vector

    bestLB = 0.0
    # Polyak step size parameter: α = θ * (ub - currentLB)/‖g‖²
    θ = 2.0

    # Pre-compute mean for diminishing step size if needed
    use_ub = (ub < typemax(Int))
    c0 = use_ub ? 0.0 : Float64(mean(D))

    # Pre-allocate subgradient vector
    g = zeros(Float64, n)

    for t in 1:max_iter
        # Compute 1-tree cost and degrees (deg is modified in-place)
        tot_mod_cost = onetree_cost_and_degrees!(D, S, depot, pi, deg)

        # Compute Held-Karp lower bound
        # Compute sum(pi[i] for i in S) efficiently
        sumpi = 0.0
        @inbounds for i in S
            sumpi += pi[i]
        end
        currLB = tot_mod_cost - 2.0 * sumpi

        if currLB > bestLB
            bestLB = currLB
        end

        # Compute subgradient: g[i] = deg[i] - 2, and gnorm2
        gnorm2 = 0.0
        @inbounds for i in S
            g[i] = Float64(deg[i] - 2)
            gnorm2 += g[i] * g[i]
        end
        g[depot] = Float64(deg[depot] - 2)
        gnorm2 += g[depot] * g[depot]

        # Check convergence
        if gnorm2 < 1e-12
            break
        end

        # Compute step size
        α = if use_ub
            # Polyak step with upper bound
            max(θ * (Float64(ub) - currLB) / gnorm2, 0.0)
        else
            # Diminishing step size
            c0 / Float64(t)
        end

        # Update Lagrangian multipliers
        @inbounds for i in S
            pi[i] += α * g[i]
        end
        pi[depot] += α * g[depot]
    end

    return bestLB
end

"""
    lowerbound(dist_mtx::AbstractMatrix{Int}; max_iter::Int=100, ub::Union{Int,Nothing}=nothing)

Compute the Held-Karp lower bound for the TSP defined by `dist_mtx`.

The Held-Karp bound is a well-known lower bound for the TSP obtained through
Lagrangian relaxation of the degree constraints. It is typically very tight
(often within a few percent of the optimal tour length).

The implementation uses subgradient optimization with 1-tree relaxations.
It assumes the tour visits all cities and returns to the starting city.

# Arguments
- `dist_mtx::AbstractMatrix{Int}`: Symmetric distance matrix (n × n), must be of Int type
- `max_iter::Int=100`: Maximum number of subgradient iterations
- `ub::Union{Int,Nothing}=nothing`: Upper bound for Polyak step size. If `nothing`, a 2-opt heuristic is used to compute one.

# Returns
- `Int`: Lower bound on the optimal tour length

# Example
```julia
M = [
    0  16   7  14;
   16   0   3   5;
    7   3   0  16;
   14   5  16   0
]
lb = lowerbound(M)
tour, cost = solve_tsp(M; algorithm="Concorde")
@assert lb <= cost

# With custom parameters
lb = lowerbound(M; max_iter=200, ub=35)
```
"""
function lowerbound(dist_mtx::Matrix{Int}; max_iter::Int=100, ub::Union{Int,Nothing}=nothing)
    n, m = size(dist_mtx)
    @assert n == m "Distance matrix must be square"
    @assert issymmetric(dist_mtx) "Distance matrix must be symmetric"

    # Visit all cities except the first (depot)
    S = collect(2:n)
    depot = 1

    # Compute upper bound using 2-opt if not provided
    ub_val = if isnothing(ub)
        _, ub_computed = solve_tsp(dist_mtx; algorithm="TwoOpt")
        ub_computed
    else
        ub
    end

    lb = heldkarp_bound(dist_mtx, S; depot=depot, max_iter=max_iter, ub=ub_val)

    return ceil(Int, lb)
end
