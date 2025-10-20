"""
    lowerbound(dist_mtx::AbstractMatrix{Int})

Compute the Held-Karp lower bound for the TSP defined by `dist_mtx`.
The implementation uses the classic dynamic programming formulation
and assumes the tour starts and ends at the first city.
"""
function lowerbound(dist_mtx::AbstractMatrix{Int})
    nrows, ncols = size(dist_mtx)
    @assert nrows == ncols "Distance matrix must be square."
    n = nrows
    if n <= 1
        return 0
    end

    full_mask = (1 << (n - 1)) - 1
    dp = Dict{Tuple{Int, Int}, Int}()

    for j in 2:n
        mask = 1 << (j - 2)
        dp[(mask, j)] = dist_mtx[1, j]
    end

    for mask in 1:full_mask
        if count_ones(mask) <= 1
            continue
        end
        nodes = _nodes_from_mask(mask)
        for j in nodes
            prev_mask = mask & ~(1 << (j - 2))
            best = typemax(Int)
            for k in _nodes_from_mask(prev_mask)
                candidate = dp[(prev_mask, k)] + dist_mtx[k, j]
                if candidate < best
                    best = candidate
                end
            end
            @assert best < typemax(Int)
            dp[(mask, j)] = best
        end
    end

    best = typemax(Int)
    for j in 2:n
        candidate = dp[(full_mask, j)] + dist_mtx[j, 1]
        if candidate < best
            best = candidate
        end
    end

    return best
end

function _nodes_from_mask(mask::Int)
    nodes = Int[]
    bit = 0
    m = mask
    while m != 0
        if (m & 0x1) == 1
            push!(nodes, bit + 2)
        end
        bit += 1
        m >>= 1
    end
    return nodes
end


# Minimum 1-tree cost 계산:
function onetree_cost_and_degrees(D::Matrix{Float64}, S::Vector{Int},
                                  depot::Int, pi::Vector{Float64})
    # c̃_ij = D[i,j] + pi[i] + pi[j]
    
    # 1) S에 대해 MST 계산 (Prim)
    nS = length(S)
    if nS == 0
        return 0.0, Dict{Int,Int}(depot => 0)  # 엣지 없음
    end
    inT = falses(nS)
    key = fill(Inf, nS)
    parent = fill(-1, nS)
    key[1] = 0.0
    idx = Dict{Int,Int}()  # node -> idx
    for (k,v) in pairs(S); idx[v] = k; end

    function ctil(i::Int, j::Int)
        return D[i,j] + pi[i] + pi[j]
    end

    total = 0.0
    deg = Dict{Int,Int}()  # degree in 1-tree
    for v in S; deg[v] = 0; end
    deg[depot] = 0

    for _ = 1:nS
        # extract-min
        u = -1; best = Inf
        for i in 1:nS
            if !inT[i] && key[i] < best
                best = key[i]; u = i
            end
        end
        inT[u] = true
        if parent[u] != -1
            iu = S[u]; ip = S[parent[u]]
            total += ctil(iu, ip)
            deg[iu] += 1; deg[ip] += 1
        end
        # relax
        iu = S[u]
        for v = 1:nS
            if !inT[v]
                jv = S[v]
                cij = ctil(iu, jv)
                if cij < key[v]
                    key[v] = cij
                    parent[v] = u
                end
            end
        end
    end

    # 2) depot에서 S로 가는 간선 중 최저 c_ij 두 개 선택
    best1 = (node=-1, cost=Inf)
    best2 = (node=-1, cost=Inf)
    for j in S
        cij = ctil(depot, j)
        if cij < best1.cost
            best2 = best1
            best1 = (node=j, cost=cij)
        elseif cij < best2.cost
            best2 = (node=j, cost=cij)
        end
    end

    if best1.node == -1 || best2.node == -1
        # S 크기가 1일 때 등 특수 처리
        if best1.node != -1
            total += best1.cost; deg[best1.node] += 1; deg[depot] += 1
        end
    else
        total += best1.cost + best2.cost
        deg[best1.node] += 1; deg[best2.node] += 1
        deg[depot] += 2
    end

    return total, deg
end

"""
HK 하한 계산: S(=depot 제외) 방문 TSP의 Held–Karp lower bound.
- Tmax를 넘는 순간 조기 중단
"""
function hk_lower_bound(D::Matrix{Float64}, S::Vector{Int};
                        depot::Int=1, max_iter::Int=100,
                        Tmax::Real=Inf, ub::Real=Inf)

    n = size(D,1)
    pi = zeros(Float64, n)
    
    bestLB = 0.0
    # ub(상계) 있으면 Polyak step: α = θ * (ub - currentLB)/‖g‖²
    θ = 2.0

    for t in 1:max_iter
        # 1-tree on (S) with depot
        tot_mod_cost, deg = onetree_cost_and_degrees(D, S, depot, pi)

        # HK lower bound
        sumpi = 0.0
        pi[depot] = 0.0

        for i in S
            sumpi += pi[i]
        end
        
        pi[depot] = 0.0  
        sumpi = sum(pi[i] for i in S) # depot 제외
        currLB = tot_mod_cost - 2.0 * sumpi

        if currLB > bestLB
            bestLB = currLB
            if bestLB > Tmax + 1e-9
                return bestLB
            end
        end

        # subgradient g_i = deg(i) - 2
        gnorm2 = 0.0
        function g(i); return (get(deg, i, 0) - 2); end

        for i in S; gnorm2 += g(i)^2; end
        gnorm2 += (g(depot))^2

        if gnorm2 < 1e-12
            break
        end

        # step size
        if isfinite(ub)
            α = θ * (ub - currLB) / gnorm2
            α = max(α, 0.0)
        else
            # diminishing step
            c0 = mean(D) 
            α = c0 / t
        end

        # multipliers update
        for i in S
            pi[i] += α * g(i)
        end
        pi[depot] += α * g(depot)
    end
    return bestLB
end
