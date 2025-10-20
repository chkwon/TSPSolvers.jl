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
