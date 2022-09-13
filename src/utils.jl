

function evaluate_tour(tour::Vector{Int}, dist_mtx::Matrix{Int})
    tour_len = 0
    for i in 1:length(tour)
        if i < length(tour)
            tour_len += dist_mtx[tour[i], tour[i+1]]
        else
            tour_len += dist_mtx[tour[end], tour[1]]
        end
    end
    return tour_len
end

function dist_mat(x::Vector{T}, y::Vector{T}) where T <: Real
    @assert length(x) == length(y)
    n = length(x)
    M = Matrix{Float64}(undef, n, n)
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

