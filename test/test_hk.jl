using TSPSolvers

function generate(n)
    x = rand(n) .* 1000
    y = rand(n) .* 1000

    M = zeros(Int, n, n)
    for i in 1:n
        for j in 1:n
            if i == j   
                continue
            end
            M[i, j] = round(sqrt((x[i] - x[j])^2 + (y[i] - y[j])^2))
        end
    end
    return M
end


for k in 1:100 
    M = generate(100)
    lb = lowerbound(M)
    tour, cost = solve_tsp(M; algorithm="Concorde")
    @show lb, cost
    @assert lb <= cost
end
