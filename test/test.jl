using TravelingSalesmanHeuristics
using LKH
using Infiltrator
using Statistics 
using Random


function dist_mat(x, y)
    @assert length(x) == length(y)
    n = length(x)
    M = zeros(n, n)
    for i in 1:n 
        for j in 1:n 
            M[i, j] = sqrt((x[i]-x[j])^2 + (y[i]-y[j])^2)
        end
    end

    return M
end

function rand_tour(n; firstcity=1)
    v = shuffle(1:n)
    idx = findfirst(x -> x==firstcity, v)
    circshift!(v, idx - 1)
    push!(v, firstcity)
    return v
end

function test(n)
    n_rep = 100 

    lkh_sol = zeros(n_rep)
    far_sol = zeros(n_rep)
    far2_sol = zeros(n_rep)
    two_sol = zeros(n_rep)
    sim_sol = zeros(n_rep)

    lkh_time = 0.0
    far_time = 0.0
    far2_time = 0.0
    two_time = 0.0
    sim_time = 0.0

    for i in 1:n_rep
        x = rand(n) .* 1000
        y = rand(n) .* 1000
        dist = round.(Int, dist_mat(x, y))

        far_time += @elapsed tour, far_sol[i] = farthest_insertion(dist; firstcity=1, do2opt=false)
        far2_time += @elapsed tour, far2_sol[i] = farthest_insertion(dist; firstcity=1, do2opt=true)
        two_time += @elapsed tour, two_sol[i] = two_opt(dist, rand_tour(n, firstcity=1))
        lkh_time += @elapsed tour, lkh_sol[i] = LKH.solve_tsp(dist)    
        sim_time += @elapsed tour, sim_sol[i] = simulated_annealing(dist)
    end

    println("-------- [ n = $n ] --------------------------------------")
    println("FAR : time = $(far_time), avg cost = $(mean(far_sol))")
    println("FAR2: time = $(far2_time), avg cost = $(mean(far2_sol))")
    println("TWO : time = $(two_time), avg cost = $(mean(two_sol))")
    println("LKH : time = $(lkh_time), avg cost = $(mean(lkh_sol))")
    println("SIM : time = $(sim_time), avg cost = $(mean(sim_sol))")

end

test(10)
test(100)
