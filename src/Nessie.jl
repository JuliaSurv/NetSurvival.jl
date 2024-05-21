struct Nessie
    expected_sample_size::Vector{Float64}
    expected_life_time::Float64
    grid::Vector{Float64}
    function Nessie(T, Δ, age, year, rate_preds, ratetable)
        annual_grid = 0:RateTables.RT_DAYS_IN_YEAR:maximum(T)
        exp_spl_size = zeros(length(annual_grid))
        exp_life_time = 0.0
        for i in eachindex(age)
            Lᵢ = Life(ratetable[rate_preds[i,:]...], age[i], year[i])
            for j in eachindex(annual_grid)
                exp_spl_size[j] += ccdf(Lᵢ, annual_grid[j])
            end
            exp_life_time += expectation(Lᵢ)
        end
        return new(exp_spl_size, exp_life_time / RateTables.RT_DAYS_IN_YEAR / length(age), annual_grid)
    end
end

"""
    nessie(formula, data, ratetable)

The Nessie function estimates the sample size by yearly intervals as well as averages an estimated lifespan left for a given group.  
"""
function nessie(args...)
    r = fit(Nessie,args...)
    if (typeof(r)<:Nessie)
        return r
    end
    transform!(r, :estimator => ByRow(x-> (x.grid, x.expected_life_time, x.expected_sample_size)) => [:grid, :expected_life_time,:expected_sample_size])
    select!(r, Not(:estimator))

    lt = deepcopy(r)
    select!(lt, Not([:expected_sample_size, :grid]))

    select!(r, Not(:expected_life_time))
    return lt, r
end

# Maybe not necessary ? No need to clutter the interface too much.. 
expected_life_time(x::Nessie) = x.expected_life_time
expected_sample_size(x::Nessie) = x.expected_sample_size
