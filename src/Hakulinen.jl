function Λₕₖ(T, Δ, age, year, rate_preds, ratetable, grid)
    # Initialize vectors: 
    num_excess   = zero(grid)
    num_pop      = zero(grid)
    den_pop      = zero(grid)
    den_excess   = zero(grid)
    
    Tmax= Int(maximum(T))
    # Loop over individuals
    for i in eachindex(age)
        Tᵢ = searchsortedlast(grid, T[i]) # index of the time of event (death or censored) in the grid
        wₚ = 1.0
        sΛₚ = 0.0
        rtᵢ = ratetable[rate_preds[i,:]...] # other predictors for this individual have to go here.
        for j in 1:Tᵢ
            λₚ           = daily_hazard(rtᵢ, age[i] + grid[j], year[i] + grid[j])
            Λₚ           = λₚ * (grid[j+1]-grid[j]) # λₚ * ∂t 
            sΛₚ         += Λₚ
            wₚ           = exp(sΛₚ)
            num_pop[j] += Λₚ * wₚ
            den_pop[j]    += wₚ
            den_excess[j] += 1
        end
        num_excess[Tᵢ]   += Δ[i]  
    end
    return (num_excess ./ den_excess) .- (num_pop ./ den_pop), num_excess ./ (den_excess.^2)
end


"""
    Hakulinen

To call this function: 

    fit(Hakulinen, @formula(Surv(time,status)~covariable1 + covariable2), data, ratetable)
"""
struct Hakulinen <: NonparametricEstimator
    Sₑ::Vector{Float64}
    ∂Λₑ::Vector{Float64}
    σₑ::Vector{Float64}
    grid::Vector{Float64}
    function Hakulinen(T, Δ, age, year, rate_preds, ratetable)
        grid = mk_grid(T,1)
        ∂Λₑ, ∂σₑ = Λₕₖ(T, Δ, age, year, rate_preds, ratetable, grid)
        Sₑ = cumprod(1 .- ∂Λₑ)
        σₑ = sqrt.(cumsum(∂σₑ))
        return new(Sₑ, ∂Λₑ, σₑ)
    end
end
