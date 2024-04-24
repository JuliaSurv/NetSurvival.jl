function Λₕₖ(T, Δ, age, year, rate_preds, ratetable, grid)
    # Initialize vectors: 
    num_excess   = zero(grid)
    den_excess   = zero(grid)
    num_pop      = zero(grid)
    den_pop      = zero(grid)
    Tmax= Int(maximum(T))
    Date_max = maximum(T .+ year)
    # Loop over individuals
    for i in eachindex(age)
        Tᵢ = searchsortedlast(grid, T[i]) # index of the time of event (death or censored) in the grid
        if Δ[i] == 1
            # T2 = searchsortedlast(grid,Tmax-year[i])-1
            T2 = searchsortedlast(grid,Date_max-year[i])-1
        else
            T2 = Tᵢ 
            # T2 = searchsortedlast(grid,Tmax)-1
        end
        Λₚ = 0.0
        rtᵢ = ratetable[rate_preds[i,:]...] # other predictors for this individual have to go here.
        for j in 1:T2
            λₚ          = daily_hazard(rtᵢ, age[i] + grid[j], year[i] + grid[j])
            ∂Λₚ         = λₚ * (grid[j+1]-grid[j]) 
            Λₚ         += ∂Λₚ
            Sₚ          = exp(-Λₚ)
            num_pop[j] += (Sₚ * ∂Λₚ)
            den_pop[j] += Sₚ
        end
        for j in 1:Tᵢ
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
        return new(Sₑ, ∂Λₑ, σₑ, grid)
    end
end
