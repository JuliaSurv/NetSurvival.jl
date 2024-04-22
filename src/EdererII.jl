function Λₑ₂(T, Δ, age, year, rate_preds, ratetable, grid)
    # Initialize vectors: 
    num_excess   = zero(grid)
    num_pop      = zero(grid)
    den_excess   = zero(grid)
    
    # Loop over individuals
    for i in eachindex(age)
        Tᵢ = searchsortedlast(grid, T[i])
        
        rtᵢ = ratetable[rate_preds[i,:]...] # other predictors for this individual have to go here.

        for j in 1:Tᵢ
            λₚ          = daily_hazard(rtᵢ, age[i] + grid[j], year[i] + grid[j])
            ∂Λₚ         = λₚ * (grid[j+1]-grid[j])
            den_excess[j] += 1
            num_pop[j] += ∂Λₚ
        end

        num_excess[Tᵢ]   += Δ[i]  
    end
    return (num_excess ./ den_excess) .- (num_pop ./ den_excess), num_excess ./ (den_excess.^2)
end


"""
    EdererII

To call this function: 

    fit(EdererII, @formula(Surv(time,status)~covariable1 + covariable2), data, ratetable)
"""
struct EdererII <: NonparametricEstimator
    Sₑ::Vector{Float64}
    ∂Λₑ::Vector{Float64}
    σₑ::Vector{Float64}
    grid::Vector{Float64}
    function EdererII(T, Δ, age, year, rate_preds, ratetable)
        grid = mk_grid(T,1)
        ∂Λₑ, ∂σₑ = Λₑ₂(T, Δ, age, year, rate_preds, ratetable, grid)
        Sₑ = cumprod(1 .- ∂Λₑ)
        σₑ = sqrt.(cumsum(∂σₑ))
        return new(Sₑ, ∂Λₑ, σₑ, grid)
    end
end