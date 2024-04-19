function Λₑ₁(T, Δ, age, year, rate_preds, ratetable, grid)
    # Initialize vectors: 
    num_excess   = zero(grid)
    num_pop      = zero(grid)
    num_variance = zero(grid)
    den          = zero(grid)
    Yᵢ           = length(T)
    
    # Loop over individuals
    for i in eachindex(age)
        Tᵢ = searchsortedlast(grid, T[i])
        sΛₚ = 0.0
        Yᵢ           = length(T)
        
        rtᵢ = ratetable[rate_preds[i,:]...] # other predictors for this individual have to go here.
        for j in Tᵢ
            λₚ           = daily_hazard(rtᵢ, age[i] + grid[i], year[i] + grid[i])
            Λₚ           = λₚ * (grid[j+1]-grid[j]) 
            sΛₚ         += Λₚ
            wₚ           = exp(Λₚ)
            den[j]      += wₚ
            num_pop[j]  = (den[j] * sΛₚ) / den[j]
            Yᵢ -= sum(Δ[i]==1)
        end
        
        num_excess[Tᵢ]   += Δ[i] / Yᵢ
        num_variance[Tᵢ] += Δ[i] / Yᵢ^2 
    end
    return num_excess .- num_pop, num_variance
end

"""
    EdererI

To call this function: 

    fit(EdererI, @formula(Surv(time,status)~covariable1 + covariable2), data, ratetable)
"""
struct EdererI <: NonparametricEstimator
    Sₑ::Vector{Float64}
    ∂Λₑ::Vector{Float64}
    σₑ::Vector{Float64}
    grid::Vector{Float64}
    function EdererI(T, Δ, age, year, rate_preds, ratetable)
        grid = mk_grid(T,1)
        ∂Λₑ, ∂σₑ = Λₑ₁(T, Δ, age, year, rate_preds, ratetable, grid)
        Sₑ = cumprod(1 .- ∂Λₑ)
        σₑ = sqrt.(cumsum(∂σₑ))
        return new(Sₑ, ∂Λₑ, σₑ)
    end
end