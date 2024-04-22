function Λₑ₁(T, Δ, age, year, rate_preds, ratetable, grid)
    # Initialize vectors: 
    num_excess   = zero(grid)
    num_pop      = zero(grid)
    den_pop      = zero(grid)
    den_excess   = zero(grid)
    
    Tmax= Int(maximum(T))
    # Loop over individuals
    for i in eachindex(age)
        Tᵢ = searchsortedlast(grid, T[i])
        Λₚ = 0.0
        
        rtᵢ = ratetable[rate_preds[i,:]...] # other predictors for this individual have to go here.
        for j in 1:Tmax
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
        return new(Sₑ, ∂Λₑ, σₑ, grid)
    end
end