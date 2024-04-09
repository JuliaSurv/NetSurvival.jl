

# Make the final function: 
function mk_grid(times,prec)
    M = maximum(times)+1
    return unique(sort([(1:prec:M)..., times..., M]))
end
function _Λ(T, Δ, age, year, sex, ratetable, grid)
    # Initialize vectors: 
    num_excess   = zero(grid)
    num_pop      = zero(grid)
    num_variance = zero(grid)
    den          = zero(grid)

    # Loop over individuals
    for i in eachindex(age)
        Tᵢ = searchsortedlast(grid, T[i]) # index of the time of event (death or censored) in the grid
        wₚ = 1.0
        sΛₚ = 0.0
        rtᵢ = ratetable[sex[i]] # other predictors for this individual have to go here.
        for j in 1:Tᵢ
            λₚ           = daily_hazard(rtᵢ, age[i] + grid[j], year[i] + grid[j])
            Λₚ           = λₚ * (grid[j+1]-grid[j]) # λₚ * ∂t 
            sΛₚ         += Λₚ
            wₚ           = exp(sΛₚ)
            num_pop[j] += Λₚ * wₚ
            den[j]     += wₚ
        end
        num_excess[Tᵢ]   += wₚ * Δ[i]
        num_variance[Tᵢ] += wₚ^2 * Δ[i]
    end
    return num_excess .- num_pop, num_variance, den
end
function Λ(T, Δ, age, year, sex, ratetable, grid)
    num_hazard, num_variance, den = _Λ(T, Δ, age, year, sex, ratetable, grid)
    ∂Λₑ = num_hazard ./ den
    ∂σₑ = num_variance ./ den.^2
    return ∂Λₑ, ∂σₑ
end

"""
    PoharPerme

This method was newly introduced in 2012 and it has been globally recognized since as the more effective method out of the net survival estimators.
To call this function : 

    fit(PoharPerme, @formula(Surv(time,status)~covariable1 + covariable2), data, ratetable)
"""
struct PoharPerme <: NonparametricEstimator
    Sₑ::Vector{Float64}
    ∂Λₑ::Vector{Float64}
    σₑ::Vector{Float64}
    grid::Vector{Float64}
    function PoharPerme(T, Δ, age, year, sex, ratetable)
        grid = mk_grid(T,1) # precision is always 1 ? 
        ∂Λₑ, ∂σₑ = Λ(T, Δ, age, year, sex, ratetable, grid)
        Sₑ = cumprod(1 .- ∂Λₑ)
        σₑ = sqrt.(cumsum(∂σₑ))
        return new(Sₑ, ∂Λₑ, σₑ, grid)
    end
end

"""
    EdererII

To call this function : 

    fit(EdererII, @formula(Surv(time,status)~covariable1 + covariable2), data, ratetable)
"""

struct EdererII <: NonparametricEstimator
    Sₑ::Vector{Float64}
    ∂Λₑ::Vector{Float64}
    σₑ::Vector{Float64}
    grid::Vector{Float64}
    function EdererII(T, Δ, age, year, sex, ratetable)
        grid = mk_grid(T,1) # precision is always 1 ? 
        ∂Λₑ, ∂σₑ = Λ(T, Δ, age, year, sex, ratetable, grid)
        Sₑ = cumprod(1 .- ∂Λₑ)
        σₑ = sqrt.(cumsum(∂σₑ))
        return new(Sₑ, ∂Λₑ, σₑ, grid)
    end
end

function StatsAPI.confint(npe::E; level::Real=0.05) where E <: NonparametricEstimator
    χ = sqrt(quantile(Chisq(1),1-level))
    return map(npe.Sₑ, npe.σₑ) do Sₑ,σₑ
        ci_low = exp.(log.(Sₑ) - σₑ * χ)
        ci_up = exp.(log.(Sₑ) + σₑ * χ)
        ci_low, ci_up
    end
end

function Base.show(io::IO, npe::E) where E <: NonparametricEstimator
    lower_bounds = [lower[1] for lower in confint(npe; level = 0.05)]
    upper_bounds = [upper[2] for upper in confint(npe; level = 0.05)]
    df = DataFrame(Sₑ = npe.Sₑ, ∂Λₑ = npe.∂Λₑ, σₑ=npe.σₑ, lower_95_CI = lower_bounds, upper_95_CI = upper_bounds)
    show(io, df)
end