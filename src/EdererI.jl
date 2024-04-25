struct EdererIMethod end

"""
    EdererI

To call this function: 

    fit(EdererI, @formula(Surv(time,status)~covariable1 + covariable2), data, ratetable)
"""
const EdererI = NPNSEstimator{EdererIMethod}

function Λ!(::Type{EdererIMethod}, num_excess, den_excess, num_pop, den_pop, num_variance, T, Δ, age, year, rate_preds, ratetable, grid)
    Tmax= Int(maximum(T))
    for i in eachindex(age)
        Tᵢ = searchsortedlast(grid, T[i])
        Λₚ = 0.0
        rtᵢ = ratetable[rate_preds[i,:]...]
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
        num_variance[Tᵢ]   += Δ[i]    
    end
end


