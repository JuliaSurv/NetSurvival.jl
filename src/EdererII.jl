struct EdererIIMethod end

"""
    EdererII

To call this function: 

    fit(EdererII, @formula(Surv(time,status)~covariable1 + covariable2), data, ratetable)
"""
const EdererII = NPNSEstimator{EdererIIMethod}

function Λ!(::Type{EdererIIMethod}, num_excess, den_excess, num_pop, den_pop, num_variance, T, Δ, age, year, rate_preds, ratetable, grid)
    for i in eachindex(age)
        Tᵢ = searchsortedlast(grid, T[i])
        rtᵢ = ratetable[rate_preds[i,:]...]
        for j in 1:Tᵢ
            λₚ          = daily_hazard(rtᵢ, age[i] + grid[j], year[i] + grid[j])
            ∂Λₚ         = λₚ * (grid[j+1]-grid[j])
            den_excess[j] += 1
            den_pop[j] += 1
            num_pop[j] += ∂Λₚ
        end
        num_excess[Tᵢ]   += Δ[i]  
        num_variance[Tᵢ]   += Δ[i]  
    end
end
