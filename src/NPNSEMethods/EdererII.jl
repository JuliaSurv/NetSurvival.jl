struct EdererIIMethod<:NPNSMethod end

"""
    EdererII

This method estimates net survival probabilities by applying the following estimation:

```math
\\partial\\hat{\\Lambda}_E(t) = \\frac{\\sum_i N_i(u)}{\\sum_i Y_i(u)} - \\frac{\\sum_i Y_i(u)d\\Lambda_{P_i}(u)}{\\sum_i Y_i(u)}
```

To call this function: 

    fit(EdererII, @formula(Surv(time,status)~covariable1 + covariable2), data, ratetable)
"""
const EdererII = NPNSEstimator{EdererIIMethod}

function Λ!(::EdererIIMethod, ∂Nₒ, Yₒ, ∂Nₚ, Yₚ, ∂V, T, Δ, age, year, rate_preds, ratetable, grid, ∂t)
    for i in eachindex(age)
        Tᵢ = searchsortedlast(grid, T[i])
        rtᵢ = ratetable[rate_preds[i,:]...]
        for j in 1:Tᵢ
            λₚ          = daily_hazard(rtᵢ, age[i] + grid[j], year[i] + grid[j])
            ∂Λₚ         = λₚ * ∂t[j]
            Yₒ[j] += 1
            Yₚ[j] += 1
            ∂Nₚ[j] += ∂Λₚ
        end
        ∂Nₒ[Tᵢ]   += Δ[i]  
        ∂V[Tᵢ]   += Δ[i]  
    end
end

