struct PoharPermeMethod<:NPNSMethod end

"""
    PoharPerme

This method estimates net survival probabilities by applying the following estimation:

```math
\\partial\\hat{\\Lambda}_E(t) = \\frac{\\sum_i \\frac{dN_i(u)}{S_{P_i}(u)} - \\sum_i \\frac{Y_i(u)}{S_{P_i}(u)}d\\Lambda_{P_i}(u)}{\\sum_i \\frac{Y_i(u)}{S_{P_i}(u)}}
```

To fit the Pohar Perme to your data based on a certain rate table, apply the example below to your code : 

    fit(PoharPerme, @formula(Surv(time,status)~covariable1 + covariable2), data, ratetable)

References: 
* [PoharPerme2012](@cite) Perme, Maja Pohar and Stare, Janez and Estève, Jacques (2012). On Estimation in Relative Survival.
"""
const PoharPerme = NPNSEstimator{PoharPermeMethod}

function Λ!(::Type{PoharPermeMethod}, num_excess, den_excess, num_pop, den_pop, num_variance, T, Δ, age, year, rate_preds, ratetable, grid, ∂t)
    for i in eachindex(age)
        Tᵢ = searchsortedlast(grid, T[i])
        Λₚ = 0.0
        wₚ = 1.0
        rtᵢ = ratetable[rate_preds[i,:]...]
        for j in 1:Tᵢ
            λₚ          = daily_hazard(rtᵢ, age[i] + grid[j], year[i] + grid[j])
            ∂Λₚ         = λₚ * ∂t[j]#(grid[j+1]-grid[j]) # λₚ * ∂t 
            Λₚ         += ∂Λₚ
            wₚ          = exp(Λₚ)
            num_pop[j] += ∂Λₚ * wₚ
            den_excess[j] += wₚ
        end
        num_excess[Tᵢ]   += wₚ * Δ[i]
        num_variance[Tᵢ] += wₚ^2 * Δ[i]
    end
    den_pop .= den_excess
end