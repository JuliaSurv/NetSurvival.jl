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

function Λ!(::Type{PoharPermeMethod}, ∂Nₑ, Yₑ, ∂Nₚ, Yₚ, ∂V, T, Δ, a, d, rate_preds, ratetable, grid, ∂t)
    for i in eachindex(age)
        Tᵢ = searchsortedlast(grid, T[i])
        Λₚ, wₚ, rtᵢ = 0.0, 1.0, ratetable[rate_preds[i,:]...]
        for j in 1:Tᵢ
            λₚ      = daily_hazard(rtᵢ, a[i] + grid[j], d[i] + grid[j])
            ∂Λₚ     = λₚ * ∂t[j]
            Λₚ     += ∂Λₚ
            wₚ      = exp(Λₚ)
            ∂Nₚ[j] += ∂Λₚ * wₚ
            Yₑ[j] += wₚ
        end
        ∂Nₑ[Tᵢ] += wₚ * Δ[i]
        ∂V[Tᵢ]  += wₚ^2 * Δ[i]
    end
    Yₚ .= Yₑ
end