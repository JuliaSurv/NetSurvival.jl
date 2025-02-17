struct HakulinenMethod<:NPNSMethod end

"""
    Hakulinen

To call this function: 

    fit(Hakulinen, @formula(Surv(time,status)~covariable1 + covariable2), data, ratetable)
"""
const Hakulinen = NPNSEstimator{HakulinenMethod}

function Λ!(::HakulinenMethod, ∂Nₒ, Yₒ, ∂Nₚ, Yₚ, ∂V, T, Δ, age, year, rate_preds, ratetable, grid, ∂t)
    Date_max = maximum(T .+ year)
    for i in eachindex(age)
        Tᵢ = searchsortedlast(grid, T[i])
        if Δ[i] == 1
            # T2 = searchsortedlast(grid,Tmax-year[i])-1
            T2 = searchsortedlast(grid,Date_max-year[i])-1
        else
            T2 = Tᵢ 
            # T2 = searchsortedlast(grid,Tmax)-1
        end
        Λₚ = 0.0
        rtᵢ = ratetable[rate_preds[i,:]...]
        for j in 1:T2
            λₚ          = daily_hazard(rtᵢ, age[i] + grid[j], year[i] + grid[j])
            ∂Λₚ         = λₚ * ∂t[j]
            Λₚ         += ∂Λₚ
            Sₚ          = exp(-Λₚ)
            ∂Nₚ[j] += (Sₚ * ∂Λₚ)
            Yₚ[j] += Sₚ
        end
        for j in 1:Tᵢ
            Yₒ[j] += 1
        end
        ∂Nₒ[Tᵢ]   += Δ[i]  
        ∂V[Tᵢ]   += Δ[i]  
    end
end