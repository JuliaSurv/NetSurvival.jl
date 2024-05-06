struct CrudeMortality
    Λₑ::Vector{Float64}
    Λₚ::Vector{Float64}
end

function StatsBase.fit(::Type{E}, formula::FormulaTerm, df::DataFrame, rt::RateTables.AbstractRateTable) where {E <: CrudeMortality}
    formula_applied = apply_schema(formula,schema(df))
    resp = modelcols(formula_applied.lhs, df)

    rate_preds = select(df, String.([RateTables.predictors(rt)...]))

    eII = NPNSEstimator{EdererIIMethod}(resp[:,1], resp[:,2], df.age, df.year, rate_preds, rt)
    num_excess      = zero(eII.grid)
    den_excess      = zero(eII.grid)
    Sₑ = 0.0
    Sₚ = 0.0 

    for i in 1:nrow(df)
        Tᵢ = searchsortedlast(eII.grid, resp[:,1][i])
        for j in 1:Tᵢ
            den_excess[j] += 1
        end
        num_excess[Tᵢ] += resp[:,2][i]
    end

    ∂λₒ = num_excess ./ den_excess
    Sₒ = cumprod(1 .- ∂λₒ)

    causeSpec = zero(eII.grid)
    population = zero(eII.grid)
    num_pop = zero(eII.grid)
    den_pop = zero(eII.grid)

    for i in 1:nrow(df)
        Tᵢ = searchsortedlast(eII.grid, resp[:,1][i])
        Λₚ = 0.0
        Λₑ = 0.0
        rtᵢ = rt[rate_preds[i,:]...]
        for j in 1:Tᵢ
            λₚ          = daily_hazard(rtᵢ, df.age[i] + eII.grid[j], df.year[i] + eII.grid[j])
            ∂Λₚ         = λₚ * (eII.grid[j+1]-eII.grid[j])
            den_pop[j] += 1
            num_pop[j] += ∂Λₚ
            Λₑ        += eII.∂Λₑ[j] * Sₒ[j]
            causeSpec[j] = Λₑ
            population[j] += (num_pop[j]/den_pop[j]) * Sₒ[j]
        end
    end

    return CrudeMortality(causeSpec, population)
end

function Base.show(io::IO, crud::CrudeMortality) 
    df = DataFrame(Λₑ = crud.Λₑ, Λₚ = crud.Λₚ)
    show(io, df)
end