struct CrudeMortality
    Λₑ::Vector{Float64}
    Λₚ::Vector{Float64}
end

function StatsBase.fit(::Type{E}, formula::FormulaTerm, df::DataFrame, rt::RateTables.AbstractRateTable) where {E <: CrudeMortality}
    formula_applied = apply_schema(formula,schema(df))
    resp = modelcols(formula_applied.lhs, df)

    rate_preds = select(df, String.([RateTables.predictors(rt)...]))

    grid            = mk_grid(resp[:,1],1)
    num_pop         = zero(grid)
    den_pop         = zero(grid)
    num_excess      = zero(grid)
    den_excess      = zero(grid)
    num_variance    = zero(grid)

    Λ!(EdererIIMethod, num_excess, den_excess, num_pop, den_pop, num_variance, resp[:,1], resp[:,2], df.age, df.year, rate_preds, rt, grid)
    
    ∂λₒ = num_excess ./ den_excess
    Sₒ = cumprod(1 .- ∂λₒ)

    causeSpec = zero(grid)
    population = zero(grid)

    for i in 1:nrow(df)
        Tᵢ = searchsortedlast(grid, resp[:,1][i])
        Λₑ = 0.0
        ∂λₚ = 0.0

        for j in 1:Tᵢ
            ∂λₑ = (num_excess[j] ./ den_excess[j]) .- (num_pop[j] ./ den_pop[j])
            Λₑ        += ∂λₑ * Sₒ[j]
            ∂λₚ = num_pop[j] / den_pop[j]
            causeSpec[j+1] = Λₑ
            population[j+1] = population[j] + ∂λₚ * Sₒ[j]
        end
    end

    return CrudeMortality(causeSpec, population)
end

function Base.show(io::IO, crud::CrudeMortality) 
    df = DataFrame(Λₑ = crud.Λₑ, Λₚ = crud.Λₚ)
    show(io, df)
end