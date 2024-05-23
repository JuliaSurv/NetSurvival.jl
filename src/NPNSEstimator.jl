abstract type NPNSMethod end
struct NPNSEstimator{Method} <: StatisticalModel
    Sₑ::Vector{Float64}
    ∂Λₑ::Vector{Float64}
    ∂Λₒ::Vector{Float64}
    ∂Λₚ::Vector{Float64}
    σₑ::Vector{Float64}
    grid::Vector{Float64}
    function NPNSEstimator{Method}(T, Δ, age, year, rate_preds, ratetable) where Method<:NPNSMethod
        grid = mk_grid(T,1) # precision is always 1 ? 
        ∂Λₒ, ∂Λₚ, ∂σₑ = Λ(Method, T, Δ, age, year, rate_preds, ratetable, grid)
        ∂Λₑ = ∂Λₒ .- ∂Λₚ
        Sₑ = cumprod(1 .- ∂Λₑ)
        σₑ = sqrt.(cumsum(∂σₑ))
        return new(Sₑ, ∂Λₑ, ∂Λₒ, ∂Λₚ, σₑ, grid)
    end
end

function mk_grid(times,prec)
    M = maximum(times)
    return unique(sort([(1:prec:M)..., times..., M]))
end
function Λ(::Type{M}, T, Δ, age, year, rate_preds, ratetable, grid) where M<:NPNSMethod
    num_excess   = zero(grid)
    num_pop      = zero(grid)
    num_variance = zero(grid)
    den_pop      = zero(grid)
    den_excess   = zero(grid)
    ∂t = [diff(grid)...,1.0]
    Λ!(M, num_excess, den_excess, num_pop, den_pop, num_variance, T, Δ, age, year, rate_preds, ratetable, grid, ∂t)
    return num_excess ./ den_excess, num_pop ./ den_pop, num_variance ./ (den_excess.^2)
end

function _get_rate_predictors(rt,df)
    prd = [RateTables.predictors(rt)...]
    cl = Symbol.(names(df))
    if !(all(prd .∈ Ref(cl)) && (:age ∈ cl) && (:year ∈ cl))
        throw(ArgumentError("Missing columns in data : the chosen ratetable expects colums :age, :year and $(prd) to be present in the dataset."))
    end
    return prd
end

function StatsBase.fit(::Type{E}, formula::FormulaTerm, df::DataFrame, rt::RateTables.AbstractRateTable) where {E<:Union{NPNSEstimator, Nessie}}
    rate_predictors = _get_rate_predictors(rt,df)
    formula_applied = apply_schema(formula,schema(df))

    if isa(formula.rhs, ConstantTerm) # No predictors
        resp = modelcols(formula_applied.lhs, df)
        return E(resp[:,1], resp[:,2], df.age, df.year, select(df,rate_predictors), rt)
    else
        # we could simply group by the left side and apply fit() again, that would make sense. 

        gdf = groupby(df, StatsModels.termnames(formula.rhs))
        return rename(combine(gdf, dfᵢ -> begin
                resp2 = modelcols(formula_applied.lhs, dfᵢ)
                E(resp2[:,1], resp2[:,2], dfᵢ.age, dfᵢ.year, select(dfᵢ, rate_predictors), rt)
            end
        ), :x1 => :estimator)
    end
end

function StatsAPI.confint(npe::E; level::Real=0.05) where E <: NPNSEstimator
    χ = sqrt(quantile(Chisq(1),1-level))
    return map(npe.Sₑ, npe.σₑ) do Sₑ,σₑ
        ci_low = exp.(log.(Sₑ) - σₑ * χ)
        ci_up = exp.(log.(Sₑ) + σₑ * χ)
        ci_low, ci_up
    end
end

function Base.show(io::IO, npe::E) where E <: NPNSEstimator
    compact = get(io, :compact, false)
    if !compact
        print(io, "$(E)(t ∈ $(extrema(npe.grid))) with summary stats:\n ")
        lower_bounds = [lower[1] for lower in confint(npe; level = 0.05)]
        upper_bounds = [upper[2] for upper in confint(npe; level = 0.05)]
        df = DataFrame(Sₑ = npe.Sₑ, ∂Λₑ = npe.∂Λₑ, σₑ=npe.σₑ, lower_95_CI = lower_bounds, upper_95_CI = upper_bounds)
        show(io, df)
    else
        print(io, "$(E)(t ∈ $(extrema(npe.grid)))")
    end
end
