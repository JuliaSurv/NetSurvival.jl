abstract type NPNSMethod end
struct NPNSEstimator{Method} <: StatisticalModel
    Sₑ::Vector{Float64}
    ∂Λₑ::Vector{Float64}
    ∂Λₒ::Vector{Float64}
    ∂Λₚ::Vector{Float64}
    σₑ::Vector{Float64}
    grid::Vector{Float64}
    function NPNSEstimator(m::Method, T, Δ, age, year, rate_preds, ratetable) where Method<:NPNSMethod
        grid = mk_grid(T,1) # precision is always 1 ? 
        ∂Λₒ, ∂Λₚ, ∂σₑ = Λ(m, T, Δ, age, year, rate_preds, ratetable, grid)
        ∂Λₑ = ∂Λₒ .- ∂Λₚ
        Sₑ = clamp.(cumprod(1 .- ∂Λₑ),0,Inf)
        σₑ = sqrt.(cumsum(∂σₑ))
        return new{Method}(Sₑ, ∂Λₑ, ∂Λₒ, ∂Λₚ, σₑ, grid)
    end
end
NPNSEstimator{Method}() where Method = Method()
NPNSEstimator{Method}(T, Δ, age, year, rate_preds, ratetable) where Method<:NPNSMethod = NPNSEstimator(Method(), T, Δ, age, year, rate_preds, ratetable)
NPNSEstimator(::Type{Method}, T, Δ, age, year, rate_preds, ratetable) where Method<:NPNSMethod = NPNSEstimator(Method(), T, Δ, age, year, rate_preds, ratetable)

function mk_grid(times,prec)
    M = maximum(times)
    return unique(sort([(1:prec:M)..., times..., M]))::Vector{Float64}
end
function Λ(m::M, T, Δ, age, year, rate_preds, ratetable, grid) where M<:NPNSMethod
    ∂Nₒ   = zero(grid)
    ∂Nₚ      = zero(grid)
    ∂V = zero(grid)
    Yₚ      = zero(grid)
    Yₒ   = zero(grid)
    ∂t = [diff(grid)...,1.0]
    Λ!(m, ∂Nₒ, Yₒ, ∂Nₚ, Yₚ, ∂V, T, Δ, age, year, rate_preds, ratetable, grid, ∂t)
    return ∂Nₒ ./ Yₒ, ∂Nₚ ./ Yₚ, ∂V ./ (Yₒ.^2)
end

function _get_rate_predictors(rt,df)
    prd = [RateTables.predictors(rt)...]
    cl = Symbol.(names(df))
    if !(all(prd .∈ Ref(cl)) && (:age ∈ cl) && (:year ∈ cl))
        throw(ArgumentError("Missing columns in data : the chosen ratetable expects colums :age, :year and $(prd) to be present in the dataset."))
    end
    return prd
end

function StatsBase.fit(::Type{NPNSEstimator{M}}, formula::FormulaTerm, df, rt) where {M<:NPNSMethod}
    return StatsBase.fit(NPNSEstimator{M}(), formula, df, rt)
end
function StatsBase.fit(m::M, formula::FormulaTerm, df::DataFrame, rt::RateTables.AbstractRateTable) where {M<:NPNSMethod}
    rate_predictors = _get_rate_predictors(rt,df)
    formula_applied = apply_schema(formula,schema(df))

    if isa(formula.rhs, ConstantTerm) # No predictors
        resp = modelcols(formula_applied.lhs, df)::Matrix{Float64}
        return NPNSEstimator(m, resp[:,1], resp[:,2], df.age, df.year, select(df,rate_predictors), rt)
    else
        gdf = groupby(df, StatsModels.termnames(formula.rhs))
        return rename(combine(gdf, dfᵢ -> begin
                resp2 = modelcols(formula_applied.lhs, dfᵢ)::Matrix{Float64}
                NPNSEstimator(m, resp2[:,1], resp2[:,2], dfᵢ.age, dfᵢ.year, select(dfᵢ, rate_predictors), rt)
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

function (S::NPNSEstimator)(t)
    i_t = findlast(S.grid .<= t)
    isnothing(i_t) && return one(t)
    return S.Sₑ[i_t]
end
function variance(S::NPNSEstimator, t)
    i_t = findlast(S.grid .<= t)
    isnothing(i_t) && return zero(t)
    return S.σₑ[i_t]
end