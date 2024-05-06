abstract type NonParametricEstimator <: StatisticalModel end # maybe this one is superfluous now. 


struct NPNSEstimator{Method} <: NonParametricEstimator
    Sₑ::Vector{Float64}
    ∂Λₑ::Vector{Float64}
    ∂Λₒ::Vector{Float64}
    ∂Λₚ::Vector{Float64}
    σₑ::Vector{Float64}
    grid::Vector{Float64}
    T::Vector{Float64}
    function NPNSEstimator{Method}(T, Δ, age, year, rate_preds, ratetable) where Method
        grid = mk_grid(T,1) # precision is always 1 ? 
        ∂Λₒ, ∂Λₚ, ∂σₑ = Λ(Method, T, Δ, age, year, rate_preds, ratetable, grid)
        ∂Λₑ = ∂Λₒ .- ∂Λₚ
        Sₑ = cumprod(1 .- ∂Λₑ)
        σₑ = sqrt.(cumsum(∂σₑ))
        return new(Sₑ, ∂Λₑ, ∂Λₒ, ∂Λₚ, σₑ, grid, T)
    end
end

function mk_grid(times,prec)
    M = maximum(times)+1
    return unique(sort([(1:prec:M)..., times..., M]))
end
function Λ(::Type{M}, T, Δ, age, year, rate_preds, ratetable, grid) where M
    num_excess   = zero(grid)
    num_pop      = zero(grid)
    num_variance = zero(grid)
    den_pop      = zero(grid)
    den_excess   = zero(grid)
    Λ!(M, num_excess, den_excess, num_pop, den_pop, num_variance, T, Δ, age, year, rate_preds, ratetable, grid)
    return num_excess ./ den_excess, num_pop ./ den_pop, num_variance ./ (den_excess.^2)
end

function StatsBase.fit(::Type{E}, formula::FormulaTerm, df::DataFrame, rt::RateTables.AbstractRateTable) where {E<:NPNSEstimator}
    column_names = names(df)
    rate_predictors = String.([RateTables.predictors(rt)...])

    expected_columns = [rate_predictors...,"age","year"]
    missing_columns = filter(name -> !(name in column_names), expected_columns)
    if !isempty(missing_columns)
        throw(ArgumentError("Missing columns in data: $missing_columns"))
    end

    formula_applied = apply_schema(formula,schema(df))

    if isa(formula.rhs, ConstantTerm)
        # then there is no predictors.
        resp = modelcols(formula_applied.lhs, df)
        return E(resp[:,1], resp[:,2], df.age, df.year, select(df,rate_predictors), rt)
    else
        nms = StatsModels.termnames(formula.rhs)
        if isa(nms, String)
            pred_names = [nms]
        else
            pred_names = nms
        end

        new_df = groupby(df, pred_names)
        pp = Vector{E}()
        for i in 1:nrow(unique(df[!,pred_names]))
            resp2 = modelcols(formula_applied.lhs, new_df[i])
            push!(pp,E(resp2[:,1], resp2[:,2], new_df[i].age, new_df[i].year, select(new_df[i],rate_predictors), rt))
        end
        return pp
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
    lower_bounds = [lower[1] for lower in confint(npe; level = 0.05)]
    upper_bounds = [upper[2] for upper in confint(npe; level = 0.05)]
    df = DataFrame(Sₑ = npe.Sₑ, ∂Λₑ = npe.∂Λₑ, σₑ=npe.σₑ, lower_95_CI = lower_bounds, upper_95_CI = upper_bounds)
    show(io, df)
end
