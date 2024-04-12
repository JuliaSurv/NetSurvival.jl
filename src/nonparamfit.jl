#####
##### `fit` for non-parametric estimators
#####


abstract type NonparametricEstimator <: StatisticalModel end

function StatsModels.apply_schema(t::FunctionTerm{typeof(Surv)},
    sch::StatsModels.Schema,
    Mod::Type{<:Any})
    return apply_schema(SurvTerm(t.args...), sch, Mod) 
end

function StatsModels.apply_schema(t::FunctionTerm{typeof(Strata)},
    sch::StatsModels.Schema,
    Mod::Type{<:Any})
    return apply_schema(StrataTerm(t.args...), sch, Mod) 
end

function StatsModels.apply_schema(t::SurvTerm,
    sch::StatsModels.Schema,
    Mod::Type{<:Any})
    T = apply_schema(t.T, sch, Mod)
    Δ = apply_schema(t.Δ, sch, Mod)
    isa(T, ContinuousTerm) || throw(ArgumentError("Surv only works with continuous terms (got $T)"))
    isa(Δ, ContinuousTerm) ||  throw(ArgumentError("Surv only works with discrete terms (got $Δ)"))
    return SurvTerm(T, Δ)
end

function StatsModels.apply_schema(t::StrataTerm,
    sch::StatsModels.Schema,
    Mod::Type{<:Any})
    X = apply_schema(t.Covariable, sch, Mod)
    return StrataTerm(X)
end

function StatsModels.modelcols(t::SurvTerm, d::NamedTuple)
    T = modelcols(t.T, d)
    Δ = modelcols(t.Δ, d)
    return hcat(T,Δ)
end

function StatsModels.modelcols(t::StrataTerm, d::NamedTuple)
    return modelcols(t.Covariable, d)
end

function StatsBase.fit(::Type{E}, formula::FormulaTerm, df::DataFrame, rt::RateTables.AbstractRateTable) where {E<:NonparametricEstimator}
    column_names = names(df)
    rate_predictors = String.([RateTables.predictors(rt)...])

    expected_columns = [rate_predictors...,"age","year"]
    missing_columns = filter(name -> !(name in column_names), expected_columns)
    if !isempty(missing_columns)
        throw(ArgumentError("Missing columns in data: $missing_columns"))
    end

    formula = apply_schema(formula,schema(df))
    pred_names = StatsModels.termvars(formula)
    
    if nrow(unique(df[!,String.(pred_names)])) == 0 
        resp = modelcols(formula.lhs, df)
        return PoharPerme(resp[:,1], resp[:,2], df.age, df.year, select(df,rate_predictors), rt)
    else
        pp = Vector{PoharPerme}()
        new_df = groupby(df, pred_names)
        for i in 1:nrow(unique(df[!,String.(pred_names)]))      
            resp = modelcols(formula.lhs, new_df[i])
            push!(pp,PoharPerme(resp[:,1], resp[:,2], new_df[i].age, new_df[i].year, select(new_df[i],rate_predictors), rt))
        end
    end 
    return pp
end

