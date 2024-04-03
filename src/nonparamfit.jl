#####
##### `fit` for non-parametric estimators
#####


abstract type NonparametricEstimator <: StatisticalModel end 


function StatsModels.apply_schema(t::FunctionTerm{typeof(Surv)},
    sch::StatsModels.Schema,
    Mod::Type{<:Any})
return apply_schema(SurvTerm(t.args...), sch, Mod)
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

function StatsModels.modelcols(t::SurvTerm, d::NamedTuple)
T = modelcols(t.T, d)
Δ = modelcols(t.Δ, d)
return hcat(T,Δ) ############ <<<<<----
end

function StatsBase.fit(::Type{E}, formula::FormulaTerm, df::DataFrame, rt::RateTableV2) where {E<:NonparametricEstimator}
    column_names = names(df)
    expected_columns = String.(keys(rt.axes))

    missing_columns = filter(name -> !(name in column_names), expected_columns)
    if !isempty(missing_columns)
        throw(ArgumentError("Missing columns in data: $missing_columns"))
    end

    formula = apply_schema(formula,schema(df))
    pred_names = StatsModels.termvars(formula)
    new_df = groupby(df, pred_names)
    pred_names = String.(pred_names)

    g = term(1) ~ foldl(+,term.(keys(rt.axes)))
    g = apply_schema(g,schema(df))
    pred_g = modelcols(g.rhs, df)

    temp = [x == 1 ? :male : :female for x in pred_g[:,3]]

    pp = Vector{PoharPerme}()

    if nrow(unique(df[!,pred_names])) == 0 
        resp = modelcols(formula.lhs, df)
        pred_g = modelcols(g.rhs, df)
        push!(pp,PoharPerme(resp[:,1], resp[:,2], pred_g[:,1], pred_g[:,2], temp, rt))
    else
        for i in 1:nrow(unique(df[!,pred_names]))        
            resp = modelcols(formula.lhs, new_df[i])
            pred_g = modelcols(g.rhs, new_df[i])
            push!(pp,PoharPerme(resp[:,1], resp[:,2], pred_g[:,1], pred_g[:,2], temp, rt))
        end
    end 
    return pp
end

