# struct for behavior
# struct Surv{X, Y}
#     T::X
#     Δ::Y
# end
struct SurvTerm{X, Y} <: AbstractTerm
    T::X
    Δ::Y
end
# Base.show(io::IO, t::Surv) = print(io, string(t.T,  t.Δ == 1 ? "+" : ""))
Base.show(io::IO, t::SurvTerm) = print(io, "Surv($((t.T, t.Δ)))")

Surv(T::Symbol, Δ::Symbol) = SurvTerm(term(T), term(Δ))
Surv(T::Float64,Δ::Bool) = (T,Δ)

Strata(x) = x
struct StrataTerm{X} <: AbstractTerm
    Covariable::X
end
Base.show(io::IO, t::StrataTerm) = print(io, "Strata($((t.Covariable)))")

Strata(Covariables::Vector) = StrataTerm(term(Covariables))

StatsModels.termvars(p::StrataTerm) = StatsModels.termvars(p.Covariable)

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






