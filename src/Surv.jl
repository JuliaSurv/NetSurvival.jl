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






