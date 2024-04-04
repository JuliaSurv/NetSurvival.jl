# struct for behavior
struct Surv{X, Y}
    T::X
    Δ::Y
end
struct SurvTerm{X, Y} <: AbstractTerm
    T::X
    Δ::Y
end
Base.show(io::IO, t::Surv) = print(io, string(t.T,  t.Δ == 1 ? "+" : ""))
Base.show(io::IO, t::SurvTerm) = print(io, "Surv($((t.T, t.Δ)))")

Surv(T::Symbol, Δ::Symbol) = SurvTerm(term(T), term(Δ))







