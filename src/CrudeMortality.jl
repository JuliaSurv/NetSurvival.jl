"""
    CrudeMortality

This method calculates the crude mortality and presents both the excess mortality and population mortality rates.

The default Cronin-Feuer estimator can be fitted to data with the following interface: 

    fit(CrudeMortality, args...)

where the `args` are passed to `fit(EdererII,args...)` to compute the excess hazard.

A more direct syntax can be used, specifying directly the estimator for the excess hazard:

    CrudeMortality(fit(EdererII,args...))

To use another excess hazard estimator, simply replace `EdererII` with the method of your choice (`PoharPerme`, `EdererI`, `Hakulinen`).

References: 
* [cronin2000cumulative](@cite) Cronin, Kathleen A and Feuer, Eric J (2000). Cumulative cause-specific mortality for cancer patients in the presence of other causes: a crude analogue of relative survival.
"""
struct CrudeMortality
    Mₒ::Vector{Float64}
    Mₑ::Vector{Float64}
    Mₚ::Vector{Float64}
    function CrudeMortality(npe::NPNSEstimator{Method}) where Method
        Sₒ = cumprod(1 .- npe.∂Λₒ)
        Mₑ = cumsum(npe.∂Λₑ .* Sₒ)
        Mₚ = cumsum(npe.∂Λₚ .* Sₒ)
        return new(Mₑ .+ Mₚ, Mₑ, Mₚ)
    end
end

StatsBase.fit(::Type{CrudeMortality}, args...) = CrudeMortality(fit(EdererII, args...)) 

function Base.show(io::IO, crud::CrudeMortality) 
    df = DataFrame(Mₒ = crud.Mₒ, Mₑ = crud.Mₑ, Mₚ = crud.Mₚ)
    show(io, df)
end