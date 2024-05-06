"""
    CrudeMortality

This method calculates the crude mortality and presents both the excess mortality and population mortality rates.

To apply the Cronin-Feuer estimator to your data based on a certain rate table, apply the example below to your code : 

    CrudeMortality(EdererII, data, ratetable)

To check for the other methods, simply replace "EdererII" with the method of your choice (PoharPerme, EdererI, Hakulinen).

References: 
* [cronin2000cumulative](@cite) Cronin, Kathleen A and Feuer, Eric J (2000). Cumulative cause-specific mortality for cancer patients in the presence of other causes: a crude analogue of relative survival.
"""
struct CrudeMortality
    Λₒ::Vector{Float64}
    Λₑ::Vector{Float64}
    Λₚ::Vector{Float64}
    function CrudeMortality(npe::Method, df::DataFrame) where Method
        Λₑ = zero(npe.grid)
        Λₚ = zero(npe.grid)
        
        Sₒ = cumprod(1 .- npe.∂Λₒ)

        for i in 1:nrow(df)
            ∂Λₑ = 0.0 
            Tᵢ = searchsortedlast(npe.grid, df.time[i])
            for j in 1:Tᵢ
                ∂Λₑ += npe.∂Λₑ[j] * Sₒ[j]
                Λₑ[j+1] = ∂Λₑ
                Λₚ[j+1] = Λₚ[j] + npe.∂Λₚ[j]*Sₒ[j]
            end
        end
        return new(Λₑ .- Λₚ, Λₑ, Λₚ)
    end
end

StatsBase.fit(::Type{CrudeMortality}, args...) = CrudeMortality(args...) 

function Base.show(io::IO, crud::CrudeMortality) 
    df = DataFrame(Λₒ = crud.Λₒ, Λₑ = crud.Λₑ, Λₚ = crud.Λₚ)
    show(io, df)
end