struct GenPoharPermeMethod{CopulaType} <:NPNSMethod 
    C::CopulaType
    function GenPoharPermeMethod(C::Copulas.Copula{2})
        # if isa(C, IndependentCopula)
        #     return PoharPermeMethod()
        # else
            return new{typeof(C)}(C)
        # end
    end
end
GenPoharPermeMethod() = GenPoharPermeMethod(IndependentCopula(2))

"""
    GenPoharPerme

This method estimates net survival probabilities under the generalized pohar perme method, taking into account a copula for (E,P) as described in the reference paper. You may use it as: 

    fit(GenPoharPerme(copula), @formula(Surv(time, status)~ x1 + x2), data, ratetable)

Note that `GenPoharPerme(::IndependenceCopula)` simply returns `PoharPerme`.

References: 
* [Laverny2024](@cite) Laverny, O,  Grafféo, N and Giorgi, R: Non parametric estimation of net survival under dependence assumption. Working paper. 
"""
const GenPoharPerme = NPNSEstimator{GenPoharPermeMethod}
GenPoharPerme(C::Copulas.Copula{2}) = GenPoharPermeMethod(C::Copulas.Copula{2})

function Λ!(m::GenPoharPermeMethod{TC},  ∂Nₒ, Yₒ, ∂Nₚ, Yₚ, ∂V, T, Δ, age, year, rate_preds, ratetable, grid, ∂t) where {TC}

    I, J = length(T), length(grid) 
    ∂Nⱼ, Yⱼ, ∂Sₚ, Sₚ  = zeros(I), trues(I), zeros(I), ones(I) # size I = nb of indivs
    ∂Λₚ, ∂Λₑ, Sₑ = 0.0, 0.0, 1.0
    rt = [ratetable[rate_preds[i,:]...] for i in 1:I]::Vector{RateTables.BasicRateTable}

    for j in 1:J
        tⱼ, ∂tⱼ = grid[j], ∂t[j]
        for i in 1:I
            Yⱼ[i] = T[i] >= tⱼ
            if Yⱼ[i]
                ∂Nⱼ[i]  = T[i] == tⱼ && Δ[i]
                ∂Λₚ = ∂tⱼ * daily_hazard(rt[i], age[i] + tⱼ, year[i] + tⱼ)
                Sₚ[i]   = Sₚ[i] * exp(- ∂Λₚ) 
                ∂Sₚ[i] = - ∂Λₚ * Sₚ[i]
            end
        end
        # Maybe we could pop out elements of the vectors when they get out instead orf re-filtering each time ? Would be faster ? 
        # If observations are sorted, might indeed be faster ! 
        ∂Λₑ             = _mk_∂Λₑ(m.C, Sₑ, Sₚ[Yⱼ], ∂Sₚ[Yⱼ], ∂Nⱼ[Yⱼ])
        Sₑ              = clamp(Sₑ * exp(- ∂Λₑ), 0, 1)
        ∂V[j], Yₒ[j], r = _eq_var(m.C, Sₑ, Sₚ[Yⱼ], ∂Nⱼ[Yⱼ])
        ∂Nₚ[j]          = ∂Λₑ * Yₒ[j] - r
        ∂Nₒ[j]          = ∂Λₑ * Yₒ[j] + ∂Nₚ[j]
        Yₚ[j]           = Yₒ[j]
    end
end
function _𝒞ₐₗₗ(C,u,v)
    # Computes the cdf of the bivariate copula C and its partial derivatives all at once using forwarddiff. 
    _u = ForwardDiff.Dual(u,(1.,0.))
    _v = ForwardDiff.Dual(v,(0.,1.))
    r = cdf(C,[_u,_v])
    return r.value, r.partials[1], r.partials[2]
end
function _mk_∂Λₑ(C, args...)
    return Roots.find_zero(x -> _equation(x, C, args...), 0.0, Roots.Order1())
end
function _equation(∂Λₑ, C, Sₑ₋, Sₚ, ∂Sₚ, ∂Nⱼ)
    # This function expect only indivs with Yⱼ == 1 to be passed
    Sₑ = Sₑ₋ * exp(- ∂Λₑ)
    numer, denom = zero(Sₑ), zero(Sₑ)
    for i in eachindex(Sₚ)
        c, c₁, c₂ = _𝒞ₐₗₗ(C, Sₑ, Sₚ[i])
        numer += ∂Nⱼ[i] / c₁ + ∂Sₚ[i] * c₂ / (c₁ * c)
        denom += Sₑ / c
    end
    return ∂Λₑ - numer/denom
end
function _eq_var(C, Sₑ, Sₚ, ∂Nⱼ)
    ∂Vₑ, Dₑ, r = zero(Sₑ), zero(Sₑ), zero(Sₑ)
    for i in eachindex(Sₚ)
        c, c₁, _ = _𝒞ₐₗₗ(C, Sₑ, Sₚ[i])
        ∂Vₑ += ∂Nⱼ[i] / c₁^2
        r += ∂Nⱼ[i] / c₁
        Dₑ += 1 / c
    end
    return ∂Vₑ, Dₑ*Sₑ, r
end

# Specialisation for the IndependentCopula
function _𝒞ₐₗₗ(::IndependentCopula,u,v)
    return u*v, v, u
end
function _mk_∂Λₑ(::IndependentCopula, Sₑ₋, Sₚ, ∂Sₚ, ∂Nⱼ)
    # In this case we know the solution, no need to rootfind. 
    numer, denom = zero(Sₑ₋), zero(Sₑ₋)
    for i in eachindex(Sₚ)
        wᵢ = 1/Sₚ[i]
        numer += ∂Nⱼ[i] * wᵢ + ∂Sₚ[i] * wᵢ * wᵢ
        denom += wᵢ
    end
    return numer/denom
end
function _eq_var(::IndependentCopula, Sₑ, Sₚ, ∂Nⱼ)
    ∂Vₑ, Dₑ, r = zero(Sₑ), zero(Sₑ), zero(Sₑ)
    for i in eachindex(Sₚ)
        wᵢ = 1/Sₚ[i]
        ∂Vₑ += ∂Nⱼ[i] * wᵢ^2
        r += ∂Nⱼ[i] * wᵢ
        Dₑ += wᵢ
    end
    return ∂Vₑ, Dₑ, r
end

# Specialisation for the ClaytonCopula
function _𝒞ₐₗₗ(C::ClaytonCopula,u,v)
    θ = C.G.θ
    iu = 1/u
    iv = 1/v
    u2 = iu^θ
    v2 = iv^θ
    rr = u2+v2-1
    rr < 0 && return 0.0, 0.0, 0.0 # julia go brrrrrr
    A = 1/rr
    val = A^(1/θ)
    du = A * u2 * val * iu
    dv = A * v2 * val * iv
    return val, du, dv
end
function _equation(∂Λₑ, C::ClaytonCopula, Sₑ₋, Sₚ, ∂Sₚ, ∂Nⱼ)
    # This function expect only indivs with Yⱼ == 1 to be passed
    Sₑ = Sₑ₋ * exp(- ∂Λₑ) # NON-KM-Like updates ? Maybe we coudl change that. 
    numer, denom = zero(Sₑ₋), zero(Sₑ₋)
    θ = C.G.θ
    iθ = 1/θ
    u2 = Sₑ^(-θ)
    B = Sₑ/u2
    for i in eachindex(Sₚ)
        v2 = Sₚ[i]^(-θ)
        rr = u2+v2-1
        if rr > 0
            A = rr^iθ
            numer += A * B * (∂Nⱼ[i] * rr + ∂Sₚ[i] * v2 / Sₚ[i])
            denom += A
        end
    end
    return Sₑ * ∂Λₑ * denom - numer
end
function _eq_var(C::ClaytonCopula, Sₑ, Sₚ, ∂Nⱼ)
    ∂Vₑ, Dₑ, r = zero(Sₑ), zero(Sₑ), zero(Sₑ)
    θ = C.G.θ
    iu = 1/Sₑ
    u2 = iu^θ
    for i in eachindex(Sₚ)
        rr = u2 + (Sₚ[i])^(-θ) - 1
        if rr > 0
            A = 1/rr
            val = A^(1/θ)
            du = A * u2 * val * iu
            ∂Vₑ += ∂Nⱼ[i] / du^2
            r += ∂Nⱼ[i] / du
            Dₑ += 1/val
        end
    end
    return ∂Vₑ, Dₑ*Sₑ, r
end

# Specialisation for the FrankCopula
function _𝒞ₐₗₗ(C::FrankCopula,u,v)
    θ = C.G.θ
    eu = expm1(-u*θ)
    ev = expm1(-v*θ)
    eθ = expm1(-θ)
    euv = eu*ev
    val = -log1p(euv/eθ)/θ
    d = 1/(eθ + euv)
    du = d*(euv+ev)
    dv = d*(euv+eu)
    return val, du, dv
end
function _equation(∂Λₑ, C::FrankCopula, Sₑ₋, Sₚ, ∂Sₚ, ∂Nⱼ)
    # This function expect only indivs with Yⱼ == 1 to be passed
    Sₑ = Sₑ₋ * exp(- ∂Λₑ)
    numer, denom = zero(Sₑ₋), zero(Sₑ₋)
    θ = C.G.θ
    eu = expm1(-Sₑ*θ)
    eθ = expm1(-θ)
    eu_over_eθ = eu/eθ
    for i in eachindex(Sₚ)
        ev = expm1(-Sₚ[i]*θ)
        euv = eu*ev
        ival = -θ/log1p(ev*eu_over_eθ)
        r = 1 / (euv+ev)
        numer += r * (∂Nⱼ[i] * (eθ + euv) + ∂Sₚ[i] * ival * (euv+eu))
        denom += ival
    end
    return Sₑ * ∂Λₑ * denom - numer
end
function _eq_var(C::FrankCopula, Sₑ, Sₚ, ∂Nⱼ)
    ∂Vₑ, Dₑ, r = zero(Sₑ), zero(Sₑ), zero(Sₑ)
    θ = C.G.θ
    eu = expm1(-Sₑ*θ)
    eθ = expm1(-θ)
    for i in eachindex(Sₚ)
        ev = expm1(-Sₚ[i]*θ)
        euv = eu*ev
        du = (euv + ev)/(eθ + euv)
        ∂Vₑ += ∂Nⱼ[i] / du^2
        r += ∂Nⱼ[i] / du
        Dₑ += 1/log1p(euv/eθ)
    end
    return ∂Vₑ, -θ*Dₑ*Sₑ, r
end

# Specialisation for the GumbelCopula
function _𝒞ₐₗₗ(C::GumbelCopula,u,v)
    θ = C.G.θ
    
    lu = -log(u)
    luθ1 = lu^(θ-1)
    luθ = luθ1 * lu

    lv = -log(v)
    lvθ1 = lv^(θ-1)
    lvθ = lvθ1 * lv
    
    r = (luθ+lvθ)
    rγ1 = r^(1/θ-1)
    rγ = rγ1 * r

    val = exp(-rγ)
    du = val * rγ1 * luθ1 / u 
    dv = val * rγ1 * lvθ1 / v 
    return val, du, dv
end