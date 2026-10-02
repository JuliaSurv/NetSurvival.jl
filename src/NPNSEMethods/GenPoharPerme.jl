struct GenPoharPermeMethod{CopulaType} <:NPNSMethod 
    C::CopulaType
    function GenPoharPermeMethod(C::Copulas.Copula{2})
        if isa(C, IndependentCopula)
            return PoharPermeMethod()
        else
            return new{typeof(C)}(C)
        end
    end
end
GenPoharPermeMethod() = GenPoharPermeMethod(IndependentCopula(2))

"""
    GenPoharPerme

This method estimates net survival probabilities under the generalized pohar perme method, taking into account a copula for (E,P) as described in the reference paper. You may use it as: 

    fit(GenPoharPerme(copula), @formula(Surv(time, status)~ x1 + x2), data, ratetable)

Note that `GenPoharPerme(::IndependenceCopula)` simply returns `PoharPerme`.

References: 
* [Laverny2025](@cite) Laverny, O., Grafféo, N., & Giorgi, R. (2025). Non-parametric estimation of net survival under dependence between death causes. arXiv preprint arXiv:2502.09273.
"""
const GenPoharPerme = NPNSEstimator{GenPoharPermeMethod}
GenPoharPerme(C::Copulas.Copula{2}) = GenPoharPermeMethod(C::Copulas.Copula{2})

function _𝒞ₐₗₗ(C,u,v)
    # Computes the cdf of the bivariate copula C and its partial derivatives all at once using forwarddiff. 
    _u = ForwardDiff.Dual(u,(1.,0.))
    _v = ForwardDiff.Dual(v,(0.,1.))
    r = cdf(C,[_u,_v])
    return r.value, r.partials[1], r.partials[2]
end

function _𝒞ₐₗₗ(::IndependentCopula,u,v)
    return u*v, v, u
end

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
