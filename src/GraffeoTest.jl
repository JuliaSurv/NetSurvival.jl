"""
    GraffeoTest

The Grafféo test is a log-rank type test and is typically used in net survival analysis to determine the impact of certain covariates in the study.

The null ``(H_0)`` hypothesis tests the following assumption:

```math
\\forall t \\in [0,T], \\Lambda_{E,g_1}(t) = \\Lambda_{E,g_2}(t) = ... = \\Lambda_{E,g_k}(t)
```

where ``G = {g_1,...,g_k}`` is a partition of ``1,...,n`` consisting of disjoint groups of individuals that we wish to compare to each other. 
For all group ``g \\in G``, let's denote the numerator and denominator of the Pohar Perme (partial) excess hazard estimators, restricted to individuals in the group, by: 

* ``\\partial N_{E,g}(s) = \\sum_{i \\in g} \\frac{\\partial N_i(s)}{S_{P_i}(s)} - \\frac{Y_i(s)}{S_{P_i}(s)}\\partial\\Lambda_{P_i}(s)``
* ``Y_{E,g}(s) = \\sum_{i \\in g} \\frac{Y_i(s)}{S_{P_i}(s)}``
* ``R_{g}(s) = \\frac{Y_{E,g}(s)}{\\sum_{g\\in G} Y_{E,g}(s)}``

Then, define the vector ``\\mathbf Z = \\left(Z_{g_r}: r \\in 1,...,k-1 \\right)`` with entries: 

```math
Z_g(T) = N_{E,g}(s) - \\int_{0}^T Y_{E,g}(s) \\partial\\hat{\\Lambda}_E(s)
```

The test statistic is then given by:

```math
U(T) = \\mathbf Z(T)'\\hat{\\Sigma}_Z^{-1} \\mathbf Z(T)
```

where the entries of the ``\\hat{\\Sigma}_Z`` matrix are given by: 

```math
\\sigma_{g,h}(T) = \\int_0^T \\sum_{\\ell \\in G} \\left(\\delta_{g,\\ell} - R_g(t) \\right)\\left(\\delta_{h,\\ell} - R_h(t)\\right) \\left(\\sum_{i\\in\\ell} \\frac{\\partial N_i(s)}{S^2_{P_i}}\\right)
```

Under ``H_0``, the statistic ``U(T)`` is asymptotically ``\\chi^2(k-1)``-distributed.

To apply the test to your data based on a certain rate table, apply the example below to your code : 

    fit(GraffeoTest, @formula(Surv(time,status)~covariable1 + covariable2), data, ratetable)

If you wish to stratify a covariate:

    fit(GraffeoTest, @formula(Surv(time,status)~covariable1 + Strata(covariable2)), data, ratetable)

The produced test statistics is supposed to follow a chi squared distribution under ``(H_0)``. You can fetch the results using the `.stat`, `.df` and `.pval` fields of the returned object. 

References: 
* [Graffeo2016](@cite) Grafféo, Nathalie and Castell, Fabienne and Belot, Aurélien and Giorgi, Roch (2016). A Log-Rank-Type Test to Compare Net Survival Distributions.  
"""
struct GraffeoTest{Method}
    method::Method
    t::Array{Float64,1}
    ∂N::Array{Float64, 2}
    ∂V::Array{Float64, 2}
    ∂Z::Array{Float64, 2}
    D::Array{Float64, 2}
    R::Array{Float64, 2}
    ∂VZ::Array{Float64, 3}
    stat::Float64
    df::Int64
    pval::Float64
    function GraffeoTest(method::M, T, Δ, age, year, rate_preds, strata, group, ratetable)  where M<:NPNSMethod

        # This version of the test is HIGHLY INNEFICIENT. 
        # We should avoid allocating that much memory. 
        # But for the moment we simply want to check that it is working correctly. 

        grid = mk_grid(T,1) # precision is always 1 ? 

        # get stratas and groups, count them.  
        stratas = unique(strata)
        groups  = unique(group)
        ngroups = length(groups)

        # Allocate: 
        ∂N  = zeros(ngroups, length(grid))
        ∂V  = zeros(ngroups, length(grid))
        ∂Z  = zeros(ngroups, length(grid))
        D   = zeros(ngroups, length(grid))
        R   = zeros(ngroups, length(grid))
        ∂VZ = zeros(ngroups, ngroups, length(grid))

        num_excess   = zero(grid)
        num_pop      = zero(grid)
        num_variance = zero(grid)
        den_pop      = zero(grid)
        den_excess   = zero(grid)
        ∂t = [diff(grid)...,1.0]

        # Compute Pohar Perme numerator and denominators on each strata&group (s,g)
        for s in eachindex(stratas)
            for g in eachindex(groups)
                idx = (group .== groups[g]) .&& (strata .== stratas[s])
                num_excess   .= 0
                num_pop      .= 0
                num_variance .= 0
                den_pop      .= 0
                den_excess   .= 0
                Λ!(method, num_excess, den_excess, num_pop, den_pop, num_variance, T[idx], Δ[idx], age[idx], year[idx], rate_preds[idx,:], ratetable, grid, ∂t)
                ∂N[g, :] .= num_excess .- num_pop
                ∂V[g, :] .= num_variance
                D[g, :]  .= den_excess
            end
            
            R .= ifelse.(sum(D,dims=1) .== 0, 0, D ./ sum(D,dims=1))
            ∂Z .= ∂N .- R .* sum(∂N,dims=1)
        
            # Compute test variance
            for ℓ in eachindex(groups)
                for g in eachindex(groups)
                    for h in eachindex(groups)
                        for t in eachindex(grid)
                            ∂VZ[g, h,t] += ((g==ℓ) - R[g, t]) * ((h==ℓ) - R[h, t]) .* ∂V[ℓ, t]
                        end
                    end
                end
            end
        end

        # Cumulate accross time
        Z =  dropdims(sum(∂Z, dims=2), dims=2)
        VZ = dropdims(sum(∂VZ, dims=3), dims=3)

        # Finally compute the stat and p-values:
        stat = dot(Z[1:end-1],(VZ[1:end-1,1:end-1] \ Z[1:end-1])) # test statistic
        df = ngroups-1 # number of degree of freedom of the chi-square test
        pval = ccdf(Chisq(df), stat[1]) # Obtained p-value. 
        return new{M}(method,grid,∂N, ∂V, ∂Z, D, R, ∂VZ, stat[1], df, pval)
    end
end

struct GraffeoTestHolder{T}
    m::T
end

GraffeoTest(T, Δ, age, year, rate_preds, strata, group, ratetable) = GraffeoTest(PoharPermeMethod(),T, Δ, age, year, rate_preds, strata, group, ratetable)

GraffeoTest()                                    = GraffeoTestHolder(PoharPermeMethod())
GraffeoTest(m::M)      where M<:NPNSMethod       = GraffeoTestHolder(m)
GraffeoTest{M}()       where M<:NPNSMethod       = GraffeoTestHolder(M())
GraffeoTest(::Type{M}) where M<:NPNSMethod       = GraffeoTestHolder(M())
GraffeoTest(::Type{NPNSEstimator{M}}) where M<:NPNSMethod       = GraffeoTestHolder(M())
GraffeoTest(C::Cop)    where Cop<:Copulas.Copula = GraffeoTest(GenPoharPermeMethod(C))

StatsBase.fit(X::GraffeoTestHolder{M}, formula, df, rt) where M<:NPNSMethod   = StatsBase.fit(GraffeoTest, X.m,                formula, df, rt)
StatsBase.fit(::Type{GraffeoTest},     formula, df, rt)                       = StatsBase.fit(GraffeoTest(PoharPermeMethod()), formula, df, rt)
StatsBase.fit(::Type{GraffeoTest{M}},  formula, df, rt) where {M<:NPNSMethod} = StatsBase.fit(GraffeoTest, M(),                formula, df, rt)

function StatsBase.fit(GTH::GraffeoTestHolder{M}, formula::FormulaTerm, df::DataFrame, rt::RateTables.AbstractRateTable) where {M<:NPNSMethod}
    terms = StatsModels.termvars(formula.rhs)
    tf = typeof(formula.rhs)
    types = (tf <: AbstractTerm) ? [tf] : typeof.(formula.rhs)
    are_strata = [t <: FunctionTerm{typeof(Strata)} for t in types]

    strata = groupindices(groupby(df,terms[are_strata]))
    group  = groupindices(groupby(df,terms))

    resp = modelcols(apply_schema(formula,schema(df)).lhs,df)
    rate_predictors = _get_rate_predictors(rt,df)
    
    return GraffeoTest(GTH.m, resp[:,1], resp[:,2], df.age, df.year, select(df,rate_predictors), strata, group, rt)
end

# A show function: 
function Base.show(io::IO, test::GraffeoTest)
    println(io, "Grafféo's log-rank-type-test (Method: $(test.method))")
    df = DataFrame(test_statistic = test.stat, degrees_of_freedom = test.df, p_value = test.pval)
    show(io, df)
end

# Potential other extraction methods: 
function statistic(X::GraffeoTest, at_time_T)
    # Cumulate accross time
    i_T = findlast(X.t .<= at_time_T)
    Z =  dropdims(sum(X.∂Z[:,1:i_T], dims=2), dims=2)
    Γ = dropdims(sum(X.∂Γ[:,:,1:i_T], dims=3), dims=3)
    stat = dot(Z[1:end-1],Γ[1:end-1,1:end-1] \ Z[1:end-1]) # test statistic
    return stat
end 
function pvalue(X::GraffeoTest, at_time_T)
    stat = statistic(X, at_time_T)
    df = size(X.∂Z,1)-1 # number of degree of freedom of the chi-square test
    return ccdf(Chisq(df), stat[1]) # Obtained p-value.
end  
