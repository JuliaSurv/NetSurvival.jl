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

References: 
* [Graffeo2016](@cite) Grafféo, Nathalie and Castell, Fabienne and Belot, Aurélien and Giorgi, Roch (2016). A Log-Rank-Type Test to Compare Net Survival Distributions.  
"""
struct GraffeoTest
    ∂N::Array{Float64, 3}
    ∂V::Array{Float64, 3}
    ∂Z::Array{Float64, 3}
    D::Array{Float64, 3}
    R::Array{Float64, 3}
    ∂VZ::Array{Float64, 4}
    stat::Float64
    df::Int64
    pval::Float64
    function GraffeoTest(T, Δ, age, year, rate_preds, strata, group, ratetable)

        # This version of the test is HIGHLY INNEFICIENT. 
        # We should avoid allocating that much memory. 
        # But for the moment we simply want to check that it is working correctly. 

        grid = mk_grid(T,1) # precision is always 1 ? 

        # get stratas and groups, count them.  
        stratas = unique(strata)
        groups  = unique(group)
        nstrata = length(stratas)
        ngroups = length(groups)

        # Allocate: 
        ∂N  = zeros(nstrata, ngroups, length(grid))
        ∂V  = zeros(nstrata, ngroups, length(grid))
        ∂Z  = zeros(nstrata, ngroups, length(grid))
        D   = zeros(nstrata, ngroups, length(grid))
        R   = zeros(nstrata, ngroups, length(grid))
        ∂VZ = zeros(nstrata, ngroups, ngroups, length(grid))

        # Compute Pohar Perme numerator and denominators on each strata&group (s,g)
        for s in eachindex(stratas)
            for g in eachindex(groups)
                idx = (group .== groups[g]) .&& (strata .== stratas[s])
                ∂N[s, g, :], ∂V[s, g, :], D[s, g, :] = _Λ(T[idx], Δ[idx], age[idx], year[idx], rate_preds[idx,:], ratetable, grid)
            end
        end

        # renormalize on groups, be carefull for zeros. 
        R .= ifelse.(sum(D,dims=2) .== 0, 0, D ./ sum(D,dims=2))
        ∂Z .= ∂N .- R .* sum(∂N,dims=2)
        
        # Compute test variance on each strata
        for s in eachindex(stratas)
            for ℓ in eachindex(groups)
                for g in eachindex(groups)
                    for h in eachindex(groups)
                        for t in eachindex(grid)
                            ∂VZ[s, g, h,t] += ((g==ℓ) - R[s, g, t]) * ((h==ℓ) - R[s, h, t]) .* ∂V[s, ℓ, t]
                        end
                    end
                end
            end
        end

        # Cumulate accross time and stratas
        Z =  dropdims(sum(∂Z, dims=(1,3)), dims=(1,3))
        VZ = dropdims(sum(∂VZ, dims=(1,4)), dims=(1,4))

        # Finally compute the stat and p-values:
        stat = dot(Z[1:end-1],(VZ[1:end-1,1:end-1] \ Z[1:end-1])) # test statistic
        df = ngroups-1 # number of degree of freedom of the chi-square test
        pval = ccdf(Chisq(df), stat[1]) # Obtained p-value. 
        return new(∂N, ∂V, ∂Z, D, R, ∂VZ, stat[1], df, pval)
    end
end

function Base.show(io::IO, test::GraffeoTest)
    println(io, "Grafféo's log-rank-type-test")
    df = DataFrame(test_statistic = test.stat, degrees_of_freedom = test.df, p_value = test.pval)
    show(io, df)
end

# The fitting and formula interfaces should be here. 

function StatsBase.fit(::Type{E}, formula::FormulaTerm, df::DataFrame, rt::RateTables.AbstractRateTable) where {E<:GraffeoTest}
    rate_predictors = String.([RateTables.predictors(rt)...])

    expected_columns = [rate_predictors...,"age","year"]
    missing_columns = filter(name -> !(name in names(df)), expected_columns)
    if !isempty(missing_columns)
        throw(ArgumentError("Missing columns in data: $missing_columns"))
    end

    strata = ones(nrow(df))
    group = ones(nrow(df))
    strata_terms = []
    group_terms = []

    if typeof(formula.rhs) == Term
        group = select(df, StatsModels.termvars(formula.rhs))
        group = [join(row, " ") for row in eachrow(group)]
    elseif typeof(formula.rhs) <: FunctionTerm{typeof(Strata)}
        strata = select(df, StatsModels.termvars(formula.rhs))
        strata = [join(row, " ") for row in eachrow(strata)]
    else
        for myterm in formula.rhs
            is_strata = typeof(myterm) <: FunctionTerm{typeof(Strata)}
            if is_strata
                append!(strata_terms, StatsModels.termvars(myterm))
            else
                push!(group_terms, Symbol(myterm))
            end
        end
    end

    if !isempty(group_terms)
        group = select(df, group_terms)
        group = [join(row, " ") for row in eachrow(group)]
    end

    if !isempty(strata_terms)
        strata = select(df, strata_terms)
        strata = [join(row, " ") for row in eachrow(strata)]
    end

    formula = apply_schema(formula,schema(df))
    resp = modelcols(formula.lhs,df)

    return GraffeoTest(resp[:,1], resp[:,2], df.age, df.year, select(df,rate_predictors), strata, group, rt)
end
