"""
    GraffeoTest

The Grafféo test is a log-rank type test and is typically used in net survival analysis to determine the impact of certain covariates in the study.

The null hypothesis ``H_0`` tests the following assumption:

```math
\\forall t \\in [0,T], \\Lambda_{E,1}(t) = \\Lambda_{E,2}(t) = ... = \\Lambda_{E,k}(t)$$
```

For this test, we first define the number of deaths caused by the event studied for a time ``s`` within the group ``h`` noted ``N_{E,h}^w(s)`` and the process of individuals at risk within the same group $h$ at time $s$ noted $Y_h^w(s)$. Both of these values are weighted with the populational estimated survival for the given patient, same as in Pohar Perme. 

Thus, we define ``\\forall h \\in [1;k]``:

```math
Z_h^w(T) = \\int_0^T \\mathbf{1}(Y_.^w(s) > 0)dN_{E,h}^w(s) - \\int_0^T \\mathbf{1}(Y_.^w(s) > 0)\\frac{Y_h^w(s)}{Y_.^w(s)}dN_{E,.}^w(s)
```

With ``Y_.^w(s) = \\sum_{h=i}^k Y_h^w(s)`` and ``dN_{E,.}^w(s) = \\sum_{h=i}^k dN_{E,h}^w(s)``.

The test statistic, defined under ``H_0``, is then given by:

```math
U^w(T) = (Z_1^w(T), ..., Z_{k-1}^w(T)) \\hat{\\Sigma}_0^w(T)^{-1} (Z_1^w(T), ..., Z_{k-1}^w(T))^t
```

Knowing that ``\\hat{\\Sigma}_0^w(T) = (\\hat{\\sigma}^w_{hj}(T))_{hj}`` for ``(h,j) \\in [1;k-1]^2``

Finally, the asymptotic distribution of the test statistic under ``H_0``:

```math
U^w(t) \\sim \\chi^2(k-1)
```

To apply the test to your data based on a certain rate table, apply the example below to your code : 

    fit(GraffeoTest, @formula(Surv(time,status)~covariable1 + covariable2), data, ratetable)

If you wish to stratify a covariate:

    fit(GraffeoTest, @formula(Surv(time,status)~covariable1 + Strata(covariable2)), data, ratetable)

References: 
* [GraffeoTest](@cite) Grafféo, Nathalie and Castell, Fabienne and Belot, Aurélien and Giorgi, Roch (2016). A Log-Rank-Type Test to Compare Net Survival Distributions.  
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
        
        # Compute test variance on each strata&group (s,g)
        for s in eachindex(stratas)
            for g in eachindex(groups)
                for h in eachindex(groups)
                    for ℓ in eachindex(groups)
                        ∂VZ[s, g, h,:] .+= ((g==ℓ) .- R[s, g, :]) .* ((h==ℓ) .- R[s, h, :]) .* ∂V[s, ℓ, :]
                    end
                end
            end
        end

        # Cumulate accross time and stratas: (but not last time value)
        Z =  dropdims(sum(∂Z[:,:,1:end-1], dims=(1,3)), dims=(1,3))
        VZ = dropdims(sum(∂VZ[:,:,:,1:end-1], dims=(1,4)), dims=(1,4))

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
