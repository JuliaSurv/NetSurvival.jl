struct Nessie
    expected_sample_size::Vector{Float64}
    expected_life_time::Float64
    grid::Vector{Float64}
    function Nessie(T, Δ, age, year, rate_preds, ratetable)
        annual_grid = 0:RateTables.RT_DAYS_IN_YEAR:maximum(T)
        exp_spl_size = zeros(length(annual_grid))
        exp_life_time = 0.0
        for i in eachindex(age)
            Pᵢ = Life(ratetable[rate_preds[i,:]...], age[i], year[i])
            for j in eachindex(annual_grid)
                exp_spl_size[j] += ccdf(Pᵢ, annual_grid[j])
            end
            exp_life_time += expectation(Pᵢ)
        end
        return new(exp_spl_size, exp_life_time / RateTables.RT_DAYS_IN_YEAR / length(age), annual_grid)
    end
end

"""
    nessie 

To call this function, use the formula below: 

    nessie(@formula(Surv(time,status)~covariate), data, ratetable)
"""
function nessie(args...)
    r = fit(Nessie,args...)
    if (typeof(r)<:Nessie)
        return r
    end
    transform!(r, :estimator => ByRow(x-> (x.grid, x.expected_life_time, x.expected_sample_size)) => [:grid, :expected_life_time,:expected_sample_size])
    select!(r, Not(:estimator))

    lt = deepcopy(r)
    select!(lt, Not([:expected_sample_size, :grid]))

    select!(r, Not(:expected_life_time))
    return lt, r
end

# Maybe not necessary ? No need to clutter the interface too much.. 
expected_life_time(x::Nessie) = x.expected_life_time
expected_sample_size(x::Nessie) = x.expected_sample_size

function StatsBase.fit(::Type{Nessie}, formula::FormulaTerm, df::DataFrame, rt::RateTables.AbstractRateTable)
    rate_predictors = _get_rate_predictors(rt,df)
    formula_applied = apply_schema(formula,schema(df))

    if isa(formula.rhs, ConstantTerm) # No predictors
        resp = modelcols(formula_applied.lhs, df)
        return Nessie(resp[:,1], resp[:,2], df.age, df.year, select(df,rate_predictors), rt)
    else
        gdf = groupby(df, StatsModels.termnames(formula.rhs))
        return rename(
            combine(
                gdf,
                dfᵢ -> begin
                    resp2 = modelcols(formula_applied.lhs, dfᵢ)
                    Nessie(resp2[:,1], resp2[:,2], dfᵢ.age, dfᵢ.year, select(dfᵢ, rate_predictors), rt)
                end
            ),
            :x1 => :estimator
        )
    end
end
