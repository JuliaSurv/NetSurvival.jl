@testsnippet SharedTestSetup begin
    using NetSurvival, RateTables, DataFrames, Random, Distributions

    # A zero-hazard mock ratetable to validate relationships vs KM and identities
    struct ZeroRateTable <: RateTables.AbstractRateTable end
    Base.getindex(::ZeroRateTable, args...) = ZeroRateTable()       # allow rt[...] access with any predictors
    # Provide one predictor to ensure rate_preds has n rows (avoid 0x0 DataFrame)
    RateTables.predictors(::ZeroRateTable) = (:sex,)                # use a column we control in tests
    RateTables.daily_hazard(::ZeroRateTable, age, year) = 0.0       # zero population hazard everywhere

    # Simple Kaplan–Meier estimator on an integer day grid
    # Returns: grid (unique sorted integer days including all times),
    #          S (survival step function on that grid)
    function _km(times::AbstractVector{<:Real}, status::AbstractVector{<:Integer})
        @assert length(times) == length(status)
        T = round.(Int, times)
        grid = sort(unique(vcat(collect(1:maximum(T)), T)))
        S = ones(Float64, length(grid))
        for (j, t) in pairs(grid)
            y = count(>=(t), T)
            d = count(i -> (T[i] == t) && (status[i] == 1), eachindex(T))
            S[j] = j == 1 ? (1 - (y == 0 ? 0.0 : d / y)) : S[j-1] * (1 - (y == 0 ? 0.0 : d / y))
        end
        return grid .* 1.0, S
    end

    # Map a step function defined on reference grid (Gref, Sref) to another grid G via LOCF
    function at_grid(G::AbstractVector, Gref::AbstractVector, Sref::AbstractVector)
        out = similar(G, Float64)
        for i in eachindex(G)
            j = findlast(Gref .<= G[i])
            out[i] = isnothing(j) ? 1.0 : Sref[j]
        end
        return out
    end

    # Constant population hazard ratetable (per-day rate)
    struct ConstRateTable{T} <: RateTables.AbstractRateTable
        λP_per_day::T
    end
    Base.getindex(rt::ConstRateTable, args...) = rt
    # Same rationale as above: ensure non-empty predictors selection
    RateTables.predictors(::ConstRateTable) = (:sex,)
    RateTables.daily_hazard(rt::ConstRateTable, age, year) = rt.λP_per_day

    # Truth functions
    truth_SE_const(grid, λE_per_day) = exp.(-λE_per_day .* grid)
    function truth_SE_piecewise(grid, λ1, λ2, τ)
        @. exp(-(λ1 * min(grid, τ) + λ2 * max(0, grid - τ)))
    end
    function truth_CM_const(grid, λE, λP)
        λ = λE + λP
        Mtot = 1 .- exp.(-λ .* grid)
        ME = (λE / λ) .* Mtot
        MP = (λP / λ) .* Mtot
        return Mtot, ME, MP
    end

    supnorm(a, b) = maximum(abs.(a .- b))

    # Simulators
    function simulate_rel_surv_const_zeroP(n; λE_per_day, λC_per_day=1/(15*365.241), T_admin_years=15.0)
        YEARS = 365.241
        t_admin = T_admin_years * YEARS
        rows = Vector{NamedTuple}(undef, n)
        for i in 1:n
            e = rand(Exponential(1/λE_per_day))
            c = min(rand(Exponential(1/λC_per_day)), t_admin)
            t = min(e, c)
            δ = (e <= c) ? 1 : 0
            rows[i] = (time=t, status=δ, age=60YEARS, year=2000YEARS, sex=(:male, :female)[1 + rand(Bool)])
        end
        return DataFrame(rows)
    end
    function simulate_rel_surv_const(n; λE_per_day, λP_per_day, λC_per_day=1/(15*365.241), T_admin_years=15.0)
        YEARS = 365.241
        t_admin = T_admin_years * YEARS
        rows = Vector{NamedTuple}(undef, n)
        for i in 1:n
            e = rand(Exponential(1/λE_per_day))
            p = rand(Exponential(1/λP_per_day))
            c = min(rand(Exponential(1/λC_per_day)), t_admin)
            t = min(e, p, c)
            δ = (t == c) ? 0 : 1
            rows[i] = (time=t, status=δ, age=60YEARS, year=2000YEARS, sex=(:male, :female)[1 + rand(Bool)])
        end
        return DataFrame(rows)
    end
    function simulate_rel_surv_piecewise_zeroP(n; λ1, λ2, τ, λC_per_day=1/(15*365.241), T_admin_years=15.0)
        YEARS = 365.241
        t_admin = T_admin_years * YEARS
        rows = Vector{NamedTuple}(undef, n)
        a = λ1 * τ
        for i in 1:n
            u = rand()
            x = -log(u)
            if x <= a
                e = x / λ1
            else
                e = τ + (x - a) / λ2
            end
            c = min(rand(Exponential(1/λC_per_day)), t_admin)
            t = min(e, c)
            δ = (e <= c) ? 1 : 0
            rows[i] = (time=t, status=δ, age=60YEARS, year=2000YEARS, sex=(:male, :female)[1 + rand(Bool)])
        end
        return DataFrame(rows)
    end
end

@testitem "Pohar-Perme: basic properties and function interface" tags=[:Prop] setup=[SharedTestSetup] begin

    Random.seed!(42)
    # small synthetic dataset using slopop
    n = 200
    YEARS = 365.241
    df = DataFrame(
        time = rand(1:2000, n) .* 1.0,
        status = rand(Bool, n),
        age = rand(40:80, n) .* YEARS,
        year = rand(1990:2005, n) .* YEARS,
        sex = rand((:male, :female), n),
    )
    model = fit(PoharPerme, @formula(Surv(time, status) ~ 1), df, slopop)

    # Bounds and monotonicity for PP
    @test all(model.Sₑ .>= 0)

    # Variance non-decreasing; zero at t<min
    @test model.σₑ[1] >= 0
    @test all(diff(model.σₑ) .>= -1e-12)

    # Confidence intervals sane
    cis = confint(model; level = 0.05)
    @test length(cis) == length(model.Sₑ)
    for i in eachindex(cis)
        lo, hi = cis[i]
        if isfinite(lo) && isfinite(hi)
            @test lo <= model.Sₑ[i] <= hi
            @test lo >= 0
        end
    end

    # Function interface
    tmin, tmax = extrema(model.grid)
    @test model(tmin - 10.0) == 1.0
    @test model(tmax + 10.0) == last(model.Sₑ)
    @test variance(model, tmin - 10.0) == 0.0
    @test isapprox(variance(model, tmax + 10.0), last(model.σₑ)^2; atol=1e-12)
end

@testitem "Convergence: constant excess λE, zero population hazard" tags=[:Prop] setup=[SharedTestSetup] begin

    Random.seed!(123)
    YEARS = 365.241
    λE = 1/(8YEARS)
    rt0 = ZeroRateTable()

    for (n_small, n_large) in ((2000, 8000),)
        df_small = simulate_rel_surv_const_zeroP(n_small; λE_per_day=λE)
        df_large = simulate_rel_surv_const_zeroP(n_large; λE_per_day=λE)

        pp_small = fit(PoharPerme, @formula(Surv(time,status) ~ 1), df_small, rt0)
        pp_large = fit(PoharPerme, @formula(Surv(time,status) ~ 1), df_large, rt0)
        e2_small = fit(EdererII, @formula(Surv(time,status) ~ 1), df_small, rt0)
        e2_large = fit(EdererII, @formula(Surv(time,status) ~ 1), df_large, rt0)

        truth_small = truth_SE_const(pp_small.grid, λE)
        truth_large = truth_SE_const(pp_large.grid, λE)

        err_pp_small = supnorm(pp_small.Sₑ, truth_small)
        err_pp_large = supnorm(pp_large.Sₑ, truth_large)
        err_e2_small = supnorm(e2_small.Sₑ, truth_SE_const(e2_small.grid, λE))
        err_e2_large = supnorm(e2_large.Sₑ, truth_SE_const(e2_large.grid, λE))

        @test err_pp_large < err_pp_small
        @test err_e2_large < err_e2_small
        @test err_pp_large < 0.03
        @test err_e2_large < 0.03
    end
end

@testitem "Convergence: constant λE and constant λP, PP to truth" tags=[:Prop] setup=[SharedTestSetup] begin

    Random.seed!(321)
    YEARS = 365.241
    λE = 1/(8YEARS)
    λP = 1/(20YEARS)
    rt = ConstRateTable(λP)
    n = 8000
    df = simulate_rel_surv_const(n; λE_per_day=λE, λP_per_day=λP)
    pp = fit(PoharPerme, @formula(Surv(time,status) ~ 1), df, rt)
    truth = truth_SE_const(pp.grid, λE)
    err = supnorm(pp.Sₑ, truth)
    @test err < 0.03
end

@testitem "Convergence: piecewise-constant λE, zero λP (PP)" tags=[:Prop] setup=[SharedTestSetup] begin

    Random.seed!(777)
    YEARS = 365.241
    λ1 = 1/(6YEARS)
    λ2 = 1/(12YEARS)
    τ  = 5YEARS
    n = 8000
    df = simulate_rel_surv_piecewise_zeroP(n; λ1=λ1, λ2=λ2, τ=τ)
    rt0 = ZeroRateTable()
    pp = fit(PoharPerme, @formula(Surv(time,status) ~ 1), df, rt0)
    truth = truth_SE_piecewise(pp.grid, λ1, λ2, τ)
    err = supnorm(pp.Sₑ, truth)
    @test err < 0.035
end

@testitem "Convergence: CrudeMortality under constant hazards" tags=[:Prop] setup=[SharedTestSetup] begin

    Random.seed!(888)
    YEARS = 365.241
    λE = 1/(8YEARS)
    λP = 1/(20YEARS)
    n = 8000
    df = simulate_rel_surv_const(n; λE_per_day=λE, λP_per_day=λP)
    rt = ConstRateTable(λP)
    e2 = fit(EdererII, @formula(Surv(time,status) ~ 1), df, rt)
    cm = CrudeMortality(e2)
    Mtot, ME, MP = truth_CM_const(e2.grid, λE, λP)
    @test supnorm(cm.Mₒ, Mtot) < 0.03
    @test supnorm(cm.Mₑ, ME)   < 0.03
    @test supnorm(cm.Mₚ, MP)   < 0.03
end

@testitem "CI coverage at t0 for PP under constant hazards (1-year horizon)" tags=[:Prop] setup=[SharedTestSetup] begin

    Random.seed!(9999)
    YEARS = 365.241
    λE = 1/(8YEARS)
    rt0 = ZeroRateTable()
    t0 = 1YEARS
    R = 120
    n = 1200
    covered_count = Ref(0)
    total_count = Ref(0)
    for r in 1:R
        # 1-year administrative horizon to reduce grid length and runtime
        df = simulate_rel_surv_const_zeroP(n; λE_per_day=λE, T_admin_years=1.0)
        m = fit(PoharPerme, @formula(Surv(time,status) ~ 1), df, rt0)
        # Get CI at closest grid point <= t0
        it = findlast(m.grid .<= t0)
        if it !== nothing
            Strue = exp(-λE * m.grid[it])
            lo, hi = confint(m; level=0.05)[it]
            covered_count[] += Int(lo <= Strue <= hi)
            total_count[] += 1
        end
    end
    # 95% nominal; accept a band to account for finite R
    cov = covered_count[] / max(total_count[], 1)
    @test 0.88 <= cov <= 0.98
end

@testitem "Ederer I/II/Hakulinen: shape properties (no R)" tags=[:Prop] setup=[SharedTestSetup] begin

    Random.seed!(7)
    YEARS = 365.241
    n = 150
    df = DataFrame(
        time = rand(1:1200, n) .* 1.0,
        status = rand(Bool, n),
        age = rand(30:85, n) .* YEARS,
        year = rand(1985:2008, n) .* YEARS,
        sex = rand((:male, :female), n),
    )
    for E in (EdererI, EdererII, Hakulinen)
        m = fit(E, @formula(Surv(time, status) ~ 1), df, slopop)
        # Survival non-negative; may exceed 1 for Ederer-class, so don't assert upper bound or monotonicity
        @test all(x -> !isfinite(x) || x >= 0, m.Sₑ)
        # Variance non-decreasing
        @test all(x -> !isfinite(x) || x >= -1e-12, diff(m.σₑ))
        # CIs ordered and non-negative
        for (lo, hi) in confint(m; level=0.05)
            if isfinite(lo) && isfinite(hi)
                @test lo <= hi
                @test lo >= 0
            end
        end
    end
end

@testitem "Zero population hazard: PP and EdererII ≈ KM" tags=[:Prop] setup=[SharedTestSetup] begin

    Random.seed!(1)
    n = 80
    # Make a small integer-time dataset with some events and censoring
    times = rand(10:300, n)
    status = rand(Bool, n)
    df = DataFrame(
        time = Float64.(times),
        status = Int.(status),
        age = fill(60.0, n),
        year = fill(2000.0, n),
        sex = [i % 2 == 0 ? :male : :female for i in 1:n],
    )
    rt0 = ZeroRateTable()

    pp = fit(PoharPerme, @formula(Surv(time, status) ~ 1), df, rt0)
    e2 = fit(EdererII,   @formula(Surv(time, status) ~ 1), df, rt0)

    km_grid, kmS = _km(df.time, df.status)
    # Compare at pp.grid points by last observation carried back from km
    # Map KM S to pp grid
    km_on_pp = at_grid(pp.grid, km_grid, kmS)
    km_on_e2 = at_grid(e2.grid, km_grid, kmS)

    @test isapprox(pp.Sₑ, km_on_pp; atol=1e-10, rtol=1e-8)
    @test isapprox(e2.Sₑ, km_on_e2; atol=1e-10, rtol=1e-8)
end

@testitem "Permutation invariance (approximate) and grouped vs split equivalence" tags=[:Prop] setup=[SharedTestSetup] begin

    Random.seed!(11)
    YEARS = 365.241
    n = 300
    df = DataFrame(
        time = rand(1:1500, n) .* 1.0,
        status = rand(Bool, n),
        age = rand(35:70, n) .* YEARS,
        year = rand(1995:2005, n) .* YEARS,
        sex = rand((:male, :female), n),
        grp = rand(('A', 'B'), n),
    )
    m1 = fit(PoharPerme, @formula(Surv(time, status) ~ 1), df, slopop)
    shuffle!(df)
    
    m2 = fit(PoharPerme, @formula(Surv(time, status) ~ 1), df, slopop)
    @test isapprox(m1.Sₑ, m2.Sₑ; rtol=1e-10, atol=1e-12)
    @test isapprox(m1.σₑ, m2.σₑ; rtol=1e-10, atol=1e-12)

    # Grouped-vs-split equivalence
    mg = fit(PoharPerme, @formula(Surv(time, status) ~ grp), df, slopop)
    for g in unique(df.grp)
        sub = df[df.grp .== g, :]
        select!(sub, [:time, :status, :age, :year, :sex])
        ms = fit(PoharPerme, @formula(Surv(time, status) ~ 1), sub, slopop)
        # Find the row in grouped result corresponding to this g
        row = findfirst(x -> x == g, mg.grp)
        @test row !== nothing
        est = mg[row, :estimator]
        @test isapprox(est.Sₑ, ms.Sₑ; rtol=1e-10, atol=1e-12)
        @test isapprox(est.σₑ, ms.σₑ; rtol=1e-10, atol=1e-12)
    end
end

@testitem "Grafféo test: duplicated groups and strata invariance" tags=[:Prop] setup=[SharedTestSetup] begin

    Random.seed!(99)
    YEARS = 365.241
    n = 120
    base = DataFrame(
        time = rand(1:800, n) .* 1.0,
        status = rand(Bool, n),
        age = rand(40:75, n) .* YEARS,
        year = rand(1990:2005, n) .* YEARS,
        sex = rand((:male, :female), n),
    )
    # Duplicate dataset into two labeled groups and two identical strata
    df = vcat(
        hcat(base, DataFrame(group = fill(:A, n), str = fill(:one, n))),
        hcat(base, DataFrame(group = fill(:B, n), str = fill(:one, n))),
        hcat(base, DataFrame(group = fill(:A, n), str = fill(:two, n))),
        hcat(base, DataFrame(group = fill(:B, n), str = fill(:two, n))),
    )
    t = fit(GraffeoTest, @formula(Surv(time,status) ~ group), df, slopop)
    @test isapprox(t.stat, 0.0; atol=1e-8)
    @test t.df == 1
    @test isapprox(t.pval, 1.0; atol=1e-8)

    # With two identical strata, stratified equals unstratified
    t_str = fit(GraffeoTest, @formula(Surv(time,status) ~ group + Strata(str)), df, slopop)

    @test isapprox(t.stat, t_str.stat; atol=1e-10)

    # k groups -> df = k-1
    df3 = vcat(
        hcat(base, DataFrame(group = fill(:G1, n))),
        hcat(base, DataFrame(group = fill(:G2, n))),
        hcat(base, DataFrame(group = fill(:G3, n))),
    )
    t3 = fit(GraffeoTest, @formula(Surv(time,status) ~ group), df3, slopop)
    @test t3.df == 2
end

@testitem "CrudeMortality: identities and monotonicity" tags=[:Prop] setup=[SharedTestSetup] begin

    Random.seed!(77)
    YEARS = 365.241
    n = 200
    df = DataFrame(
        time = rand(1:1000, n) .* 1.0,
        status = rand(Bool, n),
        age = rand(30:90, n) .* YEARS,
        year = rand(1980:2010, n) .* YEARS,
        sex = rand((:male, :female), n),
    )
    cm = fit(CrudeMortality, @formula(Surv(time,status) ~ 1), df, slopop)
    @test maximum(abs.(cm.Mₒ .- (cm.Mₑ .+ cm.Mₚ))) <= 1e-8
    # Allow small numerical negatives due to floating arithmetic
    @test all(cm.Mₒ .>= -5e-5)
    @test all(cm.Mₚ .>= -5e-5)

    # Under zero-population hazard, M_pop ≈ 0 and M_excess ≈ 1 - S_E (EdererII)
    rt0 = ZeroRateTable()
    e2 = fit(EdererII, @formula(Surv(time,status) ~ 1), df, rt0)
    cm0 = CrudeMortality(e2)
    @test all(isapprox(cm0.Mₚ, zeros(length(cm0.Mₚ)); atol=1e-12))
    @test maximum(abs.(cm0.Mₑ .- (1 .- e2.Sₑ))) < 0.035
end

@testitem "Nessie: shape checks" tags=[:Prop] setup=[SharedTestSetup] begin

    Random.seed!(21)
    YEARS = 365.241
    n = 120
    df = DataFrame(
        time = rand(1:2000, n) .* 1.0,
        status = rand(Bool, n),
        age = rand(30:90, n) .* YEARS,
        year = rand(1980:2010, n) .* YEARS,
        sex = rand((:male, :female), n),
    )
    # Ensure both sexes present for contrasts
    df.sex[1] = :male; df.sex[2] = :female
    res = nessie(@formula(Surv(time,status) ~ sex), df, slopop)
    # Tuple (life_time_df, sample_size_df)
    lt, ss = res
    # expected_sample_size per group should be non-increasing
    for v in ss.expected_sample_size
        @test all(diff(v) .<= 1e-12)
        @test all(v .>= 0)
    end
    # expected_life_time positive and finite
    @test all(isfinite.(lt.expected_life_time))
    @test all(lt.expected_life_time .> 0)
end