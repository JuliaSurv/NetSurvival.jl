@testitem "Assess all NPNSEstimators" begin
    using RateTables, NetSurvival, RCall, DataFrames
    function test_surv(r_method,::Type{E},df, rt, args...) where E

        # Main instance: 
        instance = fit(E, @formula(Surv(time,status)~1), df, rt)
        ci = confint(instance; level = 0.05)

        # Check for no-nans :
        @test !any(isnan.(instance.Sₑ))
        @test !any(isnan.(instance.∂Λₑ))
        @test !any(isnan.(instance.σₑ))
        @test !any([any(isnan.(x)) for x in ci])

        # Check for correspondance of the alternative call syntax: 
        v2 = E(args..., rt)
        @test all(v2.Sₑ .== instance.Sₑ)
        @test all(v2.∂Λₑ .== instance.∂Λₑ)
        @test all(v2.σₑ .== instance.σₑ)

        # Compare with R: 
        @rput r_method
        R"""
        rez = relsurv::rs.surv(survival::Surv(time, stat) ~ 1, rmap=list(age = age, sex = sex, year = diag), data = relsurv::colrec, ratetable = relsurv::slopop, method = r_method, add.times=1:8149)
        t = rez$time
        s = rez$surv
        e = rez$std.err
        """
        r = @rget rez

        err_S = zeros(length(r[:time]))
        err_σ = zeros(length(r[:time]))
        for i in eachindex(r[:time]) 
            j = searchsortedlast(instance.grid, r[:time][i])
            err_S[i] = (r[:surv][i] - instance.Sₑ[j]) / r[:surv][i]
            err_σ[i] = (r[:std_err][i] - instance.σₑ[j]) / r[:std_err][i]
        end
        @test all(abs.(err_σ) .<= 0.01)
        # We return this last one instead of testing it directly because it is broken for Hakulinen, so we can test_broken it only for haku. 
        return all(abs.(err_S) .<= 0.01) 
    end
    args =  (colrec, slopop, colrec.time, colrec.status, colrec.age, colrec.year, colrec.sex)
    @test        test_surv("pohar-perme", PoharPerme, args...)
    @test        test_surv("ederer1",     EdererI,    args...)
    @test        test_surv("ederer2",     EdererII,   args...)

    # Hakulinen is unfortunately broken.
    @test_broken test_surv("hakulinen",   Hakulinen,  colrec, slopop)
end

@testitem "Check Pohar Perme vs the truth on a large simulated dataset." begin 
    using RateTables, NetSurvival, RCall, DataFrames
    const YEARS = 365.241
    function sampler(age_min, age_max)
        # Population mortality distribution: 
        age  = (age_min + rand() * (age_max-age_min))YEARS
        year = rand(1990:2009)YEARS
        sex  = rand((:male,:female))
        P    = Life(slopop[sex], age, year)

        # A bit of random censoring + administrative censoring at 15 years: 
        e = rand(Exponential(10YEARS))
        p = rand(P)
        c = min(rand(Exponential(20YEARS)), 15YEARS)
        t = min(e,p,c)
        δ = c > t # 0 means censored
        return (time=t, status=δ, age=age, year=year, sex=sex,)
    end

    # A test on a large sample to check we match the truth (convergence of Pohar Perme)
    n=15000
    df = DataFrame((sampler(35, 75) for i in 1:n)) 
    model = fit(PoharPerme, @formula(Surv(time,status)~1), df, slopop)
    truth = ccdf(Exponential(10YEARS), model.grid)
    @test maximum(abs.(model.Sₑ .- truth)) < 0.05 # coresponds to a KS test, to esnure that we match the truth. 
    

    # A test on a smaller sample to check again that we match relsurv on a simulated dataset. 
    n=5000
    df = DataFrame((sampler(35, 75) for i in 1:n)) 
    model = fit(PoharPerme, @formula(Surv(time,status)~1), df, slopop)
    truth = ccdf(Exponential(10YEARS), model.grid)
    @rput df
    R"""
        df$year = as.Date(paste(trunc(df$year/365.241), 01, 01), "%Y %m %d")
        df$sex = as.integer(df$sex == "male")+1
        df$status=as.integer(df$status)
        rez = relsurv::rs.surv(survival::Surv(time, status) ~ 1, 
        rmap=list(age = age, sex =sex, year = year), 
        data = df, 
        ratetable = relsurv::slopop, 
        method = "pohar-perme", 
        add.times=1:max(df$time))
        rr = rez$surv
    """
    r_est = @rget rr
    @test mean(abs.(model.Sₑ .- r_est)) < 0.01
end


@testitem "Assess GraffeoTest" begin
    
    using RCall, RateTables, DataFrames

    # Prerequisite functions: 
    function check_equal(test1, test2)
        @test all(test1.∂N .== test2.∂N)
        @test all(test1.∂V .== test2.∂V)
        @test all(test1.∂Z .== test2.∂Z)
        @test test1.stat == test2.stat
        @test test1.df == test2.df
        @test test1.pval == test2.pval
        return nothing
    end
    function check_no_nan(instance)
        @test !any(isnan.(instance.∂N))
        @test !any(isnan.(instance.∂V))
        @test !any(isnan.(instance.∂Z))
        @test !isnan(instance.stat)
        @test !isnan(instance.df)
        @test !isnan(instance.pval)
        return nothing
    end
    function compare_with_R(test,r)
        err_F = (r[:test_stat] - test.stat) / r[:test_stat]
        err_p = r[:p_value] == 0.0 ? 0.0 : (r[:p_value] - test.pval) / r[:p_value]
        err_df = (r[:df] - test.df) / r[:df]
        @test abs(err_F)   <= 0.02
        @test abs(err_p)   <= 0.001
        @test abs.(err_df) <= 0.01
        return nothing
    end

    # Make the different tests: 
    v1 = fit(GraffeoTest, @formula(Surv(time,status)~stage), colrec, slopop)
    v2 = GraffeoTest(colrec.time, colrec.status, colrec.age, colrec.year, colrec.sex, ones(length(colrec.age)), colrec.stage, slopop)

    v1_strat = fit(GraffeoTest, @formula(Surv(time,status)~stage+Strata(sex)), colrec, slopop)
    v2_strat = GraffeoTest(colrec.time, colrec.status, colrec.age, colrec.year, colrec.sex, colrec.sex, [join(row, " ") for row in eachrow(select(colrec,[:sex,:stage]))], slopop)

    R"""
    rez = relsurv::rs.diff(survival::Surv(time, stat) ~ stage, rmap=list(age = age, sex = sex, year = diag), data = relsurv::colrec, ratetable = relsurv::slopop)
    rez_strat = relsurv::rs.diff(survival::Surv(time, stat) ~ stage+survival::strata(sex), rmap=list(age = age, sex = sex, year = diag), data = relsurv::colrec, ratetable = relsurv::slopop)
    """
    vR = @rget rez
    vR_strat = @rget rez_strat

    # Check for no-nans: 
    check_no_nan(v1)
    check_no_nan(v2)
    check_no_nan(v1_strat)
    check_no_nan(v2_strat)

    # Coompare results with R: 
    compare_with_R(v1, vR)
    compare_with_R(v1_strat, vR_strat)

    # Check for equality of the two interfaces: 
    check_equal(v1,v2)
    check_equal(v1_strat,v2_strat)
end


@testitem "Assess CrudeMortality" begin

    using DataFrames
    using RateTables
    using RCall

    colrec.country = rand(keys(hmd_countries),nrow(colrec))
    instance = fit(CrudeMortality, @formula(Surv(time,status)~1), colrec, slopop)

    # Check that this verison is returning the same thing: 
    v2 = CrudeMortality(fit(EdererII, @formula(Surv(time, status)~1), colrec, slopop))
    @test all(instance.Mₒ .== v2.Mₒ)
    @test all(instance.Mₑ .== v2.Mₑ)
    @test all(instance.Mₚ .== v2.Mₚ)

    # Check for no nans: 
    @test !any(isnan.(instance.Mₒ))
    @test !any(isnan.(instance.Mₑ))
    @test !any(isnan.(instance.Mₚ))

    # Now compare with R baseline: 
    R"""
    rez = relsurv::cmp.rel(survival::Surv(time, stat) ~ 1, rmap=list(age = age, sex = sex, year = diag), data = relsurv::colrec, ratetable = relsurv::slopop)
    """
    r = @rget rez

    err_causeSpec = (r[:causeSpec][:est][2:end, :]  .- instance.Mₑ[1:end, :]) ./ r[:causeSpec][:est][2:end, :]
    err_pop =       (r[:population][:est][2:end, :] .- instance.Mₚ[1:end, :]) ./ r[:population][:est][2:end, :]
    @test all(abs.(err_causeSpec) .<= 0.01)
    @test all(abs.(err_pop)       .<= 0.01)
end

@testitem "Assess Nessie" begin
    using RateTables
    using RCall
    
    R"""
    rez = relsurv::nessie(survival::Surv(time, stat) ~ sex, data = relsurv::colrec, ratetable = relsurv::slopop, rmap = list(age = age, sex = sex, year = diag))
    mata = t(as.matrix(rez$mata))
    povp = rez$povp
    """
    r_mata = @rget mata
    r_povp = @rget povp
    r_male, r_female = r_mata[:,1], r_mata[:,2]

    instance = nessie(@formula(Surv(time,status)~sex), colrec, slopop)
    jl_male, jl_female = instance[2].expected_sample_size
    jl_povp = instance[1].expected_life_time

    err_male = (r_male[1:end-1]  .- jl_male) ./ r_male[1:end-1]
    err_female = (r_female[1:end-1]  .- jl_female) ./ r_female[1:end-1]
    err_povp = (r_povp  .- jl_povp) ./ r_povp

    @test all(abs.(err_male)   .<= 0.01)
    @test all(abs.(err_female) .<= 0.01)
    @test all(abs.(err_povp)   .<= 0.01)
end


@testitem "Asses GenPPE" begin
    using RateTables
    using Copulas


    TruePP = GenPoharPerme()
    MockPP = reinterpret(NetSurvival.GenPoharPermeMethod{IndependentCopula{2}},(IndependentCopula(2),))
    
    m0 = fit(TruePP, @formula(Surv(time,status)~1), colrec, slopop)
    m1 = fit(MockPP, @formula(Surv(time,status)~1), colrec, slopop)
    @assert m0.Sₑ ≈ m1.Sₑ
    @assert m0.σₑ ≈ m1.σₑ

    t0 = fit(GraffeoTest, @formula(Surv(time,status)~stage), colrec, slopop)
    t1 = fit(GraffeoTest(MockPP), @formula(Surv(time,status)~stage), colrec, slopop)
    @assert t0.stat ≈ t1.stat
end