@testitem "Assess all NPNSEstimators" begin
    using RateTables, NetSurvival, RCall, DataFrames
    function mk_r_model(r_method)

        return t::Vector{Float64},s::Vector{Float64},e::Vector{Float64}
    end
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
    compare_with_R(v1_strat, vR_strat) # <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<------------------- This ones fails 

    # Check for equality of the two interfaces: 
    check_equal(v1,v2)
    check_equal(v1_strat,v2_strat) # <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<------------------- This ones fails 
end


@testitem "Assess CrudeMortality" begin

    using DataFrames
    using RateTables
    using RCall

    colrec.country = rand(keys(hmd_countries),nrow(colrec))
    instance = fit(CrudeMortality, @formula(Surv(time,status)~1), colrec, slopop)

    # Check that this verison is returning the same thing: 
    v2 = CrudeMortality(fit(EdererII, @formula(Surv(time, status)~1), colrec, slopop))
    @test all(instance.Λₒ .== v2.Λₒ)
    @test all(instance.Λₑ .== v2.Λₑ)
    @test all(instance.Λₚ .== v2.Λₚ)

    # Check for no nans: 
    @test !any(isnan.(instance.Λₒ))
    @test !any(isnan.(instance.Λₑ))
    @test !any(isnan.(instance.Λₚ))

    # Now compare with R baseline: 
    R"""
    rez = relsurv::cmp.rel(survival::Surv(time, stat) ~ 1, rmap=list(age = age, sex = sex, year = diag), data = relsurv::colrec, ratetable = relsurv::slopop)
    """
    r = @rget rez

    err_causeSpec = (r[:causeSpec][:est][2:end, :]  .- instance.Λₑ[1:end, :]) ./ r[:causeSpec][:est][2:end, :]
    err_pop =       (r[:population][:est][2:end, :] .- instance.Λₚ[1:end, :]) ./ r[:population][:est][2:end, :]
    @test all(abs.(err_causeSpec) .<= 0.01)
    @test all(abs.(err_pop)       .<= 0.01)
end