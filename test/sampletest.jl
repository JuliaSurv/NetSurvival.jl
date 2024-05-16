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
    
    # R version
    using RCall
    using RateTables
    R"""
    rez = relsurv::rs.diff(survival::Surv(time, stat) ~ stage, rmap=list(age = age, sex = sex, year = diag), data = relsurv::colrec, ratetable = relsurv::slopop)
    """
    r = @rget rez

    # Julia version
    instance = fit(GraffeoTest, @formula(Surv(time,status)~stage), colrec, slopop)

    # Check concordance between the two calling syntaxes: 
    v2 = GraffeoTest(colrec.time, colrec.status, colrec.age, colrec.year, colrec.sex, ones(length(colrec.age)), colrec.stage, slopop)
    @test all(instance.∂N .== v2.∂N)
    @test all(instance.∂V .== v2.∂V)
    @test all(instance.∂Z .== v2.∂Z)
    @test instance.stat == v2.stat
    @test instance.df == v2.df
    @test instance.pval == v2.pval


    # Check for no nans: 
    @test !any(isnan.(instance.∂N))
    @test !any(isnan.(instance.∂V))
    @test !any(isnan.(instance.∂Z))
    @test !isnan(instance.stat)
    @test !isnan(instance.df)
    @test !isnan(instance.pval)

    # Check also for no nans if strata is used: 
    second_instance = fit(GraffeoTest, @formula(Surv(time,status)~stage+Strata(sex)), colrec, frpop)
    @test !any(isnan.(second_instance.∂N))
    @test !any(isnan.(second_instance.∂V))
    @test !any(isnan.(second_instance.∂Z))
    @test !isnan(second_instance.stat)
    @test !isnan(second_instance.df)
    @test !isnan(second_instance.pval)

    # Check concordance between the two calling syntaxes, stratified
    second_v2 = GraffeoTest(colrec.time, colrec.status, colrec.age, colrec.year, colrec.sex, colrec.sex, colrec.stage, slopop)
    @test all(second_instance.∂N .== second_v2.∂N)
    @test all(second_instance.∂V .== second_v2.∂V)
    @test all(second_instance.∂Z .== second_v2.∂Z)
    @test second_instance.stat == second_v2.stat
    @test second_instance.df == second_v2.df
    @test second_instance.pval == second_v2.pval
    
    err_F = (r[:test_stat] - instance.stat) / r[:test_stat]
    err_p = r[:p_value] == 0.0 ? 0.0 : (r[:p_value] - instance.pval) / r[:p_value]
    err_df = (r[:df] - instance.df) / r[:df]

    @test abs(err_F)   <= 0.01
    @test abs(err_p)   <= 0.001
    @test abs.(err_df) <= 0.01
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