@testitem "Interface test & ensure no nans" begin
    
    ################################
    # Run code: test. 
    using RateTables

    function mktest(::Type{E}, df, rt, args...) where E
        npe = fit(E, @formula(Surv(time,status)~1), df, rt)
        ci = confint(npe; level = 0.05)
        # Check for no-nans :
        rez = true 
        rez &= !any(isnan.(npe.Sₑ))
        rez &= !any(isnan.(npe.∂Λₑ))
        rez &= !any(isnan.(npe.σₑ))
        rez &= !any([any(isnan.(x)) for x in ci])

        # Check for correspondance of the alternative call syntax: 
        v2 = E(args..., rt)
        rez &= all(v2.Sₑ .== npe.Sₑ) & all(v2.∂Λₑ .== npe.∂Λₑ) & all(v2.σₑ .== npe.σₑ)
        return rez
    end
    args =  (colrec, slopop, colrec.time, colrec.status, colrec.age, colrec.year, colrec.sex)
    @test mktest(PoharPerme, args...)
    @test mktest(EdererI, args...)
    @test mktest(EdererII, args...)
    @test mktest(Hakulinen, args...)
    
    rez = fit(GraffeoTest, @formula(Surv(time,status)~stage+Strata(sex)), colrec, frpop)
    @test !isnan(rez.pval) && !isnan(rez.stat)
end

@testitem "Compare NPNSEstimator's with relsurv::rs.surv on colrec x slopop" begin
    using RateTables, NetSurvival, RCall, DataFrames
    function mk_r_model(r_method)
        @rput r_method
        R"""
        rez = relsurv::rs.surv(
            survival::Surv(time, stat) ~ 1, 
            rmap=list(age = age, sex = sex, year = diag), 
            data = relsurv::colrec, 
            ratetable = relsurv::slopop, 
            method = r_method, 
            add.times=1:8149)
        t = rez$time
        s = rez$surv
        e = rez$std.err
        """
        t = @rget t
        s = @rget s
        e = @rget e
        return t::Vector{Float64},s::Vector{Float64},e::Vector{Float64}
    end
    function test_surv(r_method,::Type{E},df, rt) where E
        jl = fit(E, @formula(Surv(time,status)~1), df, rt)
        r_t, r_S, r_σ = mk_r_model(r_method)

        err_S = zeros(length(r_t))
        err_σ = zeros(length(r_t))
        for i in eachindex(r_t) 
            j = searchsortedlast(jl.grid, r_t[i])
            err_S[i] = (r_S[i] - jl.Sₑ[j]) / r_S[i]
            err_σ[i] = (r_σ[i] - jl.σₑ[j]) / r_σ[i]
        end
        return all(abs.(err_S) .<= 0.01) && all(abs.(err_σ) .<= 0.01)
    end

    @test        test_surv("pohar-perme", PoharPerme, colrec, slopop)
    @test        test_surv("ederer1",     EdererI,    colrec, slopop)
    @test        test_surv("ederer2",     EdererII,   colrec, slopop)
    @test_broken test_surv("hakulinen",   Hakulinen,  colrec, slopop)
end


@testitem "Comparing log rank test with R" begin
    
    # R version
    using RCall
    using RateTables
    R"""
    rez = relsurv::rs.diff(survival::Surv(time, stat) ~ stage, rmap=list(age = age, sex = sex, year = diag), data = relsurv::colrec, ratetable = relsurv::slopop)
    """
    R_model = @rget rez
    R_test = R_model[:test_stat]
    R_pvalue = R_model[:p_value]
    R_df = R_model[:df]

    # Julia version
    graffeo = fit(GraffeoTest, @formula(Surv(time,status)~stage), colrec, slopop)

    other_calling_syntax = GraffeoTest(colrec.time, colrec.status, colrec.age, colrec.year, colrec.sex, colrec.sex, colrec.stage, slopop)
    
    err_F = (R_test - graffeo.stat) / R_test
    err_p = R_pvalue == 0.0 ? 0.0 : (R_pvalue - graffeo.pval) / R_pvalue
    err_df = (R_df - graffeo.df) / R_df

    @test all(abs.(err_F) .<= 0.01)
    @test all(abs.(err_p) .<= 0.001)
    @test all(abs.(err_df) .<= 0.01)

end

@testitem "crude mortality interface" begin
    
    ################################
    # Run code: test. 
    using DataFrames
    using RateTables

    colrec.country = rand(keys(hmd_countries),nrow(colrec))
    fit(CrudeMortality, @formula(Surv(time,status)~1), colrec, slopop)
    CrudeMortality(fit(EdererII, @formula(Surv(time, status)~1), colrec, slopop))

    @test true
end

@testitem "comparing crude mortality" begin
    
    ################################
    # Run code: test. 
    using DataFrames
    using RateTables

    colrec.country = rand(keys(hmd_countries),nrow(colrec))
    instance = fit(CrudeMortality, @formula(Surv(time,status)~1), colrec, slopop)

    # R version
    using RCall
    using RateTables
    R"""
    rez = relsurv::cmp.rel(survival::Surv(time, stat) ~ 1, rmap=list(age = age, sex = sex, year = diag), data = relsurv::colrec, ratetable = relsurv::slopop)
    """
    R_model = @rget rez
    R_grid = R_model[:causeSpec][:time]
    R_causeSpec = R_model[:causeSpec][:est]
    R_population = R_model[:population][:est]

    err_causeSpec = (R_causeSpec[2:end, :] .- instance.Λₑ[1:end, :]) ./ R_causeSpec[2:end, :]
    err_pop = (R_population[2:end, :] .- instance.Λₚ[1:end, :]) ./ R_population[2:end, :]

    @test all(abs.(err_causeSpec) .<= 0.01)
    @test all(abs.(err_pop) .<= 0.01)
end