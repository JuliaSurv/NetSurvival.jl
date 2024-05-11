@testitem "Basic functionality test" begin
    using RateTables
    ################################
    instance = PoharPerme(colrec.time, colrec.status, colrec.age, colrec.year, colrec.sex, slopop)
    conf_int = confint(instance; level = 0.05)
    mytest   = GraffeoTest(colrec.time, colrec.status, colrec.age, colrec.year, colrec.sex, colrec.sex, colrec.stage, slopop)
    fit(PoharPerme, @formula(Surv(time,status)~sex), colrec, frpop)

    @test true
end

@testitem "Interface test & No-Nan" begin
    
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

@testitem "Compare NPNSEstimator's with R::relsurv::rs.surv on colrec x slopop" begin
    using RateTables, NetSurvival, RCall, DataFrames
    function test_surv(r_method,::Type{E}) where E
        @rput r_method
        jl = fit(E, @formula(Surv(time,status)~1), colrec, slopop)
        R"""
            rez = relsurv::rs.surv(
                survival::Surv(time, stat) ~ 1, 
                rmap=list(age = age, sex = sex, year = diag), 
                data = relsurv::colrec, 
                ratetable = relsurv::slopop, 
                method = r_method, 
                add.times=1:8149)
        """
        R_model = @rget rez
        err_S = (R_model[:surv] .- jl.Sₑ[1:end-1]) ./ R_model[:surv]
        err_σ = (R_model[:std_err] .- jl.σₑ[1:end-1]) ./ R_model[:std_err]
        return all(abs.(err_S) .<= 0.01) && all(abs.(err_σ) .<= 0.01)
    end

    @test        test_surv("pohar-perme", PoharPerme)
    @test        test_surv("ederer1",     EdererI)
    @test        test_surv("ederer2",     EdererII)
    @test_broken test_surv("hakulinen",   Hakulinen)
end


@testitem "Comparing PoharPerme, Ederer1, Ederer2 and Hakulinen with R" begin
    
    # R version
    using RCall
    using RateTables
    R"""
    rez = relsurv::rs.surv(
        survival::Surv(time, stat) ~1, 
        rmap=list(age = age, sex = sex, year = diag), 
        data = relsurv::colrec, 
        ratetable = relsurv::slopop, 
        method = "pohar-perme", 
        add.times=1:8149)
    """
    R_model = @rget rez
    R_grid = R_model[:time]
    R_Sₑ = R_model[:surv]
    R_σ = R_model[:std_err]
    R_low = R_model[:lower]
    R_upp = R_model[:upper]


    # Julia version
    instance = fit(PoharPerme, @formula(Surv(time,status)~1), colrec, slopop)
    conf_int = confint(instance; level = 0.05)
    lower_bounds = [lower[1] for lower in conf_int]
    upper_bounds = [upper[2] for upper in conf_int]
    
    err_S = (R_Sₑ .- instance.Sₑ[1:end-1]) ./ R_Sₑ
    err_σ = (R_σ .- instance.σₑ[1:end-1]) ./ R_σ

    @test all(abs.(err_S) .<= 0.01)
    @test all(abs.(err_σ) .<= 0.01)

end

@testitem "Comparing EdererI with R" begin
    
    # R version
    using RCall
    using RateTables
    R"""
    rez = relsurv::rs.surv(
        survival::Surv(time, stat) ~1, 
        rmap=list(age = age, sex = sex, year = diag), 
        data = relsurv::colrec, 
        ratetable = relsurv::slopop, 
        method = "ederer1", 
        add.times=1:8149)
    """
    R_model = @rget rez
    R_grid = R_model[:time]
    R_Sₑ = R_model[:surv]
    R_σ = R_model[:std_err]
    R_low = R_model[:lower]
    R_upp = R_model[:upper]


    # Julia version
    instance = fit(EdererI, @formula(Surv(time,status)~1), colrec, slopop)
    conf_int = confint(instance; level = 0.05)
    lower_bounds = [lower[1] for lower in conf_int]
    upper_bounds = [upper[2] for upper in conf_int]

    err_S = zeros(length(R_grid))
    err_σ = zeros(length(R_grid))
    
    for i in eachindex(R_grid) 
        for j in eachindex(instance.grid)
            if R_grid[i] == instance.grid[j]
                err_S[i] += (R_Sₑ[i] - instance.Sₑ[j]) / R_Sₑ[i]
                err_σ[i] += (R_σ[i] - instance.σₑ[j]) / R_σ[i]
            end
        end
    end

    @test all(abs.(err_S) .<= 0.01)
    @test all(abs.(err_σ) .<= 0.01)

end

@testitem "Comparing EdererII with R" begin
    
    # R version
    using RCall
    using RateTables
    R"""
    rez = relsurv::rs.surv(
        survival::Surv(time, stat) ~1, 
        rmap=list(age = age, sex = sex, year = diag), 
        data = relsurv::colrec, 
        ratetable = relsurv::slopop, 
        method = "ederer2", 
        add.times=1:8149)
    """
    R_model = @rget rez
    R_grid = R_model[:time]
    R_Sₑ = R_model[:surv]
    R_σ = R_model[:std_err]
    R_low = R_model[:lower]
    R_upp = R_model[:upper]


    # Julia version
    instance = fit(EdererII, @formula(Surv(time,status)~1), colrec, slopop)
    conf_int = confint(instance; level = 0.05)
    lower_bounds = [lower[1] for lower in conf_int]
    upper_bounds = [upper[2] for upper in conf_int]

    err_S = zeros(length(R_grid))
    err_σ = zeros(length(R_grid))
    
    for i in eachindex(R_grid) 
        for j in eachindex(instance.grid)
            if R_grid[i] == instance.grid[j]
                err_S[i] += (R_Sₑ[i] - instance.Sₑ[j]) / R_Sₑ[i]
                err_σ[i] += (R_σ[i] - instance.σₑ[j]) / R_σ[i]
            end
        end
    end

    @test all(abs.(err_S) .<= 0.01)
    @test all(abs.(err_σ) .<= 0.01)

end

@testitem "Comparing Hakulinen with R" begin
    
    # R version
    using RCall
    using RateTables
    R"""
    rez = relsurv::rs.surv(
        survival::Surv(time, stat) ~1, 
        rmap=list(age = age, sex = sex, year = diag), 
        data = relsurv::colrec, 
        ratetable = relsurv::slopop, 
        method = "hakulinen", 
        add.times=1:8149)
    """
    R_model = @rget rez
    R_grid = R_model[:time]
    R_Sₑ = R_model[:surv]
    R_σ = R_model[:std_err]
    R_low = R_model[:lower]
    R_upp = R_model[:upper]


    # Julia version
    instance = fit(Hakulinen, @formula(Surv(time,status)~1), colrec, slopop)
    conf_int = confint(instance; level = 0.05)
    lower_bounds = [lower[1] for lower in conf_int]
    upper_bounds = [upper[2] for upper in conf_int]

    err_S = zeros(length(R_grid))
    err_σ = zeros(length(R_grid))
    
    for i in eachindex(R_grid) 
        for j in eachindex(instance.grid)
            if R_grid[i] == instance.grid[j]
                err_S[i] += (R_Sₑ[i] - instance.Sₑ[j]) / R_Sₑ[i]
                err_σ[i] += (R_σ[i] - instance.σₑ[j]) / R_σ[i]
            end
        end
    end

    @test_broken all(abs.(err_S) .<= 0.01)
    @test all(abs.(err_σ) .<= 0.01)
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

    err_causeSpec = (R_causeSpec[2:end, :] .- instance.Λₑ[1:(end-1), :]) ./ R_causeSpec[2:end, :]
    err_pop = (R_population[2:end, :] .- instance.Λₚ[1:(end-1), :]) ./ R_population[2:end, :]

    @test all(abs.(err_causeSpec) .<= 0.01)
    @test all(abs.(err_pop) .<= 0.01)
end

