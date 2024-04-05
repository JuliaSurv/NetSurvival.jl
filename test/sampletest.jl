# This file could be removed =, but here is a sample test: 

@testitem "trivial test" begin
    @test true
end


@testitem "trivial test 2" begin
    
    ################################
    # Run code: estimator and confidence interval. 
    instance = PoharPerme(colrec.time, colrec.status, colrec.age, colrec.year, colrec.sex, slopop)
    conf_int = confint(instance; level = 0.05)
    # Run code: test. 
    strata = colrec.sex
    group = colrec.stage
    mytest = GraffeoTest(colrec.time, colrec.status, colrec.age, colrec.year, colrec.sex, strata, group, slopop)

    fit(PoharPerme, @formula(Surv(time,status)~sex), colrec, frpop)

    @test true

    # How to use R and R packages in here ? 
end


@testitem "trying RCall" begin
    
    using RCall

    # R version
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
    instance = PoharPerme(colrec.time, colrec.status, colrec.age, colrec.year, colrec.sex, slopop)
    conf_int = confint(instance; level = 0.05)
    lower_bounds = [lower[1] for lower in conf_int]
    upper_bounds = [upper[2] for upper in conf_int]
    
    err_S = (R_Sₑ .- instance.Sₑ[1:end-1]) ./ R_Sₑ
    err_σ = (R_σ .- instance.σₑ[1:end-1]) ./ R_σ

    @test all(abs.(err_S) .<= 0.01)
    @test all(abs.(err_σ) .<= 0.01)

end