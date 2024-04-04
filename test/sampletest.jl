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
