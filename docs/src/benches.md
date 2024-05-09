```@meta
CurrentModule = NetSurvival
```

# Benchmarking results

This page provides benchmark results of several standard net survival routines implemented in this package. Note that the runtime also depends on other packages, in particular on [`RateTables.jl`](https://github.com/JuliaSurv/RateTables.jl). 

!!! note "Take numbers displayed here carefully"
    All the following benchmarks are run on github action contunous integration platform, which is a very slow computing engine. Thus, these numbers may not represent your local performance. A locally ran version of these benchmarks in availiable on the github readme, and the below code blocks can be used to check performance on your own hardware.

## Benchmarks w.r.t. `relsurv`

This first set of benchmark compares standard functionalities with their implementation in `relsurv`. Below numbers gives runtime mulitpliers w.r.t. [`R::relsurv`](https://cran.r-project.org/web/packages/relsurv/index.html), computed on github action CI.

```@example 1
using RateTables, NetSurvival, RCall, DataFrames

function test_surv(r_method,::Type{E}, stratified) where E
    if stratified
        jl = @timed fit(E, @formula(Surv(time,status)~sex), colrec, slopop)
        @rput r_method
        r = @timed R"""
            rez = relsurv::rs.surv(survival::Surv(time, stat) ~ sex, rmap=list(age = age, sex = sex, year = diag), data = relsurv::colrec, ratetable = relsurv::slopop, method = r_method, add.times=1:8149)
        """
    else
        jl = @timed fit(E, @formula(Surv(time,status)~1), colrec, slopop)
        @rput r_method
        r = @timed R"""
            rez = relsurv::rs.surv(survival::Surv(time, stat) ~ 1, rmap=list(age = age, sex = sex, year = diag), data = relsurv::colrec, ratetable = relsurv::slopop, method = r_method, add.times=1:8149)
        """
    end
    return r.time / jl.time
end

function test_graffeo(stratified)
    if stratified
        jl = @timed fit(GraffeoTest, @formula(Surv(time,status)~stage+Strata(sex)), colrec, slopop)
        r = @timed R"""
        rez = relsurv::rs.diff(survival::Surv(time, stat) ~ stage + survival::strata(sex), rmap=list(age = age, sex = sex, year = diag), data = relsurv::colrec, ratetable = relsurv::slopop)
        """
    else 
        jl = @timed fit(GraffeoTest, @formula(Surv(time,status)~stage), colrec, slopop)
        r = @timed R"""
        rez = relsurv::rs.diff(survival::Surv(time, stat) ~ stage, rmap=list(age = age, sex = sex, year = diag), data = relsurv::colrec, ratetable = relsurv::slopop)
        """
    end
    return r.time / jl.time 
end
test_all(stratified) = [
    test_surv("pohar-perme", PoharPerme, stratified),
    test_surv("ederer1", EdererI, stratified),
    test_surv("ederer2", EdererII, stratified),
    test_surv("hakulinen", Hakulinen, stratified),
    test_graffeo(stratified),
]
test_all() = DataFrame(
    Algorithm = ["Pohar Perme", "EdererI", "EdererII", "Hakulinen", "Graffeo's LRT"], 
    unstratified = test_all(false), 
    stratified = test_all(true)
)
rez = test_all()
rez = test_all() # discard first run.

# note: to obtain the pretty printing from the readme, you need to install PrettyTables.jl and do : 
# using PrettyTables
# pretty_table(rez, backend = Val(:markdown))

rez
```




# Benchmarking across time

The following charts provide a glimpse of `NetSurvival.jl`'s performance along time, also ran on github CI: 

```@raw html
<iframe src="../../benchmarks/" style="height:500px;width:100%;"></iframe>
```
