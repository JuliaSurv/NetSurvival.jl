```@meta
CurrentModule = NetSurvival
```

# Benchmarking results

The following benchmarks are run on github actions continuous integration platform, which is a very slow computing engine. Local experiments suggests performances that are twice as fast on correct hardware -- note that we do not use multithreading at all, but underlying BLAS calls might. 

```@example 1
using RCall
using NetSurvival, RateTables, BenchmarkTools

R_bench = @benchmark R"""
relsurv::rs.surv(
    survival::Surv(time, stat) ~1, 
    rmap=list(age = age, sex = sex, year = diag), 
    data = relsurv::colrec, 
    ratetable = relsurv::slopop, 
    method = "pohar-perme", 
    add.times=1:8149)
"""

jl_bench = @benchmark fit(PoharPerme, @formula(Surv(time,status)~1), colrec, slopop)

ratio = time(minimum(R_bench)) / time(minimum(jl_bench))
```


# Benchmarking across time

The folloiwng charts provide a glimpse of `NetSurvival.jl`'s performance along time: 

```@raw html
<iframe src="../../benchmarks/" style="height:500px;width:100%;"></iframe>
```
