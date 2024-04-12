```@meta
CurrentModule = NetSurvival
```

# Examples

In this section, we will be showcasing an example on how the key functions in the package work (with a comparaison to the results obtained in the `relsurv` package in R). 

In this example, we will be using the dataset `colrec` which refers to patients with colon and rectal cancer diagnosed in 1994-2000. By loading the `slopop` rate table based on the Slovenian population, we will be able to apply the Pohar Perme estimator as well as the Grafféo log-rank test for net survival analysis purposes.

```@example
pp1 = fit(PoharPerme, @formula(Surv(time,status)~1), colrec, slopop)
```

```@example
pp2 = fit(PoharPerme, @formula(Surv(time,status)~sex), colrec, slopop)
```

```@example
test1 = fit(GraffeoTest, @formula(Surv(time,status)~stage), colrec, slopop)
```

```@example
test2 = fit(GraffeoTest, @formula(Surv(time,status)~stage+sex), colrec, slopop)
```

```@example
test3 = fit(GraffeoTest, @formula(Surv(time,status)~stage+Strata(sex)), colrec, slopop)
```

```@example
test4 = fit(GraffeoTest, @formula(Surv(time,status)~stage+Strata(sex)+Strata(site)), colrec, slopop)
```

For the sake of comparison, the examples below detail the difference in performance between `NetSurvival.jl` and `relsurv` on R : 

```
@time fit(PoharPerme, @formula(Surv(time,status)~1), colrec, slopop);
@time R"""
rez = relsurv::rs.surv(survival::Surv(time, stat) ~1, rmap=list(age = age, sex = sex, year = diag), data = relsurv::colrec, ratetable = relsurv::slopop, method = "pohar-perme", add.times=1:8149)
"""
```

```
@time fit(GraffeoTest, @formula(Surv(time,status)~stage), colrec, slopop);
@time R"""
rez = relsurv::rs.diff(survival::Surv(time, stat) ~ stage, rmap=list(age = age, sex = sex, year = diag), data = relsurv::colrec, ratetable = relsurv::slopop)
"""
```
