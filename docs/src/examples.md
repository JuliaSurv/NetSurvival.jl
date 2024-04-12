```@meta
CurrentModule = NetSurvival
```

# Examples

In this section, we will be showcasing an example on how the key functions in the package work (with a comparaison to the results obtained in the `relsurv` package in R). 

In this example, we will be using the dataset `colrec` which refers to patients with colon and rectal cancer diagnosed in 1994-2000. By loading the `slopop` rate table based on the Slovenian population, we will be able to apply the Pohar Perme estimator as well as the Grafféo log-rank test for net survival analysis purposes.

```@example
using NetSurvival

pp = fit(PoharPerme, @formula(Surv(time,status)~1), colrec, slopop)
```

```@example
using NetSurvival

test1 = fit(PoharPerme, @formula(Surv(time,status)~sex), colrec, slopop)
```

```@example
using NetSurvival

test2 = fit(GraffeoTest, @formula(Surv(time,status)~stage), colrec, slopop)
```

```@example
using NetSurvival

test3 = fit(GraffeoTest, @formula(Surv(time,status)~stage+sex), colrec, slopop)
```

```@example
using NetSurvival

test4 = fit(GraffeoTest, @formula(Surv(time,status)~stage+Strata(sex)), colrec, slopop)
```

```@example
using NetSurvival

test5= fit(GraffeoTest, @formula(Surv(time,status)~stage+Strata(sex)+Strata(site)), colrec, slopop)
```