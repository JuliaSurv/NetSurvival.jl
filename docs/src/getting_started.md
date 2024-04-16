```@meta
CurrentModule = NetSurvival
```

# Getting Started

## Fitting the non parametric estimators

```@docs
PoharPerme
```

We will demonstrate through the use of an example. We will be using the dataset `colrec` which refers to patients with colon and rectal cancer diagnosed in 1994-2000, taken from the Slovenia cancer registry. By loading the `slopop` rate table based on the Slovenian population, we will be able to apply the different non-parametric estimators as well as the Grafféo log-rank test for net survival analysis purposes.


```@example 1
using NetSurvival, RateTables
pp1 = fit(PoharPerme, @formula(Surv(time,status)~1), colrec, slopop)
```

## Applying the Grafféo log-rank test

```@docs
GraffeoTest
```

When applying the test to the same data as before, we get:

```@example 1
test1 = fit(GraffeoTest, @formula(Surv(time,status)~stage), colrec, slopop)
```

```@example 1
test2 = fit(GraffeoTest, @formula(Surv(time,status)~stage+sex), colrec, slopop)
```

```@example 1
test4 = fit(GraffeoTest, @formula(Surv(time,status)~stage+Strata(sex)+Strata(site)), colrec, slopop)
```