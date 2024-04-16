```@meta
CurrentModule = NetSurvival
```

# Getting Started

## Fitting the non parametric estimators

```@docs
PoharPerme
```

### Example

We will demonstrate through the use of an example. We will be using the dataset `colrec` which refers to patients with colon and rectal cancer diagnosed in 1994-2000, taken from the Slovenia cancer registry. By loading the `slopop` rate table based on the Slovenian population, we will be able to apply the different non-parametric estimators for net survival analysis purposes.


**Pohar Perme**
```@example 1
using NetSurvival, RateTables
pp1 = fit(PoharPerme, @formula(Surv(time,status)~1), colrec, slopop)
```

## Applying the Grafféo log-rank test

```@docs
GraffeoTest
```

### Example

When applying the test to the same data as before, we get:

```@example 1
test1 = fit(GraffeoTest, @formula(Surv(time,status)~stage), colrec, slopop)
```

The p-value is well under $0.05$, meaning that the different groups determined by the `stage` variable have different survival probabilities. Thus, it should be taken into consideration in the study.

```@example 1
test2 = fit(GraffeoTest, @formula(Surv(time,status)~sex), colrec, slopop)
```

For the `sex` variable, we notice that the p-value is above $0.05$ indicating that there isn't a difference between male and female.

```@example 1
test4 = fit(GraffeoTest, @formula(Surv(time,status)~stage+Strata(sex)), colrec, slopop)
```