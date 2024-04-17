```@meta
CurrentModule = NetSurvival
```

# Getting Started

## Fitting the non parametric estimators

```@docs
PoharPerme
```

### Example

We will illustrate with an example using the dataset `colrec`, which comprises $5971$ patients diagnosed with colon or rectal cancer  between 1994 and 2000. This dataset is sourced from the Slovenia cancer registry. Given the high probability that the patients are Slovenian, we will be using the Slovenian mortality table `slopop` as reference for the populational rates. Subsequently, we can apply various non-parametric estimators for net survival analysis.

!!! note "N.B." 
    Mortality tables may vary in structure, with options such as the addition or removal of specific covariates. To confirm that the mortality table is in the correct format, please refer to the documentation of `RateTables.jl`, or directly extract it from there.

By examining `slopop`, we notice it contains information regarding `age` and `year`, as expected for mortality tables. Additionally, it incorporates the covariate sex, which has two possible entries (`:male` or `:female`).

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

The p-value is well under $0.05$, meaning that the different groups identified by the `stage` variable have different survival probabilities. Thus, it should be taken into consideration in the study.

```@example 1
test2 = fit(GraffeoTest, @formula(Surv(time,status)~sex), colrec, slopop)
```

For the `sex` variable, we notice that the p-value is above $0.05$ indicating that there isn't a difference between male and female patients.

```@example 1
test4 = fit(GraffeoTest, @formula(Surv(time,status)~stage+Strata(sex)), colrec, slopop)
```