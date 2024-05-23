```@meta
CurrentModule = NetSurvival
```

## Introduction

The `NetSurvival.jl` package provides the necessary tools to perform estimations and analysis in the Net Survival field. This specialized branch of Survival Analysis focuses on estimating the probability of survival from a specific event of interest, for example a given cancer, without considering other causes of death. This is especially relevant in the (unfortunately quite common) case where the cause of death indicatrix is either unavailable or untrustworthy. Consequently, the so-called *missing indicatrix* issue forbids the use of standard competitive risks survival analysis methods on these datasets.  For that, a few standard estimators were established in the last 50 years, backed by a wide literature.

By integrating observed data from the target population with historical population mortality data (usually sourced from national census datasets), Net Survival allows the extraction of the specific mortality hazard associated with the particular disease, even under the missing indicatrix issue. The concept of relative survival analysis dates back several decades to the seminal article by Ederer, Axtell, and Cutler in 1961 [Ederer1961](@cite) and the one by Ederer and Heise in 1959 [Ederer1959](@cite).

For years, the Hakulinen estimator (1977) [Hakulinen1977](@cite) and the Ederer I and II estimators were widely regarded as the gold standard for non-parametric survival curve estimation. However, the introduction of the Pohar-Perme, Stare, and Estève estimator in 2012 [PoharPerme2012](@cite) resolved several issues inherent in previous estimators, providing a reliable and consistent non-parametric estimator for net survival analysis.

## Features

Standard tools nowadays are composed of R packages, with underlying C and C++ routines, that are hard to read, maintain, and use. This package is an attempt to bring standard relative survival analysis modeling routines to Julia, while providing an interface that is close to the `relsurv` standard, albeit significantly faster and easier to maintain in the future. Our hope is that the junction with classical modeling API in Julia will allow later extensions of the existing modeling methods, with a simple interface for the practitioners.

Some key features in `NetSurvival.jl` are:

- A panel of different non-parametric net survival estimators (Ederer I [Ederer1961](@cite), Ederer II [Ederer1959](@cite), Hakulinen [Hakulinen1977](@cite), Pohar Perme [PoharPerme2012](@cite)) with an interface compliant with Julia's standards. 
- Grafféo's log-rank test [Graffeo2016](@cite) to compare net survival curves accross groups, including stratified testing.
- Crude mortality, Expected Sample Size, and other useful metrics in net survival field.
- A 'Nessie' function that outputs the estimated sample size by yearly intervals and the average lifespan expectancy left for a given group. 
- A compact, readable and efficient codebase (up to 1000x less LOC than `relsurv` for the same functionalities), ensuring long-term maintenability.
- Significant performance improvements (up to 50x) compared to the R package `relsurv`.

## Installation

The package is not yet available on Julia's general registry, and thus can be installed through the following command:

```julia
using Pkg
Pkg.add("https://github.com/JuliaSurv/NetSurvival.jl.git")
```

See the rest of this documentation to have a glimpse of the functionalities!

# References

```@bibliography
Pages = ["index.md"]
Canonical = false
```
