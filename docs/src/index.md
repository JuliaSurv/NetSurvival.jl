```@meta
CurrentModule = NetSurvival
```

## Introduction


The `NetSurvival.jl` package provides the necessary tools to perform estimations and analysis in the Net Survival field, a specialized branch of Survival Analysis focused on estimating the probability of survival from a specific event of interest, for example a given cancer, without taking into account other causes of mortality, in the (unfortunately quite common) case where the cause of death indicatrix is *not* available (or not trustable). This so-called *missing indicatrix* issue forbids the use of standard competitive risks survival analysis methods on these datasets.  

By integrating observed data of the target population with historical population mortality data (usually grabbed from national census datasets), Net Survival allows the extraction of the specific mortality hazard associated with the particular disease, even under the missing indicatrix issue. The concept of relative survival analysis dates back several decades to the seminal article by Ederer, Axtell, and Cutler in 1961 [Ederer1961](@cite) and the one by Ederer and Heise in 1959 [Ederer1959](@cite).

For years, the Hakulinen estimator (1987) [Hakulinen1977](@cite) and the Ederer I and II estimators were widely regarded as the gold standard for non-parametric survival curve estimation. However, the introduction of the Pohar-Perme, Stare, and Estève estimator in 2012 [PoharPerme2012](@cite) resolved several issues inherent in previous estimators, providing a reliable and consistent non-parametric estimator for net survival analysis.

This package is an attempt to bring standard relative survival analysis modeling routines in Julia, while providing an interface that is close to the `relsurv` standard while being much much faster and easier to maintain in the future. 

Some key features in `NetSurvival.jl` are:

- A panel of different non-parametric net survival estimators (Ederer I [Ederer1961](@cite), Ederer II [Ederer1959](@cite), Hakulinen [Hakulinen1977](@cite), Pohar Perme [PoharPerme2012](@cite)) with an interface that comply to Julia's standards. 
- Grafféo's log-rank test [Graffeo2016](@cite) to compare net survival curves accross groups, including stratified testing.
- Small, readable and efficient codebase (up to 1000x less LOC than `relsurv` for the same functionalities), which allows long-term maintenability of the code.
- Large performance improvements (up to 50x) relative to the R package `relsurv`

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
