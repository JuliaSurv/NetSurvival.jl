```@meta
CurrentModule = NetSurvival
```

## Introduction

This package provides the necessary tools to perform net survival analysis, a specialized branch of survival analysis focused on estimating the probability of survival from a specific event of interest relative to the general population.

By integrating observed data with population mortality data, this method extracts insights into the hazard associated with a particular disease. The concept of relative survival analysis dates back several decades to the seminal article by Ederer, Axtell, and Cutler in 1961 [Ederer1961](@cite) and the one in 1959 by Ederer and Heise [Ederer1959](@cite). 

For years, the Hakulinen estimator (1987) and the Ederer I and II estimators were widely regarded as the gold standard for non-parametric survival curve estimation. However, the introduction of the Pohar-Perme, Stare, and Estève estimator in 2012 [PoharPerme2012](@cite) resolved several issues inherent in previous estimators, providing a reliable and consistent non-parametric estimator for net survival analysis.

Some key features in `NetSurvival.jl` are:

- Fitting different non-parametric estimators (Ederer II [Ederer1959](@cite), Pohar Perme [PoharPerme2012](@cite), ...)
- Applying Grafféo's log-rank test [GraffeoTest](@cite) on different groups, including stratified covariables 
- ... 

## Installation

The package is available on Julia's general registry, and can be installed with the following command : 

```julia
] add NetSurvival
```

```@bibliography
Pages = ["index.md"]
Canonical = false
```
