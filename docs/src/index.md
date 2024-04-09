```@meta
CurrentModule = NetSurvival
```

## Introduction

This package serves to provide the necessary tools to perform net survival analysis, a branch of survival analysis dedicated to estimating the probability of survival from a particular event of interest compared to the general public. Some key features in `NetSurvival.jl` are:

- Fitting different non-parametric estimators (Pohar Perme[PoharPerme2012](@cite), Ederer II, ...)
- Applying Grafféo's log-rank test[GraffeoTest](@cite) on different groups, including stratified covariables 
- ... 

## Installation

The package is available on Julia's general registry, and can be installed either with the command `Pkg.add("NetSurvival")` or via the Pkg REPL mode: 

```julia
] add NetSurvival
```

```@example
1 == 1
```

```@example 1
a = 2
```

```@example 1
a
```


```@index
```

```@autodocs
Modules = [NetSurvival]
```

```@bibliography
Pages = ["index.md"]
Canonical = false
```