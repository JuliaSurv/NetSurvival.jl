# NetSurvival

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://JuliaSurv.github.io/NetSurvival.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://JuliaSurv.github.io/NetSurvival.jl/dev/)
[![Build Status](https://github.com/JuliaSurv/NetSurvival.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/JuliaSurv/NetSurvival.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Benchmark](https://github.com/JuliaSurv/NetSurvival.jl/actions/workflows/Benchmark.yml/badge.svg?branch=main)](https://JuliaSurv.github.io/NetSurvival.jl/benchmarks)
[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![PkgEval](https://JuliaCI.github.io/NanosoldierReports/pkgeval_badges/N/NetSurvival.svg)](https://JuliaCI.github.io/NanosoldierReports/pkgeval_badges/N/NetSurvival.html)
[![Aqua QA](https://raw.githubusercontent.com/JuliaTesting/Aqua.jl/master/badge.svg)](https://github.com/JuliaTesting/Aqua.jl)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![ColPrac: Contributor's Guide on Collaborative Practices for Community Packages](https://img.shields.io/badge/ColPrac-Contributor's%20Guide-blueviolet)](https://github.com/SciML/ColPrac)
[![Code Style: Blue](https://img.shields.io/badge/code%20style-blue-4495d1.svg)](https://github.com/JuliaDiff/BlueStyle)


The `NetSurvival.jl` package provides the necessary tools to perform estimations and analysis in the Net Survival field. This specialized branch of Survival Analysis focuses on estimating the probability of survival from a specific event of interest, for example a given cancer, without considering other causes of death. This is especially relevant in the (unfortunately quite common) case where the cause of death indicatrix is either unavailable or untrustworthy. Consequently, the so-called *missing indicatrix* issue forbids the use of standard competitive risks survival analysis methods on these datasets.  For that, a few standard estimators were established in the last 50 years, backed by a wide literature.

# Features 

This package is an attempt to bring standard relative survival analysis modeling routines to Julia, while providing an interface that is close to the `relsurv` standard, albeit significantly faster and easier to maintain in the future.

Some key features in `NetSurvival.jl` are:

- A panel of different non-parametric net survival estimators (Ederer I, Ederer II, Hakulinen, Pohar Perme) with an interface compliant with Julia's standards. 
- Grafféo's log-rank test to compare net survival curves accross groups, including stratified testing.
- A compact, readable and efficient codebase (up to 1000x less LOC than `relsurv` for the same functionalities), ensuring long-term maintenability.
- Significant performance improvements (up to 50x) compared to the R package `relsurv`.

# Getting Started

The package is not yet available on Julia's general registry, and thus can be installed through the following command:

```julia
using Pkg
Pkg.add("https://github.com/JuliaSurv/NetSurvival.jl.git")
```

See the rest of this [documentation](https://juliasurv.github.io/NetSurvival.jl/dev/) to have a glimpse of the functionalities!

