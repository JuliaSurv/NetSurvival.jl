# NetSurvival.jl

*A pure-Julia take on standard net survival routines*

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


The `NetSurvival.jl` package is part of the JuliaSurv survival analysis suite. It provides the necessary tools to perform estimations and analysis in the Net Survival field. This specialized branch of Survival Analysis focuses on estimating the probability of survival from a specific event of interest, for example a given cancer, without considering other causes of death, in the (unfortunately quite common) case where the cause of death indicatrix is unavailable (or e.g. untrustworthy). Consequently, the so-called *missing indicatrix* issue forbids the use of standard competitive risks survival analysis methods on these datasets. Thus, a few standard estimators were established in the last 50 years, backed by a wide literature.

# Features 

This package is an attempt to bring standard relative survival analysis modeling routines to Julia, while providing an interface that is close to the R package `relsurv`, albeit significantly faster and easier to maintain in the future. We aim at covering the standard estimators, needed for routines and comparisons, but also to provide the most up to date state of the art. 

Some key features in `NetSurvival.jl` are:

- A panel of different non-parametric net survival estimators (Ederer I, Ederer II, Hakulinen, Pohar Perme) with an interface compliant with Julia's standards. 
- Grafféo's log-rank test to compare net survival curves accross groups, including stratified testing.
- A 'Nessie' function that outputs the estimated sample size by yearly intervals and the average lifespan expectancy left for a given group. 
- A compact, readable and efficient codebase (up to 100x less LOC than `relsurv` for the same functionalities), ensuring long-term maintenability.
- Significant performance improvements (see below) compared `relsurv`.

# Getting Started

The package is available on Julia's general registry, and thus can be installed through the following command:

```julia
using NetSurvival
```

See the rest of this [documentation](https://juliasurv.github.io/NetSurvival.jl/dev/) to have a glimpse of the functionalities! You can also take a look at [our talk at JuliaCon2024](https://www.youtube.com/watch?v=Bh3K1FHSW3A&t=2s).

# Benchmarks

`NetSurvival.jl` is *fast*. Below numbers gives runtime mulitpliers w.r.t. [`R::relsurv`](https://cran.r-project.org/web/packages/relsurv/index.html), computed on a i9-13900 processor. A version of these numbers computed on (even slower) github action's runners are availiable in [our documentation](https://juliasurv.github.io/NetSurvival.jl/dev/benches/), alongside the code needed to re-ran these numbers on your environnement. 

The comparison is done on the `colrec` dataset with the `slopop` ratetable. The first numbers compare the timing in the obtention of the net survival curve: 

|  | **Unstratified**<br>`Surv(time,status)~1` | **Stratified**<br>`Surv(time,status)~sex` |
|--------------------------:|------------------------------:|----------------------------:|
| Pohar Perme               | 20.8431                       | 20.1461                     |
| EdererI                   | 7.216                         | 4.1363                      |
| EdererII                  | 29.2397                       | 29.0399                     |
| Hakulinen                 | 23.493                        | 15.6676                     |

While the second numbers compare the implementation the Grafféo's log-rank-type test:

|  | **Unstratified**<br>`Surv(time,status)~stage` | **Stratified**<br>`Surv(time,status)~stage+Strata(sex)` |
|--------------------------:|------------------------------:|----------------------------:|
| Graffeo's LRT             | 13.1556                       | 18.156                      |

*Call to contributions* : If you have access to stata's implementation (which is not free) and want to report timings, do not hesitate to open an issue.


# Contributions are welcome

If you want to contribute to the package, ask a question, found a bug or simply want to chat, do not hesitate to open an issue on this repo. General guidelines on collaborative practices (colprac) are available at https://github.com/SciML/ColPrac.


