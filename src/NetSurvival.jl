module NetSurvival

using DataFrames
using Distributions
using LinearAlgebra
using StatsAPI
using StatsBase
using Tables
using CSV
using RateTables
using StatsModels
using SurvivalBase: Surv, Strata
using Copulas
using Roots
using ForwardDiff

include("fetch_datasets.jl")
include("Nessie.jl")
include("NPNSEstimator.jl")
include("NPNSEMethods/PoharPerme.jl")
include("NPNSEMethods/EdererI.jl")
include("NPNSEMethods/EdererII.jl")
include("NPNSEMethods/Hakulinen.jl")
include("NPNSEMethods/GenPoharPerme.jl")
include("CrudeMortality.jl")
include("GraffeoTest.jl")

export PoharPerme, EdererI, EdererII, Hakulinen, GenPoharPerme
export CrudeMortality
export Nessie, nessie
export fit, confint, variance, statistic, pvalue
export GraffeoTest
export Surv, Strata
export colrec, ccolon
export @formula

end