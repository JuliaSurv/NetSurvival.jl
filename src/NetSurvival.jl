module NetSurvival
using PrecompileTools: @setup_workload, @compile_workload    # this is a small dependency
using DataFrames
using Distributions
using LinearAlgebra
using StatsAPI
using StatsBase
using StatsModels
using Tables
using Base.Cartesian
using CSV
using RateTables

include("fetch_datasets.jl")
include("Surv_and_Strata.jl")

include("NPNSEstimator.jl")
include("PoharPerme.jl")
include("EdererI.jl")
include("EdererII.jl")
include("Hakulinen.jl")

include("CrudeMortality.jl")

include("GraffeoTest.jl")

export PoharPerme, EdererI, EdererII, Hakulinen
export CrudeMortality
export fit, confint
export GraffeoTest
export Surv, Strata
export colrec
export @formula

@setup_workload begin
    # Putting some things in `@setup_workload` instead of `@compile_workload` can reduce the size of the
    # precompile file and potentially make loading faster.
    @compile_workload begin
        fit(PoharPerme, @formula(Surv(time,status)~1), colrec, slopop)
        fit(PoharPerme, @formula(Surv(time,status)~sex), colrec, slopop)
        fit(GraffeoTest, @formula(Surv(time,status)~sex), colrec, slopop)
        fit(GraffeoTest, @formula(Surv(time,status)~sex+Strata(stage)), colrec, slopop)
        fit(CrudeMortality, @formula(Surv(time,status)~1), colrec, slopop)
    end
end


end