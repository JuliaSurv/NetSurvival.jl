using TestItemRunner
 # filter=ti->(:Prop in ti.tags)
include("generalized_inference.jl")

@run_package_tests
