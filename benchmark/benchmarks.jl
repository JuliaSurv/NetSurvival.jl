using NetSurvival
using Ratetables
using BenchmarkTools


# example here : https://github.com/JuliaCI/PkgBenchmark.jl/blob/master/benchmark/benchmarks.jl

suite = BenchmarkGroup()
# SUITE["rand"] = @benchmarkable rand(10)


suite["PoharPerme"] = BenchmarkGroup()
suite["PoharPerme"]["colrec x frpop - formula"] = @benchmarkable fit(PoharPerme, @formula(Surv(time,status)~sex), colrec, frpop)
suite["PoharPerme"]["colrec x slopop - explicit"] = @benchmarkable PoharPerme(colrec.time, colrec.status, colrec.age, colrec.year, colrec.sex, slopop)

# Next, the benchmark with hmd rates: 
# colrec.country = :fracnp
# SUITE["PoharPerme"]["colrec x hmd_rates"] = @benchmarkable fit(PoharPerme, @formula(Surv(time,status)~sex), colrec, hmd_rates)

# Then, a benchmark for GraffeoTest: 
suite["GraffeoTest"] = BenchmarkGroup()
suite["GraffeoTest"]["colrec x slopop - explicit"] = @benchmarkable GraffeoTest(colrec.time, colrec.status, colrec.age, colrec.year, colrec.sex, colrec.sex, colrec.stage, slopop)

tune!(suite)
results = run(suite, verbose = true)

BenchmarkTools.save("output.json", median(results))