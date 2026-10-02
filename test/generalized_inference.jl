using Test, Random, LinearAlgebra
using NetSurvival, Copulas, RateTables, DataFrames, Distributions



@testset "Fused derivatives versus automatic differentiation" begin
    for C in (ClaytonCopula(2, -0.4), ClaytonCopula(2, 0.1), ClaytonCopula(2, 6/7),
              ClaytonCopula(2, 2.0), FrankCopula(2, -8.0), FrankCopula(2, -0.2),
              FrankCopula(2, 0.2), FrankCopula(2, 2.0), FrankCopula(2, 8.0)),
        z in (-0.0001, 0.0, 0.4, 1.2), deaths in ([0.0, 0.0], [0.0, 1.0])
        population, increment = [0.99, 0.9], [-0.0001, -0.0002]
        cache = NetSurvival._step_workspace(C, 2)
        NetSurvival._refresh_step!(C, cache, population, increment)
        analytic = NetSurvival._evaluate_step(C, z, population, increment, deaths, cache)
        @test NetSurvival._step_value(C, z, population, increment, deaths, cache) == first(analytic)
        automatic = invoke(NetSurvival._evaluate_step, Tuple{Any, Any, Any, Any, Any, Any},
                           C, z, population, increment, deaths, nothing)
        @test all(isapprox.(analytic, automatic; rtol=1e-9, atol=1e-11))
    end
end

@testset "Clayton counterexample: corrected continuous-time limits" begin
    theta, excess, population_rate, horizon = 6 / 7, 0.1, 0.1, 10.0
    C = ClaytonCopula(2, theta)
    function limits(steps)
        dt = horizon / steps
        variance = 0.0
        covariance = zeros(4, 4)
        old_variance, old_score_variance = 0.0, 0.0
        for j in 1:steps
            t = j * dt
            u, v = exp(-excess * t), exp(-population_rate * t)
            base = u^(-theta) + v^(-theta) - 1
            q, r = base^(-1 / theta), u^(-theta) / base
            hazard = (excess * u^(-theta) + population_rate * v^(-theta)) / base
            # Expected observed increments isolate discretization and tangents;
            # the number at risk is n*q, hence the variance normalization by q.
            _, _, raw, P = NetSurvival._inference_step(C, excess * (t - dt),
                [v], [-population_rate * v * dt], [hazard * dt])
            Q = raw / q
            variance = P^2 * (variance + Q)
            covariance, _ = NetSurvival._propagate_joint_covariance(covariance,
                fill(u / 2, 2), fill(P, 2), fill(2Q, 2))
            old_variance += hazard / (q * r^2) * dt
            old_score_variance += u^2 * hazard / (q * r^2) * dt / 4
        end
        return variance, covariance[3, 3], old_variance, old_score_variance
    end
    coarse, fine = limits(1000), limits(2000)
    u, v = exp(-excess * horizon), exp(-population_rate * horizon)
    base = u^(-theta) + v^(-theta) - 1
    q, r = base^(-1 / theta), u^(-theta) / base
    exact_variance = (1 / q - 1) / r^2
    @test fine[1] ≈ exact_variance rtol=0.003
    @test abs(fine[1] - exact_variance) < abs(coarse[1] - exact_variance)
    @test fine[1] / fine[3] ≈ 1.291890 rtol=0.003
    @test fine[2] / fine[4] ≈ 1.343107 rtol=0.003
end

@testset "Implicit step and independent-copula baseline" begin
    C = IndependentCopula(2)
    population = [0.8, 0.9]
    dpopulation = [-0.008, -0.018]
    deaths = [1.0, 0.0]
    x, D, Q, P = NetSurvival._inference_step(C, 0.3, population, dpopulation, deaths)
    expected = sum(deaths ./ population .+ dpopulation ./ population.^2) / sum(1 ./ population)
    @test x ≈ expected
    @test D ≈ sum(1 ./ population)
    @test Q ≈ sum(deaths ./ population.^2) / D^2
    @test P ≈ 1.0
    @test_throws DomainError NetSurvival._inference_step(C, 0.3, Float64[], Float64[], Float64[])
    @test_throws DomainError NetSurvival._inference_step(ClaytonCopula(2, -0.5),
        5.0, [0.01], [-0.00001], [0.0])
    initial = first(NetSurvival._inference_step(ClaytonCopula(2, 1.0),
        0.0, [0.999], [-0.001], [0.0]))
    @test initial < 0 # the local extension permits survival above one
end

@testset "Derivative of the actual implicit solution" begin
    for C in (ClaytonCopula(2, 1.0), FrankCopula(2, 2.0), FrankCopula(2, -2.0),
              FGMCopula(2, 0.3))
        previous = 0.4
        population, dpopulation, deaths = [0.7, 0.8], [-0.0007, -0.0008], [0.0, 1.0]
        x, _, _, P = NetSurvival._inference_step(C, previous, population, dpopulation, deaths)
        epsilon = 1e-6
        plus = first(NetSurvival._inference_step(C, previous + epsilon, population, dpopulation, deaths))
        minus = first(NetSurvival._inference_step(C, previous - epsilon, population, dpopulation, deaths))
        @test P ≈ 1 + (plus - minus) / (2epsilon) rtol=1e-6
        @test isfinite(x)
    end
end

@testset "Covariance recursion versus explicit influence matrix" begin
    rng = MersenneTwister(42)
    k, steps = 3, 8
    covariance = zeros(2k, 2k)
    influence = zeros(2k, 0)
    for j in 1:steps
        D, P, Q = rand(rng, k) .+ 2, rand(rng, k) .* 0.2 .+ 0.9, rand(rng, k) .* 0.01
        R = D / sum(D)
        B = (Matrix{Float64}(I, k, k) .- R) * Diagonal(D)
        L = [Diagonal(P) zeros(k, k); B * Diagonal(P .- 1) Matrix{Float64}(I, k, k)]
        loading = vcat(Diagonal(P), B * Diagonal(P)) * Diagonal(sqrt.(Q))
        influence = hcat(L * influence, loading)
        covariance, _ = NetSurvival._propagate_joint_covariance(covariance, D, P, Q)
        @test covariance ≈ influence * influence'
        @test minimum(eigvals(Symmetric(covariance))) > -1e-10
        @test covariance[k+1:end, k+1:end] * ones(k) ≈ zeros(k) atol=1e-10
    end
    # P=1 recovers the original log-rank covariance exactly.
    D, Q = [10.0, 20.0], [0.01, 0.02]
    C, R = NetSurvival._propagate_joint_covariance(zeros(4, 4), D, ones(2), Q)
    B = (Matrix{Float64}(I, 2, 2) .- R) * Diagonal(D)
    @test C[3:4, 3:4] ≈ B * Diagonal(Q) * B'
end

@testset "End-to-end APIs, stopping, and horizons" begin
    rng = MersenneTwister(12)
    n = 80
    data = DataFrame(time=min.(rand(rng, Exponential(300.0), n), 500.0),
                     status=trues(n), age=fill(50 * 365.241, n),
                     year=fill(2000 * 365.241, n), sex=repeat([:male, :female], n ÷ 2),
                     group=repeat([1, 2], n ÷ 2))
    data.status = data.time .< 500
    method = GenPoharPerme(ClaytonCopula(2, 1.0))
    estimate = fit(method, @formula(Surv(time, status) ~ 1), data, slopop)
    @test all(isfinite, estimate.Sₑ)
    @test all(isfinite, estimate.σₑ)
    @test estimate.Sₑ ≈ exp.(-cumsum(estimate.∂Λₑ))
    @test variance(estimate, 200.0) >= 0
    @test variance(estimate, 0.0) == 0
    empty_tail = NetSurvival._generalized_path(method, [10.0], [false],
        [50 * 365.241], [2000 * 365.241], DataFrame(sex=[:male]), slopop,
        [5.0, 10.0, 15.0])
    @test empty_tail[1][end] == empty_tail[2][end] == empty_tail[3][end] == 0
    @test empty_tail[4][end] == 1
    test = fit(GraffeoTest(ClaytonCopula(2, 1.0)), @formula(Surv(time, status) ~ group), data, slopop)
    @test isfinite(test.stat)
    @test 0 <= pvalue(test, 200.0) <= 1
    @test statistic(test, maximum(test.t)) ≈ test.stat
    @test sum(test.∂Z; dims=1) ≈ zeros(1, length(test.t)) atol=1e-10
    independent = fit(GenPoharPerme(IndependentCopula(2)), @formula(Surv(time, status) ~ 1), data, slopop)
    pohar_perme = fit(PoharPerme, @formula(Surv(time, status) ~ 1), data, slopop)
    @test independent.Sₑ == pohar_perme.Sₑ
    @test independent.σₑ == pohar_perme.σₑ
    permutation = randperm(rng, n)
    shuffled = fit(method, @formula(Surv(time, status) ~ 1), data[permutation, :], slopop)
    @test shuffled.Sₑ ≈ estimate.Sₑ
    @test shuffled.σₑ ≈ estimate.σₑ
    three_groups = copy(data)
    three_groups.group = repeat(1:3, cld(n, 3))[1:n]
    three = fit(GraffeoTest(ClaytonCopula(2, 1.0)), @formula(Surv(time, status) ~ group), three_groups, slopop)
    @test three.df == 2
    @test 0 <= pvalue(three, 200.0) <= 1
end

@testset "Clayton population recurrence across risk-set compaction" begin
    C = ClaytonCopula(2, 0.1)
    cache = NetSurvival._step_workspace(C, 1)
    NetSurvival._initialize_path_step!(C, cache)
    @test_throws DomainError NetSurvival._advance_population!(C, cache, ones(1), 1, 1, 1, 1000.0)
    for C in (ClaytonCopula(2, -0.4), ClaytonCopula(2, 6/7), ClaytonCopula(2, 2.0))
        population = ones(3)
        cache = NetSurvival._step_workspace(C, 3)
        NetSurvival._initialize_path_step!(C, cache)
        active = [1, 2, 3]
        for step in 1:4000
            step == 1000 && deleteat!(active, 2)
            values, increments = zeros(length(active)), zeros(length(active))
            for (position, i) in enumerate(active)
                # Unequal patients and fractional intervals; no time rounding.
                h = (0.3 + mod(step, 7)/10) * i * 1e-4
                previous_position = step == 1000 && i == 3 ? 3 : position
                values[position], increments[position] = NetSurvival._advance_population!(
                    C, cache, population, i, previous_position, position, h)
            end
            if step in (1, 999, 1000, 2000, 4000)
                reference = NetSurvival._step_workspace(C, length(active))
                NetSurvival._refresh_step!(C, reference, values, increments)
                @test cache[1][1:length(active)] ≈ reference[1] rtol=1e-11
                @test cache[2][1:length(active)] ≈ reference[2] rtol=1e-11
                deaths = zeros(length(active)); deaths[end] = 1
                actual = NetSurvival._evaluate_step(C, .1, values, increments, deaths, cache)
                expected = NetSurvival._evaluate_step(C, .1, values, increments, deaths, reference)
                @test all(isapprox(a, b; rtol=1e-10, atol=1e-12) for (a, b) in zip(actual, expected))
            end
        end
    end
end
