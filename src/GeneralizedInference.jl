# Generalized estimator and test: implicit tangents and propagated covariance.
# Existing fit/variance/statistic/pvalue interfaces and result types are preserved.

_step_workspace(C, n) = nothing
_refresh_step!(C, cache, population, increment) = nothing

function _step_workspace(::ClaytonCopula, n)
    return zeros(n), zeros(n)
end
function _refresh_step!(C::ClaytonCopula, cache, population, increment)
    power, drift = cache
    theta = C.G.θ
    for i in eachindex(population)
        power[i] = population[i]^(-theta)
        drift[i] = increment[i] / population[i] * power[i]
    end
end

# Path-only recurrence; standalone step evaluations still derive powers from S_P.
_initialize_path_step!(C, cache) = nothing
_initialize_path_step!(::ClaytonCopula, cache) = fill!(first(cache), 1.0)
_refresh_path_step!(C, cache, population, increment) = _refresh_step!(C, cache, population, increment)
_refresh_path_step!(::ClaytonCopula, cache, population, increment) = nothing

function _advance_population!(C, cache, population, i, position, k, hazard_increment)
    population[i] *= exp(-hazard_increment)
    return population[i], -hazard_increment * population[i]
end
function _advance_population!(C::ClaytonCopula, cache, population, i, position, k, hazard_increment)
    power, drift = cache
    # Read the old active position before compacting into k <= position.
    power[k] = power[position] * exp(C.G.θ * hazard_increment)
    drift[k] = -hazard_increment * power[k]
    population[i] *= exp(-hazard_increment)
    population[i] > 0 || throw(DomainError(population[i], "Nonregular population survival"))
    return population[i], -hazard_increment * population[i]
end

# Hoist population-only powers, and compute F and F' in the same scalar pass.
_evaluate_step(C::ClaytonCopula, z, population, increment, deaths, cache) =
    _evaluate_step(C, z, population, increment, deaths, cache, Val(true))

function _evaluate_step(C::ClaytonCopula, z, population, increment, deaths, cache, ::Val{full}) where full
    theta = C.G.θ
    u, excess_power = exp(-z), exp(theta * z)
    power, drift = cache
    S, N, Sderivative, Nderivative, noise, observed = ntuple(_ -> zero(z), 6)
    for i in eachindex(population)
        base = excess_power + power[i] - 1
        base > 0 || throw(DomainError(base, "Nonregular copula section; estimation stops here."))
        inverse_q = base^(1 / theta)
        (isfinite(inverse_q) && inverse_q > 0 && u > 0) ||
            throw(DomainError(base, "Nonregular copula section; estimation stops here."))
        contribution = deaths[i] * base + drift[i]
        S += inverse_q
        N += inverse_q * contribution
        if full
            Sderivative += inverse_q / base
            Nderivative += inverse_q * (deaths[i] * (1 + theta) + drift[i] / base)
        end
        if full && !iszero(deaths[i])
            weight = base * inverse_q
            noise += deaths[i] * weight^2
            observed += deaths[i] * weight
        end
    end
    S > 0 || throw(DomainError(S, "Empty risk set"))
    scale = excess_power * S
    F = N / scale
    derivative = full ? Nderivative / S - F * (theta + excess_power * Sderivative / S) : zero(F)
    return F, u * S, noise / scale^2, observed / scale, derivative
end

_step_workspace(::FrankCopula, n) = zeros(n)
function _refresh_step!(C::FrankCopula, cache, population, increment)
    for i in eachindex(population)
        cache[i] = expm1(-C.G.θ * population[i])
    end
end

_evaluate_step(C::FrankCopula, z, population, increment, deaths, cache) =
    _evaluate_step(C, z, population, increment, deaths, cache, Val(true))

function _evaluate_step(C::FrankCopula, z, population, increment, deaths, cache, ::Val{full}) where full
    theta = C.G.θ
    u = exp(-z)
    eu, etheta = expm1(-theta * u), expm1(-theta)
    numerator, denominator, noise, observed, dnumerator, ddenominator = ntuple(_ -> zero(z), 6)
    for i in eachindex(population)
        ev = cache[i]
        divisor = etheta + eu * ev
        q = -log1p(eu * ev / etheta) / theta
        a = ev * (eu + 1) / divisor
        c2 = eu * (ev + 1) / divisor
        (isfinite(q) && isfinite(a) && q > 0 && a > 0) ||
            throw(DomainError((u, population[i]), "Nonregular copula section; estimation stops here."))
        numerator += deaths[i] / a + increment[i] * c2 / (a * q)
        denominator += u / q
        if full
            da = theta * u * a * (1 - a)
            dq = -u * a
            dc2 = theta * u * (eu + 1) * (ev + 1) * etheta / divisor^2
            dnumerator += -deaths[i] * da / a^2 + increment[i] / (a * q) * (dc2 - c2 * (da / a + dq / q))
            ddenominator += -u / q - u * dq / q^2
            noise += deaths[i] / a^2
            observed += deaths[i] / a
        end
    end
    denominator > 0 || throw(DomainError(denominator, "Empty risk set"))
    F = numerator / denominator
    derivative = full ? (dnumerator - F * ddenominator) / denominator : zero(F)
    return F, denominator, noise / denominator^2, observed / denominator, derivative
end

# Other copulas retain automatic differentiation, but share its value pass.
function _evaluate_step(C, z, population, population_increment, deaths, cache)
    # The copula sections and their derivatives are evaluated by the upstream
    # analytic specializations (Clayton/Frank) or automatic differentiation.
    function terms(z)
        u = exp(-z)
        inner_numerator, inner_denominator, inner_noise, inner_observed = zero(z), zero(z), zero(z), zero(z)
        for i in eachindex(population)
            q, a, c2 = _𝒞ₐₗₗ(C, u, population[i])
            (isfinite(q) && isfinite(a) && q > 0 && a > 0) ||
                throw(DomainError((u, population[i]), "Nonregular copula section; estimation stops here."))
            inner_numerator += deaths[i] / a + population_increment[i] * c2 / (a * q)
            inner_observed += deaths[i] / a
            inner_denominator += u / q
            inner_noise += deaths[i] / a^2
        end
        inner_denominator > 0 || throw(DomainError(inner_denominator, "Empty risk set"))
        return inner_numerator / inner_denominator, inner_denominator, inner_noise / inner_denominator^2, inner_observed / inner_denominator
    end
    tag = typeof(ForwardDiff.Tag(_evaluate_step, typeof(z)))
    value = terms(ForwardDiff.Dual{tag}(z, one(z)))
    return (ForwardDiff.value.(value)..., ForwardDiff.partials(value[1])[1])
end

# Value-only evaluations preserve the same summation order and section checks.
_step_value(C, z, population, increment, deaths, cache) =
    first(_evaluate_step(C, z, population, increment, deaths, cache))
_step_value(C::Union{ClaytonCopula,FrankCopula}, z, population, increment, deaths, cache) =
    first(_evaluate_step(C, z, population, increment, deaths, cache, Val(false)))

function _inference_step(C, previous, population, increment, deaths)
    cache = _step_workspace(C, length(population))
    _refresh_step!(C, cache, population, increment)
    return _inference_step(C, previous, population, increment, deaths, cache)
end

function _inference_step(C, previous, population, population_increment, deaths, cache)
    evaluate(z) = _evaluate_step(C, z, population, population_increment, deaths, cache)
    # Follow the same local branch, tolerance and line search as the reference.
    x = _step_value(C, previous, population, population_increment, deaths, cache)
    evaluated = evaluate(previous + x)
    for iteration in 1:40
        value, denominator, noise, observed, derivative = evaluated
        slope = 1 - derivative
        (isfinite(slope) && slope > 1e-8) ||
            throw(DomainError(slope, "Implicit update has no stable local branch"))
        residual = x - value
        if abs(residual) <= 1e-10 * (1 + abs(x))
            return x, denominator, noise, inv(slope), observed
        end
        change = residual / slope
        accepted = false
        for power in 0:20
            candidate = x - change * 2.0^(-power)
            try
                # Keep the first full trial: it is reused if accepted. Later
                # backtracking trials only need the residual until acceptance.
                candidate_evaluated = (power == 0 || !(C isa Union{ClaytonCopula,FrankCopula})) ?
                    evaluate(previous + candidate) : nothing
                candidate_value = candidate_evaluated === nothing ?
                    _step_value(C, previous + candidate, population, population_increment, deaths, cache) : first(candidate_evaluated)
                candidate_residual = candidate - candidate_value
                if abs(candidate_residual) < abs(residual)
                    x = candidate
                    # Reuse the accepted evaluation at the next Newton iteration.
                    evaluated = candidate_evaluated === nothing ? evaluate(previous + candidate) : candidate_evaluated
                    accepted = true
                    break
                end
            catch exception
                exception isa DomainError || rethrow()
            end
        end
        accepted || error("Implicit update failed to reduce its residual")
    end
    error("Implicit update did not converge in 40 iterations")
end

function _generalized_path(m, T, Δ, age, year, rate_preds, ratetable, grid)
    output = ntuple(_ -> zeros(length(grid)), 5)
    _generalized_path!(output..., m, T, Δ, age, year, rate_preds, ratetable, grid)
    return output
end

function _generalized_path!(increments, denominator, noise, propagation, observed, m, T, Δ, age, year, rate_preds, ratetable, grid)
    n = length(T)
    for output in (increments, denominator, noise, propagation)
        fill!(output, 0.0)
    end
    observed === nothing || fill!(observed, 0.0)
    population = ones(n)
    active_population, population_increment, deaths = zeros(n), zeros(n), zeros(n)
    cache = _step_workspace(m.C, n)
    _initialize_path_step!(m.C, cache)
    # Compact in input order, without allocating a fresh risk index each step.
    active = collect(eachindex(T))
    remaining = n
    rates = [ratetable[rate_preds[i, :]...] for i in 1:n]
    cumulative = 0.0
    previous_time = 0.0
    for j in eachindex(grid)
        time = grid[j]
        step = time - previous_time
        # After the last at-risk subject in a group leaves, stop that path.
        if remaining == 0
            propagation[j] = 1.0
            previous_time = time
            continue
        end
        # Preserve input order for numerically stable comparison to the reference.
        k = 0
        for position in 1:remaining
            i = active[position]
            T[i] >= time || continue
            k += 1
            active[k] = i
            hazard = daily_hazard(rates[i], age[i] + previous_time, year[i] + previous_time)
            hazard >= 0 || throw(DomainError(hazard, "Negative population hazard"))
            # Backward time increment [t_{j-1},t_j], not the next interval.
            active_population[k], population_increment[k] =
                _advance_population!(m.C, cache, population, i, position, k, step * hazard)
            deaths[k] = T[i] == time && Δ[i]
        end
        remaining = k
        if remaining == 0
            propagation[j] = 1.0
            previous_time = time
            continue
        end
        pop = @view active_population[1:remaining]
        dpop = @view population_increment[1:remaining]
        dn = @view deaths[1:remaining]
        _refresh_path_step!(m.C, cache, pop, dpop)
        x, D, Q, P, observation_increment = _inference_step(m.C, cumulative, pop, dpop, dn, cache)
        increments[j], denominator[j], noise[j], propagation[j] = x, D, Q, P
        observed === nothing || (observed[j] = observation_increment)
        cumulative += x
        previous_time = time
    end
    return nothing
end

# More specific than the upstream Λ(::NPNSMethod,...) method: no overwrite.
function Λ(m::GenPoharPermeMethod, T, Δ, age, year, rate_preds, ratetable, grid)
    x, _, Q, P, observed = _generalized_path(m, T, Δ, age, year, rate_preds, ratetable, grid)
    variance_increments = similar(x)
    variance = 0.0
    for j in eachindex(x)
        next_variance = P[j]^2 * (variance + Q[j])
        variance_increments[j] = next_variance - variance
        variance = next_variance
    end
    # Keep the field decomposition but use the actual cumulative hazard increment.
    return observed, observed .- x, variance_increments
end

function _propagate_joint_covariance(covariance, D, P, Q)
    groups = length(D)
    result = similar(covariance)
    R, B, left = zeros(groups), zeros(groups, groups), similar(covariance)
    _propagate_joint_covariance!(result, covariance, D, P, Q, R, B, left)
    return result, R
end

# Apply the block transition without allocating matrices or calling tiny GEMMs.
function _propagate_joint_covariance!(result, covariance, D, P, Q, R, B, left)
    groups = length(D)
    total = sum(D)
    for i in 1:groups
        R[i] = total > 0 ? D[i] / total : 0.0
        for ell in 1:groups
            B[i, ell] = ((i == ell) - R[i]) * D[ell]
        end
    end
    for j in 1:2groups, i in 1:groups
        left[i, j] = P[i] * covariance[i, j]
        value = covariance[groups+i, j]
        for ell in 1:groups
            value += B[i, ell] * (P[ell] - 1) * covariance[ell, j]
        end
        left[groups+i, j] = value
    end
    for j in 1:groups, i in 1:2groups
        result[i, j] = left[i, j] * P[j]
        value = left[i, groups+j]
        for ell in 1:groups
            value += left[i, ell] * B[j, ell] * (P[ell] - 1)
        end
        result[i, groups+j] = value
    end
    for j in 1:groups, i in 1:groups
        noise = P[i]^2 * Q[i]
        i == j && (result[i, j] += noise)
        result[i, groups+j] += noise * B[j, i]
        result[groups+j, i] += noise * B[j, i]
        value = 0.0
        for ell in 1:groups
            value += B[i, ell] * P[ell]^2 * Q[ell] * B[j, ell]
        end
        result[groups+i, groups+j] += value
    end
    for j in 1:2groups, i in 1:j-1
        value = (result[i, j] + result[j, i]) / 2
        result[i, j] = result[j, i] = value
    end
    return result
end

function _generalized_graffeo(method, T, Δ, age, year, rate_preds, strata, group, ratetable)
    grid = mk_grid(T, 1)
    groups, stratas = unique(group), unique(strata)
    k, jmax = length(groups), length(grid)
    k >= 2 || throw(ArgumentError("The comparison needs at least two groups"))
    ∂N, ∂V, ∂Z, Dall, Rall = (zeros(k, jmax) for _ in 1:5)
    ∂VZ = zeros(k, k, jmax)
    xs, Ds, Qs, Ps = (zeros(k, jmax) for _ in 1:4)
    covariance = zeros(2k, 2k)
    next_covariance, left = similar(covariance), similar(covariance)
    R, B = zeros(k), zeros(k, k)
    for stratum in stratas
        for output in (xs, Ds, Qs, Ps)
            fill!(output, 0.0)
        end
        fill!(covariance, 0.0)
        for g in eachindex(groups)
            idx = (strata .== stratum) .& (group .== groups[g])
            if !any(idx)
                Ps[g, :] .= 1
                continue
            end
            _generalized_path!(view(xs, g, :), view(Ds, g, :), view(Qs, g, :), view(Ps, g, :), nothing, method, T[idx], Δ[idx], age[idx], year[idx],
                                         rate_preds[idx, :], ratetable, grid)
            for j in eachindex(grid)
                ∂N[g, j] += Ds[g, j] * xs[g, j]
                ∂V[g, j] += Ds[g, j]^2 * Qs[g, j]
            end
        end
        for j in eachindex(grid)
            @views _propagate_joint_covariance!(next_covariance, covariance, Ds[:, j], Ps[:, j], Qs[:, j], R, B, left)
            score_total = 0.0
            for g in 1:k
                score_total += Ds[g, j] * xs[g, j]
            end
            for g in 1:k
                ∂Z[g, j] += Ds[g, j] * xs[g, j] - R[g] * score_total
                Dall[g, j] += Ds[g, j]
                for h in 1:k
                    ∂VZ[g, h, j] += next_covariance[k+g, k+h] - covariance[k+g, k+h]
                end
            end
            covariance, next_covariance = next_covariance, covariance
        end
    end
    for j in eachindex(grid)
        total = sum(@view Dall[:, j])
        for g in 1:k
            Rall[g, j] = total > 0 ? Dall[g, j] / total : 0.0
        end
    end
    Z = vec(sum(∂Z; dims=2))
    covariance = dropdims(sum(∂VZ; dims=3); dims=3)
    reduced = Symmetric(covariance[1:end-1, 1:end-1])
    # Rank failure is an explicit failure, not a silently regularized p-value.
    statistic = dot(Z[1:end-1], cholesky(reduced) \ Z[1:end-1])
    return grid, ∂N, ∂V, ∂Z, Dall, Rall, ∂VZ, statistic, k - 1, ccdf(Chisq(k - 1), statistic)
end
