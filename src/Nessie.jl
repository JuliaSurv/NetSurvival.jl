struct Nessie
    expected_sample_size::Vector{Float64}
    expected_life_time::Float64
    grid::Vector{Float64}
    function Nessie(T, Δ, age, year, rate_preds, ratetable)
        grid = mk_grid([1,maximum(T)],1) # mk_grid(T,1)
        exp_spl_size = zeros(length(grid))
        life_time = 0.0
        for i in eachindex(age)
            Λₚ = 0.0
            rtᵢ = ratetable[rate_preds[i,:]...]
            for j in 1:(length(grid)-1)
                Sₚ  = exp(-Λₚ)
                λₚ  = daily_hazard(rtᵢ, age[i] + grid[j], year[i] + grid[j])
                ∂Λₚ = λₚ * (grid[j+1]-grid[j]) # λₚ * ∂t 
                Λₚ += ∂Λₚ
                exp_spl_size[j] += Sₚ
                life_time += Sₚ * (1 - exp(-∂Λₚ)) / λₚ
            end
        end

        exp_life_time = life_time / 365.241 / length(age)
        annually = [searchsortedlast(grid, i) for i in (365.241 * (0:floor(maximum(T)/365.241))).+1]
        
        return new(exp_spl_size[annually], exp_life_time, grid[annually])
    end
end

"""
    nessie(formula, data, ratetable)

bla bla

"""
function nessie(args...)
    r = fit(Nessie,args...)
    transform!(r, :estimator => ByRow(x-> (x.grid, x.expected_life_time, x.expected_sample_size)) => [:grid, :expected_life_time,:expected_sample_size])
    select!(r, Not(:estimator))

    lt = deepcopy(r)
    select!(lt, Not([:expected_sample_size, :grid]))

    select!(r, Not(:expected_life_time))
    return lt, r
end


expected_life_time(x::Nessie) = x.expected_life_time
expected_sample_size(x::Nessie) = x.expected_sample_size




# function old_Nessie(formula::FormulaTerm, df::DataFrame, rt::RateTables.AbstractRateTable)
#     formula_applied = apply_schema(formula,schema(df))
#     rate_predictors = String.([RateTables.predictors(rt)...])

#     nms = StatsModels.termnames(formula_applied.rhs)
#     if isa(nms, String)
#         pred_names = [nms]
#     else
#         pred_names = nms
#     end

#     times = sort(unique(floor.(df.time ./ 365.241)))
#     times = unique([0.0; times])

#     times_d = times .* 365.241

#     new_df = groupby(df, pred_names)
#     povp = zeros(nrow(unique(df[!,pred_names])))
#     sit = zeros(length(times))
#     num_pop = zeros(nrow(unique(df[!,pred_names])), length(times))

#     for i in 1:nrow(unique(df[!,pred_names]))
#         for j in 1:nrow(new_df[i])
#             Tᵢ = searchsortedlast(times_d, new_df[i].time[j])
#             rate_preds = select(new_df[i],rate_predictors)
#             rtᵢ = rt[rate_preds[j,:]...]
#             Λₚ = 0.0

#             for m in 1:Tᵢ
#                 λₚ          = daily_hazard(rtᵢ, new_df[i].age[j] + times_d[m], new_df[i].year[j] + times_d[m])
#                 ∂Λₚ         = λₚ * 365.241
#                 Λₚ         += ∂Λₚ
#                 Sₚ          = exp(-Λₚ)
#                 num_pop[i,m] += Sₚ 
#                 sit[m]     += (1-Sₚ) / λₚ
#             end
#         end  
#         povp[i] = mean(sit ./ 365.241)     
#     end
#     return num_pop, povp
# end