function Nessie(formula::FormulaTerm, df::DataFrame, rt::RateTables.AbstractRateTable)
    formula_applied = apply_schema(formula,schema(df))
    rate_predictors = String.([RateTables.predictors(rt)...])

    nms = StatsModels.termnames(formula_applied.rhs)
    if isa(nms, String)
        pred_names = [nms]
    else
        pred_names = nms
    end

    times = sort(unique(floor.(df.time ./ 365.241)))
    times = unique([0.0; times])

    times_d = times .* 365.241

    new_df = groupby(df, pred_names)
    povp = zeros(nrow(unique(df[!,pred_names])))
    sit = zeros(length(times))
    num_pop = zeros(nrow(unique(df[!,pred_names])), length(times))

    for i in 1:nrow(unique(df[!,pred_names]))
        for j in 1:nrow(new_df[i])
            Tᵢ = searchsortedlast(times_d, new_df[i].time[j])
            rate_preds = select(new_df[i],rate_predictors)
            rtᵢ = rt[rate_preds[j,:]...]
            Λₚ = 0.0

            for m in 1:Tᵢ
                λₚ          = daily_hazard(rtᵢ, new_df[i].age[j] + times_d[m], new_df[i].year[j] + times_d[m])
                ∂Λₚ         = λₚ 
                Λₚ         += ∂Λₚ
                Sₚ          = exp(-Λₚ)
                num_pop[i,m] += Sₚ 
                sit[m]     += (1-Sₚ) / λₚ
            end
        end  
        povp[i] = mean(sit ./ 365.241)     
    end
    return num_pop, povp
end