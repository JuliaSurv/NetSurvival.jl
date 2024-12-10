"""
    colrec

The `colrec` dataset provides an example of relative survival dataset from slovenia, and should be paired with the `slopop` ratetable in the `RateTables.jl` package. It contains the following columns: 

- `time`, in days: follow-up time since diagnosis,
- `status`, boolean: observed death (1) or censored (0),
- `age`, in days: age at diagnosis,
- `year`, in days since 0000-00-00: date of diagnosis,
- `sex`, `:male` or `:female`,
- `stage`, cancer stage, can be 1,2,3 or 99 (unknown),
- `site`, cancer site, can be `:rectum` or `:colon`.

The 5971 patients were diagnosed with colon or rectal cancer in 1994-2000. The original source of the dataset is the [`relsurv` R package](https://CRAN.R-project.org/package=relsurv). Data were provided by the Slovene Cancer Registry, and the `age`,`time` and `year` variables were randomly perturbed to make the identification of patients impossible. Due to this perturbation, no medical conclusion should be made from results extracted from this dataset.

References: 
* [Pavlik2018](@cite) Pohar Perme, Maja  and Pavlic, Klemen (2018). Nonparametric relative survival analysis with the R package relsurv. Journal of Statistical Software
* [Zadnik2012](@cite) Zadnik V, Primic Žakelj M, Krajc M (2012). Cancer Burden in Slovenia in Comparison with the Burden in Other European Countries. Zdravniški Vestnik
* [Zadnik2016](@cite) Zadnik V, Žagar T, Primic Žakelj M (2016). Cancer Patients’ Survival: Standard Calculation Methods and Some Considerations Regarding Their Interpretation. Zdravstveno Varstvo
"""
colrec = let 
    colrec = CSV.read(
        joinpath(@__DIR__,"..","data","colrec.csv"), 
        DataFrames.DataFrame
    )
    DataFrame(
        time = colrec.time,
        status = colrec.stat.==1,
        age = colrec.age,
        year = trunc.((1960*365.241) .+ colrec.diag),
        sex = ifelse.(colrec.sex .== 1, :male, :female),
        stage = colrec.stage,
        site = Symbol.(colrec.site),
    )
end

"""
    ccolon

The `ccolon` dataset provides an example of relative survival dataset from france, and should be paired with the `survexp_fr` ratetable in the `RateTables.jl` package. It contains the following columns: 

- `time`, in days: follow-up time since diagnosis,
- `status`, boolean: observed death (1) or censored (0),
- `age`, in days: age at diagnosis,
- `year`, in days since 0000-00-00: date of diagnosis,
- `sex`, `:male` or `:female`,
- `stage`, interger from 0 to 3, correspnding to Cancer TNM (tumor node metastatis) stages at diagnosis, either I, II, III, IIIb or IV. 
- `side`, the primary tumor location, either `:right` or `:left` of the colon.

This dataset has been studied in [Giorgi2003](@cite), [Wolski2020](@cite) and [Laverny2024](@cite). It contains population-based survival data on cases of colorectal cancer from the Registry of Digestive Cancers in Burgundy, France, diagnosed between 1976 and 1990. The `year` variable was randomly perturbed, and patients were censored 1) after 10 years of followup and 2) on the 31 December 1994 based on this perturbed year information, for anonymity purposes. Hence, no medical conclusion should be made from results extracted from this dataset.

References: 
* [Giorgi2003](@cite) Giorgi, Roch and Abrahamowicz, Michal and Quantin, Catherine and Bolard, Philippe and Esteve, Jacques and Gouvernet, Joanny and Faivre, Jean. A relative survival regression model using B-spline functions to model non-proportional hazards. Statistics in medicine
* [Wolski2020](@cite) Wolski, Anna and Graff{\'e}o, Nathalie and Giorgi, Roch and {the CENSUR working survival group}. A Permutation Test Based on the Restricted Mean Survival Time for Comparison of Net Survival Distributions in Non-Proportional Excess Hazard Settings. Statistical Methods in Medical Research.
* [Laverny2024](@cite) Laverny, Oskar and Graff{\'e}o, Nathalie and Giorgi, Roch. Non-parametric estimation of net survival under dependence between death causes. To appear. 
"""
ccolon = let 
    ccolon = CSV.read(
        joinpath(@__DIR__,"..","data","ccolon.csv"), 
        DataFrames.DataFrame
    )
    DataFrame(
        time = ccolon.time,
        status = ccolon.status.==1,
        age = ccolon.age,
        year = 365.241 * ccolon.year,
        sex = Symbol.(ccolon.sex),
        stage = ccolon.stage,
        side = ifelse.(ccolon.rightside .== 1.0, :right, :left),
    )
end