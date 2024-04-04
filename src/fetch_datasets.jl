colrec = let 
    colrec = CSV.read(
        joinpath(@__DIR__,"..","data","colrec.csv"), 
        DataFrames.DataFrame
    )
    DataFrame(
        age = colrec.age,
        year = trunc.((1960*365.241) .+ colrec.diag),
        sex = ifelse.(colrec.sex .== 1, :male, :female),
        time = colrec.time,
        status = colrec.stat.==1,
        stage = colrec.stage,
        site = colrec.site,
    )
end



slopop = let
    # Load the data from R to Julia
    df = CSV.read(
        joinpath(@__DIR__,"..","data","slopop.csv"), 
        DataFrames.DataFrame
    )

    df.age = Int.(trunc.(df.age .* 365.241))
    df.year = Int.(trunc.(df.year .* 365.241))
    df.sex = Symbol.(df.sex)

    RateTableV2(df, description="This ratetable corresponds to the `relsurv::slopop` ratetable from R. It is based on data collected on the population of Slovenia. Both age and year categories are in days.")
end



frpop = let

    df = CSV.read(
        joinpath(@__DIR__,"..","data","frpop.csv"), 
        DataFrames.DataFrame
    )

    df.age = Int.(trunc.(df.age .* 365.241))
    df.year = Int.(trunc.(df.year .* 365.241))
    df.sex = Symbol.(df.sex)

    RateTableV2(df, description="This ratetable corresponds to the `frpop.r` structure defined in R. It is based on data collected on the population of France. Both age and year categories are in days.")
end