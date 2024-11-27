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
        site = Symbol.(colrec.site),
    )
end