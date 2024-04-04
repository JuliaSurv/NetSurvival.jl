using RCall, CSV
# Load the data from R to Julia
R"""
colrec = relsurv::colrec
"""
colrec = @rget colrec
CSV.write(joinpath(@__DIR__, "colrec.csv"),colrec)


R"""
slopop = relsurv::slopop
slopop1 <- expand.grid(age = as.numeric(dimnames(slopop)$age), year=as.numeric(dimnames(slopop)$year), sex = dimnames(slopop)$sex)
slopop1$value <- c(slopop[,,1],slopop[,,2])
"""
slopop = @rget slopop1
CSV.write(joinpath(@__DIR__, "slopop.csv"),slopop)


R"""
source("data/frpop.r")
frpop1 <- expand.grid(age = as.numeric(dimnames(frpop)$age), year=as.numeric(dimnames(frpop)$year), sex = dimnames(frpop)$sex)
frpop1$value <- c(frpop[,,1],frpop[,,2])
"""
frpop = @rget frpop1
CSV.write(joinpath(@__DIR__, "frpop.csv"),frpop)

