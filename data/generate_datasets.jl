using RCall, CSV
# Load the data from R to Julia
R"""
colrec = relsurv::colrec
"""
colrec = @rget colrec
CSV.write(joinpath(@__DIR__, "colrec.csv"),colrec)