
"""
fetch(v,x)

Equivalent to searchsortedlast(v,x) but with a dumb algorithm for NTuples, and a good one for FastIntegerFetchingVector
"""
fetch(v,x) = searchsortedlast(v,x)


struct FastIntegerFetchingVector
    vec::Vector{Int64}
    idx::Vector{Int64}
    m::Int64
    function FastIntegerFetchingVector(vec)
        @assert eltype(vec)<:Integer
        @assert issorted(vec)
        m, M = vec[1]-1, vec[end] # minimum and maximum integer values. 
        idx = zeros(Int, M-m)
        s = 1
        for i in 1:length(idx)
            if m + i > vec[s+1]
                s+=1
            end 
            idx[i]=s
        end
        return new(vec,idx,m)
    end
end
function fetch(v::FastIntegerFetchingVector, x)
    i = clamp(Int(trunc(x-v.m)),1,length(v.idx))
    return v.idx[i]
end

# import Base.getindex, Base.lastindex # needed in the package. 
Base.lastindex(v::FastIntegerFetchingVector) = Base.lastindex(v.vec)
Base.getindex(v::FastIntegerFetchingVector, i) = Base.getindex(v.vec,i)

@generated function fetch(v::NTuple{N,T}, x) where {N,T}
    quote
        @nexprs $(N-1) (i-> x == v[i] && return i)
        return $N
    end
end

function fast_fetching_vector(v)
    if length(v) < 10
        N = length(v)
        T = eltype(v)
        return NTuple{N,T}(v)
    elseif eltype(v) == Int64
        return FastIntegerFetchingVector(v)
    else
        return v
    end
end
