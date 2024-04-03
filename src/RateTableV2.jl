

struct RateTableV2{n_axes,axes_type,values_type}
    axes::axes_type
    values::values_type
    description::String
    function RateTableV2(df;description="")
        nms = Symbol.(names(select(df,Not(:value))))
        axes = NamedTuple(n => sort(unique(df[!,n])) for n in nms)
        axes_lengths = Tuple(length(axes[i]) for i in eachindex(axes))
        values = zeros(axes_lengths)
        for i in 1:length(df.value)
            ci = CartesianIndex(Tuple(searchsortedlast(axes[n], df[i,n]) for n in nms))
            values[ci] = df.value[i]
        end
        axes = NamedTuple(n =>fast_fetching_vector(axes[n]) for n in nms)
        return new{length(axes),typeof(axes),typeof(values)}(axes, values,description)
    end
end

@generated function λ(rt::RateTableV2{N,A,V}, args...) where {N,A,V}
    quote
        Vals = rt.values
        @nexprs $N j->(i_j = fetch(rt.axes[j],args[j]))
        @nref $N Vals i
    end
end
λ(rt::RateTableV2; kwargs...) = λ(rt, (kwargs[k] for k in keys(rt.axes))...)

function Base.show(io::IO, rt::RateTableV2)
    println(io, "$(typeof(rt)) with $(length(rt.axes)) named axes:")
    for n in keys(rt.axes)
        println(io, "- `$(n)` from $(rt.axes[n][1]) to $(rt.axes[n][end]) (This axe is a $(typeof(rt.axes[n])))")
    end
    println(io, rt.description)
end




