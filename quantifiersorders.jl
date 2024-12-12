using Combinatorics

function mysetdiff(y, x)
    res = Vector{eltype(y)}(undef, length(y) - length(x))
    i = 1
    @inbounds for el in y
        el ∈ x && continue
        res[i] = el
        i += 1
    end
    res
end

function quant_and_varsorders(n1,n2)
    indexbase=collect(1:n1)
    nt=n1+n2
    index=collect(n1+1:nt)
    base=["forall" for j in 1:n1]
    res=[]
    for i in 1:n2-1
        q1=vcat(base,["forall" for j in 1:i],["exists" for j in 1:n2-i])
        q2=vcat(base,["forall" for j in 1:n2-i],["exists" for j in 1:i])
        @show(q1)
        comb=collect(combinations(index,i))
        @show(comb)
        diff=setdiff(index,comb)
        @show(diff)
        for c in comb
            @show(c)
            ordersq1=vcat(indexbase,c,diff)
            ordersq2=vcat(indexbase,diff,c)
            push!(res,[[q1,ordersq1],[q2,ordersq2]])
        end
    end
    return res
end

quant_and_varsorders(2,3)