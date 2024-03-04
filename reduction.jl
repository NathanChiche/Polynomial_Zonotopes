function norm_combination(x)
    return norm(x,1)-norm(x,Inf)
end


function Simple_reduce_order(P::SimpleSparsePolynomialZonotope, r::Real,
    method=LazySets.GIR05();nor=2)
    @assert r ≥ 1
    
    c=P.c
    G=P.G 
    n,h=size(G) #h is the number of generators
    println("ICI L ORDRE: ",h//n)
    if h//n<=r 
        println("CAS PARTICULIER SANS REDUCTION")
        return P
    end
    E=P.E  

    a = max(0, min(h , ceil(Int, h - n * (r - 1))))
    if a==0
        println("A=0 !!!!!! BIG WARNING")
        return P
    end
    #Gbar = hcat(G, GI)
    if nor==2
        norms = [norm(g) for g in eachcol(G)]
    else
        norms = [norm_combination(g) for g in eachcol(G)]
    end
    th = sort(norms)[a]

    # TODO: case a = 0
    # TODO is constructing an array of booleans the most efficient way?
    K = [norms[i] ≤ th for i in 1:h] #here is a boolean array LE PROBLEME DE NON REDUCTION SE TROUVE ICI
    Kbar = .!K
    #=if length(Kbar)>r 
        println("NAN")
    end=#


    PZ = SimpleSparsePolynomialZonotope(c, G[:, K], E[:, K])
    Z = reduce_order(overapproximate(PZ, Zonotope), 1, method)
    Ebar = E[:, Kbar]
    
    #N = [!iszero(e) for e in eachrow(Ebar)]
    cz = Z.center
    Gz = Z.generators
    m=size(Gz)[2]
    a,b=size(Ebar)
    #println("on a réduit et voici la taille: ",b)

    return SimpleSparsePolynomialZonotope(cz,
        hcat(G[:, Kbar],Gz),
        vcat(hcat(Ebar,zeros(Int64,a,m)),hcat(zeros(Int64,m,b),Matrix(1I,m,m))))

end

function remove_unused_variables(SPZ)
    toremove=[]
    nbrows,nbcol=size(SPZ.E)
    zer=zeros(Int64,nbcol)
    for i in 1:nbrows
        if @view(SPZ.E[i,:])==zer
            push!(toremove,i)
        end
    end
    return SimpleSparsePolynomialZonotope(SPZ.c,SPZ.G,SPZ.E[setdiff(1:end,toremove),1:end])
end


function factor_identity(SPZ)
    toremove=[]
    nbrows,nbcol=size(SPZ.E)
    mi=min(nbrows,nbcol)
    for i in 1:mi
        if sum(@view(E[:,nbcol-i])) && sum(@view(E[nbrows-i,:]))
            push!(toremove,i)
        end
    end
end
