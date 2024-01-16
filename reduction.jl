function norm_combination(x)
    return norm(x,1)-norm(x,Inf)
end

norm_combination([1 ,1])

function Simple_reduce_order(P::SimpleSparsePolynomialZonotope, r::Real,
    method=LazySets.GIR05();nor=2)
    @assert r ≥ 1
    
    
    c=P.c
    G=P.G 
    E=P.E  
    n,h=size(G) #h is the number of generators


    a = max(0, min(h , ceil(Int, h - n * (r - 1))))
    #Gbar = hcat(G, GI)
    if nor==2
        norms = [norm(g) for g in eachcol(G)]
    else
        norms = [norm_combination(g) for g in eachcol(G)]
    end
    th = sort(norms)[a]

    # TODO: case a = 0
    # TODO is constructing an array of booleans the most efficient way?
    K = [norms[i] ≤ th for i in 1:h] #here is a boolean array
    Kbar = .!K


    PZ = SimpleSparsePolynomialZonotope(c, G[:, K], E[:, K])
    Z = reduce_order(overapproximate(PZ, Zonotope), 1, method)
    Ebar = E[:, Kbar]
    #N = [!iszero(e) for e in eachrow(Ebar)]
    cz = Z.center
    Gz = Z.generators
    m=size(Gz)[2]
    a,b=size(Ebar)

    return SimpleSparsePolynomialZonotope(cz,
        hcat(G[:, Kbar],Gz),
        vcat(hcat(Ebar,zeros(Int64,a,m)),hcat(zeros(Int64,m,b),Matrix(1I,m,m))))

end
