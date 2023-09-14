#THE CODE IS STRONGLY INSPIRED BY NIKLAS KOCHDUMPER'S CODE FROM CORA
using LazySets

#QUADMAP SUR UN SPZ ET QUADMAP ENTRE DEUX SPZ DIFFERENTS EST DEJA DEFINIE DANS LAZYSETS

function project(PZ,dim)
    cnew=PZ.c[dim]
    Gnew=PZ.G[dim,:]
    return SimpleSparsePolynomialZonotope(cnew,Gnew,PZ.E)
end

function polynomial_map(PZ,coeffs,expmat)
    n=size(PZ.G)[1] #dimension of the polynomial zonotope
    summands=[]
    c=zeros(size(coeffs)[1])
    nbmonomials=size(PZ.E)[2]
    #tmp=Array{Int64}(undef, 1,nbmonomials)
    tmp=sum(E,dims=1)#summin E column-wise to see if there are null monomials
    ind=[]
    for i in 1:nbmonomials
        if tmp[i]==0
            push!(ind,i)
        end
    end
    if length(ind)>1
        println("REDUNDANT GENERATORS NOT REMOVED IN THE REPRESENTATION ")
    elseif length(ind)==1
        c=sum(coeffs[:,ind[1]])
        coeffs=coeffs[:,1:end .!=ind[1]]
        expmat=expmat[:,1:end .!=ind[1]]
    end

    for i in 1:size(expmat)[2]
        E=Array{Int64}(undef,0)
        for j in 1:size(expmat)[1]
            if expmat[j,i]!=0
                append!(E,j*ones(Int64,expmat[j,i]))
            end
        end
        if length(E)%2==1
            list=[project(PZ,E[end])]
        else
            list=[]
        end


    end

end

c_ = [0.0, 0]
    G = [-1.0 -2.0 -1.0 2.0 0.01 0.4
         1.0 0.0 -1.0 1.0 0.2 0.0]

    GI = [0.2 0.01
          0.02 -0.4]

    E = [1 0 1 2 2 0
         0 1 1 0 0 2
         0 0 0 0 1 2]

c_[end]
a=sum(E,dims=1)
a
vec(a==0)
a[6]
tm=[]
length(tm)
push!(tm,3)
sum(G[:,6])
G=G[:,1:end .!=6]
G
o=Array{Int64}(undef,0)
append!(o,ones(Int64,2))
o

o=3*o
append!(o,ones(Int64,2))
o
E[1,6]

    P = SparsePolynomialZonotope(c_, G, GI, E)
    size(P.G)[1]

3%2