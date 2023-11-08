#THE CODE IS STRONGLY INSPIRED BY NIKLAS KOCHDUMPER'S CODE FROM CORA
using LazySets

#QUADMAP SUR UN SPZ ET QUADMAP ENTRE DEUX SPZ DIFFERENTS EST DEJA DEFINIE DANS LAZYSETS

function project(PZ,dim)
    cnew=PZ.c[dim]
    Gnew=Array{Float64}(undef,1,0)
    Gnew=PZ.G[dim,:]
    return SimpleSparsePolynomialZonotope(cnew,Gnew,PZ.E)
end

function exact_sum_SSPZ(PZ1,PZ2) #CELLE CI NE FAIT PAS TELLEMENT DE SENS
    cnew=PZ1.c+PZ2.c
    Gnew=vcat(PZ1.G,PZ2.G)
    Enew=vcat(PZ1.E,PZ2.E)
    return remove_redundant_generators(cnew,Gnew,Enew)
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

        for j in 2:2:length(E)
            Q=zeros(n,n)
            Q[E[j-1],E[j]]=1
            push!(list,quadratic_map([Q],PZ))
        end
        
        res=list[1]
        GI = Array{Float64}(undef,1,1)
        GI[1,1]=1
        for j in 2:length(list)
            res=quadratic_map([GI],res,list[j])
        end

        push!(summands,coeffs[i,:]*res)
    end

    if length(summands)>0
        res=summands[1]+c
        for j in 2:length(summands)
            res=exact_sum_SSPZ(res,summands[j])
        end
    else
        res=0*PZ+c
    end

    return res

end

c=[1 , 2.0]
G=[1.0 -2 1; 2 3 1]
E=[1 0 2;0 1 1]

Ptest=SimpleSparsePolynomialZonotope(c,G,E)
coef = [1 2;-1 0]
expMat = [0 1;3 2]

polynomial_map(Ptest,coef,expMat)


zeros(Int64,2,2)

c_ = [1.0, 0]
c_+c_
    G = [-1.0  0.01 0.4 -1
         1.0 0.2 0.0 1]

    GI = [0.2 0.01
          0.02 -0.4]

    E = [1 2 0 1
         0 0 2 0
         0 1 2 0]

Ptest=SimpleSparsePolynomialZonotope(c_,G,E)
remove_redundant_generators(Ptest)  
Ptest+Ptest

c_[end+1]=1
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