
function createIndexes(maxdeg)
    indx=ntuple(d->1:maxdeg[1],1)
    for i in 2:length(maxdeg)
        indx=(indx...,1:maxdeg[i])
    end
    #return indx
    res=[I for I in Iterators.product(indx...)]
    res=collect(collect(i) for i in res)
    return res
end
c
length(c)
c[12]
l::Int=3
coeffs = Vector(undef, 2)
coeffs[1]=1;coeffs[2]=1
am=[2, 1.1]
coeffs+am
lll=
function bernsteinmulti_from_univariates(monomial,maxdeg,domain,indexlist)
    ll=multivariate(monomial,maxdeg,domain)
    nbvar=length(maxdeg)
    nb_coeffs=prod(k+1 for k in maxdeg)
    coeffs=Vector{Float64}(undef,nb_coeffs)
     for i in 1:length(indexlist)
        coeffs[i]=prod(ll[j][indexlist[i][j]] for j in 1:len)
    end
end


Ma=[0 1; 1 1]
D=2
N=3
tf = [I for I in Iterators.product(ntuple(d->1:N+1, D)...)] 
tf[1,2]
Ma[collect(tf[2,2])]
tf
Mi=[2 , 3]
Nt=ntuple(d->1:N+1, D)
Mt=ntuple(d->1:N+2, D)
1:5
Nt=(Nt...,1:5)
[I for I in Iterators.product(Nt...)]


for i,j in Nt
    print(i)
end
for 
[1;1]<[2;2]

l=[]
collect(i for i<5)
collect([i,j] for i in 1:3 for j in 1:2 if i!=2)



function createloops(n)
    if n==0
        return 
    end
    for i=1:n
        println(i)
        createloops(n-1)
    end
end
