
function createIndexes(maxdeg)
    #println(maxdeg)
    indx=ntuple(d->1:maxdeg[1]+1,1)
    for i in 2:length(maxdeg)
        indx=(indx...,1:maxdeg[i]+1)
    end
    #return indx
    res=[I for I in Iterators.product(indx...)]
    res=collect(collect(i) for i in res)
    return res
end


function monobernsteinmultivar_from_univariatesbis(monomial,maxdeg,domain,indexlist)
    #on obtient dans la liste coeffs les coeffs de Bernstein multivariés pour un monome précis sur un domaine précis
    
    ll=multivariate(monomial,maxdeg,domain)
    nbvar=length(maxdeg)
    nb_coeffs=prod(k+1 for k in maxdeg)
    coeffs=Vector{Float64}(undef,nb_coeffs)
    len=length(indexlist[1])
     for i in 1:length(indexlist)
        coeffs[i]=prod(ll[j][indexlist[i][j]] for j in 1:len)
    end
    return coeffs
end

function monobernsteinmultivar_from_univariates(monomial,maxdeg,domain)
    #on obtient dans la liste coeffs les coeffs de Bernstein multivariés pour un monome précis sur un domaine précis
    
    ll=multivariate(monomial,maxdeg,domain)
    
    coeffs=everyproducts(ll,maxdeg)
    return coeffs
end


function findmaxdeg_listmono(Expo)
    #@show(Expo)
    n1,n2=size(Expo)
    maxdeg=Vector{Int64}(undef,n1)
    for i in 1:n1
        maxdeg[i]=maximum(@view Expo[i,:])
    end

    @polyvar x[1:n1]
    listmonomials=Array{Monomial}(undef,n2)
    for i in 1:n2
        listmonomials[i]=Monomial(x,Expo[:,i])
    end

    return maxdeg,listmonomials
end


function every_multivariate_bernsteincoeff(Expo,domain)
    
    maxdegr,listmono=findmaxdeg_listmono(Expo)
    #index=createIndexes(maxdegr)
    nbmonomials=size(Expo)[2]
    All=Array{Vector{Float64}}(undef,nbmonomials)
    for i in 1:nbmonomials
        All[i]=monobernsteinmultivar_from_univariates(listmono[i],maxdegr,domain)
    end
    return All
end


function ranges_from_Bernsteincoeff(G,E,domain)
    dim=size(G)[1]
    coeffs=every_multivariate_bernsteincoeff(E,domain)
    nbmonomials=length(coeffs)
    polynomials_coeffs=Vector{Vector{Float64}}(undef,dim)
    ranges=Vector{Vector{Float64}}(undef,dim)

    for i in 1:dim
        #for j in 1:length(coeffs)
        polynomials_coeffs[i]=sum(G[i,j]*coeffs[j] for j in 1:nbmonomials)   
        ranges[i]=[minimum(polynomials_coeffs[i]),maximum(polynomials_coeffs[i])]
    end
    #@show(polynomials_coeffs)
    
    return ranges
end

listest=[[1,2,1],[1,2,4]#=,[0,1,3]=#]
listest[1]
a=[1,2,1]
b=[2,2,1]
[a[i]*b[j] for i in 1:3 for j in 1:3]
prod(a[i]*b[j] for i in 1:3 for j in 1:3)
length(listest)
everyproductsbis(listest)


function everyproducts(listlist)
    num=prod(length(k) for k in listlist)
    res=zeros(Float64,num)
    temp=zeros(Float64,num)
    len=length(listlist[1])
    for c in 1:len
        temp[c]=listlist[1][c]
    end
    for i in 2:length(listlist)
        t=length(listlist[i])
        for j in 1:len
            for h in 1:t
                res[(j-1)*t+h]=temp[j]*listlist[i][h]
            end
        end
        len=len*t
        temp=deepcopy(res)
    end
    return res
end

function everyproductsbis(listlist)
    res=[]
    temp=listlist[1]
    len=length(temp)
    for i in 2:length(listlist)
        t=length(listlist[i])
        res=[temp[j]*listlist[i][h] for j in 1:len for h in 1:t]
        len=len*t
        temp=res
    end
    return res
end

listest=[[1.0,2,1],[1,2,4],[0,1,2],[1,1,2],[111,123,4134],[0,1,2],[0,1,2],[0,1,2],[0,1,2],[0,1,2],[0,1,2],[0,1,2],[0,1,2]]
@time (everyproducts(listest))
3^12
"""everyproducts(listest)

 #TESTS POUR VERIFIER SI LES COEFFS DE BERNSTEIN SONT BIEN CALCULES
dom=IntervalBox(-1..1,2)
G=[1 -1 -1 1; 1 2 3 2]
length(G)
E=[0 0 2 2; 0 2 0 2]
G[1,2]
a=@view G[1,2] 

G2=[1 2 1 2 4 2 1 2 1;1 2 1 2 4 2 1 2 1]
E2=[2 2 2 1 1 1 0 0 0;2 1 0 2 1 0 2 1 0]

rang=ranges_from_Bernsteincoeff(G2,E2,dom)
rang[1][1]

findmaxdeg_listmono(E)

#ranges_from_Bernsteincoeff(G,E,dom)
dom
inte=IntervalBox(0..1,2)
@polyvar x[1:2]
test=Monomial(x,[1,1])
ber=multivariate(test,[1,1],inte)
everyproducts(ber)

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
[I for I in Iterators.product(Nt...)]"""


