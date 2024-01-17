
function createIndexes(maxdeg)
    indx=ntuple(d->1:maxdeg[1]+1,1)
    for i in 2:length(maxdeg)
        indx=(indx...,1:maxdeg[i]+1)
    end
    #return indx
    res=[I for I in Iterators.product(indx...)]
    res=collect(collect(i) for i in res)
    return res
end


function monobernsteinmultivar_from_univariates(monomial,maxdeg,domain,indexlist)
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


function findmaxdeg_listmono(Expo)
    n1,n2=size(Expo)
    maxdeg=Vector{Int64}(undef,n2)
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
    index=createIndexes(maxdegr)
    println(maxdegr)
    println(listmono)
    nbmonomials=size(Expo)[2]
    All=Array{Vector{Float64}}(undef,nbmonomials)
    for i in 1:nbmonomials
        All[i]=monobernsteinmultivar_from_univariates(listmono[i],maxdegr,domain,index)
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
    
    return ranges
end

G=[1 0; 1 2]
E=[0 1; 1 1]

ranges_from_Bernsteincoeff(G,E,dom)
dom

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


