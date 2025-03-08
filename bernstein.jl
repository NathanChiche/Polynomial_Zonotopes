
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

ar=createIndexes([2,2])

function list_bernstein_polyn(geners,maxdeg,domain)
    listlist=[]
    temp=[]
    for i in 1:length(maxdeg)
        l=maxdeg[i]
        up=domain[i].hi
        lo=domain[i].lo
        for j in 0:l
            push!(temp,(binomial(l,j)/(up-lo)^l)*(up-geners[i])^(l-j)*(geners[i]-lo)^j)
        end
        push!(listlist,temp)
        temp=[]
    end
    return listlist
end

function list_nemo_unibernstein(geners,maxdeg,domain)
    """on trouve tous les polynomes de bernstein univariés selon leurs degrés et domaines respectifs"""
    listlist=[]
    temp=[]
    for i in 1:length(maxdeg)
        l=maxdeg[i]
        up=domain[i].hi
        lo=domain[i].lo
        for j in 0:l
            push!(temp,(binomial(l,j)/(up-lo)^l)*(up-geners[i])^(l-j)*(geners[i]-lo)^j)
        end
        push!(listlist,temp)
        temp=[]
    end
    return listlist
end

function list_nemo_multivbernstein(geners,maxdeg,domain)

    ll=list_nemo_unibernstein(geners,maxdeg,domain)
    res=everyproductsbis(ll)
    return res
end


function polynomial_from_bernstein_coeffs(Anneau,maxdeg,domain,listcoeffs)
    geners=gens(Anneau)
    nemobernstein=list_nemo_multivbernstein(geners,maxdeg,domain)
    polynomeres=Anneau(0)
    for i in 1:length(listcoeffs)
        polynomeres+=listcoeffs[i]*nemobernstein[i]
    end
    return polynomeres
end

function monobernsteinmultivar_from_univariatesbis(monomial,maxdeg,domain,indexlist)
    #on obtient dans la liste coeffs les coeffs de Bernstein multivariés pour un monome précis sur un domaine précis
    println("bis")
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
    
    coeffs=everyproducts(ll)
    return coeffs
end

domain

function maxdeg(Expo)
    n1=size(Expo)[1]
    maxdeg=Vector{Int64}(undef,n1)
    for i in 1:n1
        maxdeg[i]=maximum(@view Expo[i,:])
    end
    return maxdeg
end

function list_monomials(Expo)
    n1,n2=size(Expo)
    @polyvar x[1:n1]
    listmonomials=Array{any}(undef,n2)
    for i in 1:n2
        listmonomials[i]=MultivariatePolynomials.monomial(x,Expo[:,i])
    end

    return listmonomials
end

function make_monomial(vars, exponents::AbstractVector{Int})
    prod(vars[i]^exponents[i] for i in eachindex(vars))
end

function findmaxdeg_listmono(Expo)
    #@show(Expo)
    n1,n2=size(Expo)
    maxdeg=Vector{Int64}(undef,n1)
    for i in 1:n1
        maxdeg[i]=maximum(@view Expo[i,:])
    end

    @polyvar x[1:n1]
    #listmonomials=Array{any}(undef,n2)
    listmonomials=[]
    for i in 1:n2
        #listmonomials[i]=MultivariatePolynomials.monomial(x,Expo[:,i])
        push!(listmonomials,make_monomial(x,Expo[:,i]))
    end

    return maxdeg,listmonomials
end

function every_multivariate_bernsteincoeff(Expo,domain,maxdegr,listmono)
    
    #maxdegr,listmono=findmaxdeg_listmono(Expo)
    #index=createIndexes(maxdegr)
    nbmonomials=size(Expo)[2]
    #@show(Expo)
    All=Array{Vector{Float64}}(undef,nbmonomials)
    #@show(listmono)
    for i in 1:nbmonomials
        All[i]=monobernsteinmultivar_from_univariates(listmono[i],maxdegr,domain)
    end
    return All
end


function convergencebetweenlists(P1,P2)

    """Need to give an alternative with the maxdeg as an argument to compare """
    domain=IntervalBox(-1..1,size(P1.E)[2])
    ll1=ranges_from_Bernsteincoeff(P1.G,P1.E,domain)[2]
    l=length(ll1)
    for i in 1:l
        ll1[i]=ll1[i].+P1.c[i]
    end
    domain=IntervalBox(-1..1,size(P2.E)[2])
    ll2=ranges_from_Bernsteincoeff(P2.G,P2.E,domain)[2]
    for i in 1:l
        ll2[i]=ll1[i].+P2.c[i]
    end
    R=RealField()
    An,(x,y)=polynomial_ring(R,["x","y"])
    diff=Array{Vector{Float64}}(undef,l)
    absolconv=Array{Float64}(undef,l)
    for i in 1:l
        #println("VOOIR ICI")
        #println("VOOIR ICI")
        #@show(ll1[i])
        #@show(polynomial_from_bernstein_coeffs(An,[2,2],domain,ll1[i]))
        #@show(ll2[i])
        #@show(polynomial_from_bernstein_coeffs(An,[2,2],domain,ll2[i]))
        #@show(map(abs,(ll1[i]-ll2[i])))
        diff[i]=map(abs,(ll1[i]-ll2[i]))
        #@show(diff[i])
        absolconv[i]=maximum(diff[i])
    end
    return diff,absolconv
end

map(abs,[1,2,3]-[2,2,2])



function ranges_from_Bernsteincoeff(G,E,domain;maxdeg=nothing)
    #@show(E)
    if maxdeg ===nothing
        maxdegr,listmono=findmaxdeg_listmono(E)
    else
        maxdegr=maxdeg
        listmono=list_monomials(E)
    end
    #maxdegr,listmono=findmaxdeg_listmono(E)
    dim=size(G)[1]
    coeffs=every_multivariate_bernsteincoeff(E,domain,maxdegr,listmono)
    nbmonomials=length(coeffs)
    polynomials_coeffs=Vector{Vector{Float64}}(undef,dim)
    ranges=Vector{Vector{Float64}}(undef,dim)

    for i in 1:dim
        #for j in 1:length(coeffs)
        polynomials_coeffs[i]=sum(G[i,j]*coeffs[j] for j in 1:nbmonomials)   
        ranges[i]=[minimum(polynomials_coeffs[i]),maximum(polynomials_coeffs[i])]
    end
    #@show(polynomials_coeffs)
    
    return ranges,polynomials_coeffs
end




"""listest=[[1,2,1],[1,2,4]#=,[0,1,3]=#]
listest[1]
a=[1,2,1]
prod(b)
b=[2,2,1]
[a[i]*b[j] for i in 1:3 for j in 1:3]
prod(a[i]*b[j] for i in 1:3 for j in 1:3)
length(listest)
everyproductsbis(listest)"""


function everyproductscaput(listlist)
    num=[length(k) for k in listlist]
    p=prod(num)
    l=length(listlist)
    res=zeros(Float64,p)
    temp=Vector{Vector{Float64}}(undef,l)
    len=length(listlist[1])
    #=for c in 1:len
        temp[1][c]=listlist[1][c]
    end=#
    temp[1]=deepcopy(listlist[1])
    for i in 2:l
        t=length(listlist[i])
        for j in 1:len
            for h in 1:t
                res[(j-1)*t+h]=temp[i-1][j]*listlist[i][h]
            end
        end
        len=len*t
        if i<l
            temp[i]=deepcopy(res)
        end
    end
    return res
end

function everyproducts(listlist)
    num=prod(length(k) for k in listlist)
    res=zeros(Float64,num)
    temp=zeros(Float64,num)
    len=length(listlist[1])
    for c in 1:len
        temp[c]=listlist[1][c]
    end
    l=length(listlist)
    for i in 2:l
        t=length(listlist[i])
        for j in 1:len
            for h in 1:t
                res[(j-1)*t+h]=temp[j]*listlist[i][h]
            end
        end
        len=len*t
        #recopie de res dans temp
        if i<l
            for s in 1:len
                temp[s]=res[s]
            end
        end
    end
    return res
end


function everyproductsb(listlist,anneau)
    num=prod(length(k) for k in listlist)
    res=zeros(Float,num)
    type=typeof(anneau(0))
    temp=zeros(type,num)
    len=length(listlist[1])
    for c in 1:len
        temp[c]=listlist[1][c]
    end
    l=length(listlist)
    for i in 2:l
        t=length(listlist[i])
        for j in 1:len
            for h in 1:t
                res[(j-1)*t+h]=temp[j]*listlist[i][h]
            end
        end
        len=len*t
        #recopie de res dans temp
        if i<l
            for s in 1:len
                temp[s]=res[s]
            end
        end
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




"""listest=[[1.0,2,1],[1,2,4],[2,4,6,1,77,2],[0,1,2],[1,1,2],[111,123,4134],[0,1,2],[0,1,2],[0,1,2],[0,1,2],[0,1,2],[0,1,2],[0,1,2],[0,1,2],[1,3,2,2222222,1]]
prod(length(k) for k in listest)
listest=[[1.0,2,1],[1,2],[0,1]]
tem=Vector{Vector{Float64}}(undef,length(listest))
tem[1]=deepcopy(listest[1])
tem[1][2]
everyproducts(listest)==everyproductsbis(listest)==everyproductscaput(listest)
@time(everyproductsbis(listest))
"""

"""everyproducts(listest)

 #TESTS POUR VERIFIER SI LES COEFFS DE BERNSTEIN SONT BIEN CALCULES
 using IntervalArithmetic
 using DynamicPolynomials
 using BernsteinExpansions
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
test=monomial(x,[1,1])
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
