#include("plot_polyzono.jl")
using Nemo
using LazySets
using Plots
using Random 
using Symbolics
using Dates
using BenchmarkTools

using TimerOutputs
#include("test.jl")
#using StatsPlots

const to = TimerOutput();
function arg_min(a::Float64,b::Float64)
    """arg_min classique entre deux éléments a et b: réel c dans [a;b] de module minimal"""
    if sign(a)!=sign(b)
        return 0
    end
    if sign(a)>=0 && sign(b)>=0
        return min(a,b)
    end
    return max(a,b)
end


function sup_concret(CX,PX,indice)
    """ sup de la concrétisation de la forme affine de l'indice "indice" donnée par les matrices CX et PX"""
    somme=CX[1,indice]
    #println(size(CX)[1])
    for i in 2:size(CX)[1]
        somme=somme+abs(CX[i,indice])
    end
    for j in 1:size(PX)[1]
        somme=somme+abs(PX[j,indice])
    end
    return somme
end


function inf_concret(CX,PX,indice)
    """ inf de la concrétisation de la forme affine de l'indice "indice" donnée par les matrices CX et PX"""
    somme=CX[1,indice]
    #println(size(CX)[1])
    for i in 2:size(CX)[1]
        somme=somme-abs(CX[i,indice])
    end
    for j in 1:size(PX)[1]
        somme=somme-abs(PX[j,indice])
    end
    return somme
end


function mid_concret(CX,PX,CY,PY,indice)
    """ milieu de l'union des concrétisations des formes affines de l'indice "indice" données par les matrices CX,PX pour l'une et CY,PY pour l'autre"""
    return (max(sup_concret(CX,PX,indice),sup_concret(CY,PY,indice))+min(inf_concret(CX,PX,indice),inf_concret(CY,PY,indice)))/2
end

function union(CX,PX,CY,PY)
    """CX et CY sont dans M(n+1,p); PX et PY sont dans M(m,p)
    le résultat de l'union sur les deux vecteurs affines X et Y est dnné par CZ dans M(n+1,p) et PZ dans M(m+p,p)
    on parle ici d'une union zonotopique"""
    n=size(CX)[1]-1
    p=size(CX)[2]
    m=size(PX)[1]
    CZ=zeros(n+1,p)
    PZ=zeros(m+p,p)
    #il faudrait ne pas toucher lkes vecteurs qui sont identiques
    for k in 1:p
        CZ[1,k]=mid_concret(CX,PX,CY,PY,k)
        for i in 2:(n+1)
            CZ[i,k]=arg_min(CX[i,k],CY[i,k])
        end
        if m!=0
            for j in 1:m
                PZ[j,k]=arg_min(PX[j,k],PY[j,k])
            end 
        end
    end
    for j in 1:p
        if m==0
            PZ[m+j,j]=max(sup_concret(CX,PX,j),sup_concret(CY,PY,j)) - CZ[1,j] - sum(abs(CZ[i,j]) for i in 2:n+1)
        else
            PZ[m+j,j]=max(sup_concret(CX,PX,j),sup_concret(CY,PY,j)) - CZ[1,j] - sum(abs(CZ[i,j]) for i in 2:n+1) - sum(abs(PZ[i,j]) for i in 1:m)
        end
    end
    return [CZ,PZ]
end


@timeit to function get_monomials_from_expmat(expmat::Matrix{Int64},anneau)#sert probablement à rien
    monomials=[]
    for i in 1:size(expmat)[2]
        m=1
        for j in 1:size(expmat)[1]
            m=m*gens(anneau)[j]^expmat[j,i]
        end
        if m==1
            append!(monomials,[0])
        else
            append!(monomials,[m])
        end
    end
    return monomials
end    


@timeit to function get_polynomials_from_SSPZ(PZ::SimpleSparsePolynomialZonotope,field::Field)
    # on récupère les polynomes P1,...,Pn issus de la forme PZ={(P1(x1,...xp),...,Pn(x1,...xp)) pour x dans la boule unité pour la distance max}
    c=LazySets.center(PZ)
    G=genmat(PZ)
    E=expmat(PZ)
    nb_vars=size(E)[1]
    anneau,(x)=PolynomialRing(field,nb_vars) 
    #monomials=get_monomials_from_expmat(E,anneau)
    #Polynomes=Array{Any}(undef,size(c)[1])
    Polynomes=[anneau(0) for k in 1:size(c)[1]]
    for i in 1:size(G)[1]
        Polynomes[i]+=c[i]
        for j in 1:size(G)[2]
            setcoeff!(Polynomes[i],E[:,j] , G[i,j])
        end
    end
    return Polynomes
end

@timeit to function get_SSPZ_from_polynomials(Polynomes::Vector{Nemo.AbstractAlgebra.Generic.MPoly{FieldElem}})#il faut faire qqchose pour éviter les polynomes de dimesion 1
    """récupérer la forme PZ=<c,G,E> à partir de  PZ={(P1(x1,...xp),...,Pn(x1,...xp)) """

    liste_exp=Vector{Int64}[]
    n=length(Polynomes)
    deg_list=Array{Int}(undef,n)
    for i in 1:n
        for j in 1:length(Polynomes[i])
            exp=exponent_vector(Polynomes[i],j)
            push!(liste_exp,exp)
        end
    end
    unique!(liste_exp)#on enlève les doublons dans la liste qui viennent du fait que des monomes interviennent dans plusieurs polynomes
    nb_vars=length(liste_exp[1])#on veut connaitre le nombre de variables dans l'anneau de polynomes pour les lignes de E 

    zero=zeros(Int64,nb_vars)
    deleteat!(liste_exp,findall(x->x==zero,liste_exp))
    nb_monomes=length(liste_exp)#on compte le nombre de monômes pour la taille de G après avoir enlevé le monôme constant
    G=zeros(n,nb_monomes)
    #E=zeros([T={Int64}],nb_vars,nb_monomes)
    E=Array{Int64}(undef,nb_vars,nb_monomes)
    c=zeros(n)
    for i in 1:nb_monomes
        for k in 1:nb_vars
            E[k,i]=liste_exp[i][k] # remplissage de la matrice des exposants
        end
        for j in 1:n
            exp=liste_exp[i]
            G[j,i]+=Float64(coeff(Polynomes[j],exp)) #remplissage de la matrice génératrice 
        end
    end
    for i in 1:n
        c[i]=Float64(coeff(Polynomes[i],zero)) #remplissage ici du vecteur de centre 
    end
    return SimpleSparsePolynomialZonotope(c,G,E)
end

@timeit to function variables_up_to(variables,exponent::Vector{Int64})
    if length(variables)!=length(exponent)
        println("probleme de longueur")
        return False
    end
    #println(variables[1])
    return prod(variables[i]^exponent[i] for i in eachindex(variables))
end

@timeit to function evaluate2(p::Nemo.AbstractAlgebra.Generic.MPoly{FieldElem},list_poly::Vector{Nemo.AbstractAlgebra.Generic.MPoly{FieldElem}})
    """p doit être dans un anneau de polynômes dont le nombre de variables est identique à la longueur de la liste de polynoômes"""
    l=exponents=collect(exponent_vector(p,i) for i in 1:length(p))
    res=0
    for exp in l
       res+=coeff(p,exp)*variables_up_to(list_poly,exp)
    end
    return res
end

@timeit to function compose(p,list_poly)
    return evaluate(p,list_poly)
end

@timeit to function poly_apply_on_SSPZ(PZ::SimpleSparsePolynomialZonotope,list_poly::Vector{Nemo.AbstractAlgebra.Generic.MPoly{FieldElem}},field::Field)
    nb_vars=size(expmat(PZ))[1]#nombre de variables est le nombre de lignes de la matrice des exposants
    composit=[]
    #anneau,(x)=PolynomialRing(field,nb_vars)
    Poly_fromPZ=get_polynomials_from_SSPZ(PZ,field)
    #anneau=parent(Poly_fromPZ[1])
    
    for p in list_poly
        res=compose(p,Poly_fromPZ)
        push!(composit,res)
    end
    return get_SSPZ_from_polynomials(composit)
end


@timeit to function affichematrice(M::Matrix{Any})
    for k in 1:size(M)[1]
        println(M[k,:])
    end
    println(" ")
end


@timeit to function inf_and_sup(polynomial::Nemo.AbstractAlgebra.Generic.MPoly{FieldElem})
    exponents=collect(exponent_vector(polynomial,i) for i in 1:length(polynomial))
    l=length(exponents[1])
    zero=zeros(l)
    sup=Float64(evaluate(polynomial,zero))
    inf=Float64(evaluate(polynomial,zero))
    for e in exponents
        if e!=zero
            sup=sup+abs(Float64(coeff(polynomial,e)))
            inf=inf-abs(Float64(coeff(polynomial,e)))
        end
    end
    return [inf,sup]
end

@timeit to function mid_polynomials(p::Nemo.AbstractAlgebra.Generic.MPoly{FieldElem},q::Nemo.AbstractAlgebra.Generic.MPoly{FieldElem})
    return Float64((max(inf_and_sup(p)[2],inf_and_sup(q)[2])+min(inf_and_sup(p)[1],inf_and_sup(q)[1]))/2)
end

@timeit to function union_pol_zono(PZ1,PZ2,field)
    #on peut soit faire à partir des polynomes soit des matrices
    #intersection des listes des exposants des polynomes dans PZ1 et PZ2
    #p=arg_min(coeff(exp),coeff(exp)))*monomial(exp)
    #on augmente la taille des polynomes en déclarant un nouvel anneau
    """Cette partie a l'air un peu lourde parce qu'il faut tout réordonner selon les monômes pour faire marcher les choses
    G1=genmat(PZ1)'
    G2=genmat(PZ2)'
    PX1=Matrix{Float64}(undef,0,0) 
    PX2=Matrix{Float64}(undef,0,0)
    G,P=union(G1,PX1,G2,PX2)
    PU=P'
    l=size(PU)[2]
    G3=hcat(G',PU) # ici on a la nouvelle matrice des coefficients i.e genmat(union)
    Id=
    E1=expmat(PZ1)
    E2=expmat(PZ2)""" 
    Polynomes1=get_polynomials_from_SSPZ(PZ1,field)
    Polynomes2=get_polynomials_from_SSPZ(PZ2,field)
    nb_vars=size(expmat(PZ1))[1]
    nb=size(genmat(PZ1))[1]
    Anneau,(x)=PolynomialRing(field,nb_vars+nb) # on sait qu'on va ajouter autant de variables qu'il y a de formes polynomiales dans le vecteur
    poly_union=[Anneau(0) for i in 1:nb]
    l=length(exponent_vector(Polynomes1[1],1))
    zero=zeros(Int64,l)
    for i in 1:nb
        exponents=collect(exponent_vector(Polynomes1[i],j) for j in 1:length(Polynomes1[i]))
        #poly_union[i]=sum(arg_min(Flot64(coeff(Polynomes1[i],e)),Float64(coeff(Polynomes2[i],e)))*monomial(e) for e in exponents if e!=zero) #on ajoute tous les termes avec monômes non triviaux
        for e in exponents
            if e!=zero
                setcoeff!(poly_union[i],vcat(e,zeros(Int64,nb)),arg_min(Float64(coeff(Polynomes1[i],e)),Float64(coeff(Polynomes2[i],e)))) #faire hcat avec zzros(nb)
            end
        end
        mid=mid_polynomials(Polynomes1[i],Polynomes2[i])
        poly_union[i]+=mid # on ajoute le centre
        #jusqu'ici tout est bon
        exp_res=collect(exponent_vector(poly_union[i],j) for j in 1:length(poly_union[i]))
        sup=Float64(max(inf_and_sup(Polynomes1[i])[2],inf_and_sup(Polynomes2[i])[2]))
        
        sum_abs_values=0
        for e in exp_res
            if e!=vcat(zero,zeros(Int64,nb))
                sum_abs_values+=abs(Float64(coeff(poly_union[i],e)))
            end
        end
        poly_union[i]+=gens(Anneau)[nb_vars+i]*(sup - sum_abs_values - mid)
    end
    #return poly_union ----->suremenet bcp mieux 
    return get_SSPZ_from_polynomials(poly_union)
end

@timeit to function copy_poly(pol::Nemo.AbstractAlgebra.Generic.MPoly{FieldElem},Anneau::Nemo.AbstractAlgebra.Generic.MPolyRing{FieldElem})
    """on met le polynome pol dans un anneau multivarié Anneau ssi il est plus grand que parent(pol)"""
    n=length(gens(parent(pol)))
    m=length(gens(Anneau))
    @assert m>=n "dimension problem "
    p=Anneau(0)
    exponents=collect(exponent_vector(pol,j) for j in 1:length(pol))
    zero=zeros(Int64,m-n)
    for e in exponents
        setcoeff!(p,vcat(e,zero),coeff(pol,e))
    end
    return p
end


@timeit to function barycentre_union(PZ1::SimpleSparsePolynomialZonotope,PZ2::SimpleSparsePolynomialZonotope,field::Field)
    """union barycentrique entre PZ1 et PZ2"""
    nb=size(genmat(PZ1))[1]
    nb_vars=size(expmat(PZ1))[1]
    Polynomes1=get_polynomials_from_SSPZ(PZ1,field)
    Polynomes2=get_polynomials_from_SSPZ(PZ2,field)
    Anneau,(x)=PolynomialRing(field,nb_vars+1)
    Polynomes=[Anneau(0) for i in 1:nb]
    for i in 1:nb
        Polynomes[i]=(1/2*gens(Anneau)[nb_vars+1] + 1/2)*copy_poly(Polynomes1[i],Anneau) + copy_poly(Polynomes2[i],Anneau)*(1/2 - 1/2*gens(Anneau)[nb_vars+1])
    end
    return get_SSPZ_from_polynomials(Polynomes)
end

@timeit to function scale_SSPZ(factor::Float64,PZ::SimpleSparsePolynomialZonotope)
    G=factor*genmat(PZ)
    return SimpleSparsePolynomialZonotope(LazySets.center(PZ),G,expmat(PZ))
end


@timeit to function iterate_polynomials_over_PZ(Polynomes::Vector{Nemo.AbstractAlgebra.Generic.MPoly{FieldElem}},PZ::SimpleSparsePolynomialZonotope,nb_iter::Int64,borne_union::Int64,field::Field,max_order::Int64,toreduce::Int64,maxdegree::Int64,scale_factor::Float64,choice::Bool=true)
    """il faudrait quand même trouver un moyen efficace de tester l'inclusion entre polynomial zonotopes"""
    i=0
    liste=[PZ]
    if choice==true
        f=union_pol_zono
    else
        f=barycentre_union
    end
    while i<nb_iter
        println(i)
        PZ_interm=poly_apply_on_SSPZ(PZ,Polynomes,field)
        PZ_interm=reduce_order_SSPZ(PZ_interm,max_order,toreduce,maxdegree,field)
        #plot_sampling(PZ_interm,field,filename*string(i)*".png")
        if i>= borne_union
            PZ=f(PZ,PZ_interm,field)
        else
            PZ=PZ_interm
        end
        println("number of terms ",size(expmat(PZ))[2])
        #=if inclusion_test(get_polynomials_from_SSPZ(PZ,field),get_polynomials_from_SSPZ(PZ_interm,field))
            println("on a trouvé notre invariant")
            return PZ_interm
        end=#

        """PZ=reduce_order_SSPZ(PZ,max_order,toreduce,maxdegree,field)"""

        #=if i%3==0 && i>0
            PZ=scale_SSPZ(scale_factor,PZ)
            println("scale passée")
        end=#
        i+=1
        push!(liste,PZ)
    end
    return liste
end

@timeit to function evaluate_polynomials_on_vector(polynomials::Vector{Nemo.AbstractAlgebra.Generic.MPoly{FieldElem}},vector::Vector{Float64})
    return [Float64(evaluate(polynomials[i],vector)) for i in 1:length(polynomials)]
end

@timeit to function plot_sampling(PZ::SimpleSparsePolynomialZonotope,field::Field,filename::String,nbpoints=300000,xlim=nothing,ylim=nothing)
    """enregistre dans filename le tracé de PZ avec nbpoints différents"""
    nb_vars=size(expmat(PZ))[1]
    nb=size(genmat(PZ))[1]
    polynomes=get_polynomials_from_SSPZ(PZ,field)
    step=1/nbpoints
    liste=collect(rand(-1.0:step:1.0,nb_vars) for i in 1:nbpoints)
    points = [evaluate_polynomials_on_vector(polynomes,v) for v in liste]
    #pl=plot(first.(points),last.(points),legend=false)
    x=collect(points[i][1] for i in 1:length(points))
    y=collect(points[i][2] for i in 1:length(points))
    pl=scatter(x,y,#=xlimits=xlim,ylimits=ylim,=#legend=false)
    savefig(pl,filename)
    return pl
end    

@timeit to function strict_iterates(Polynomes::Vector{Nemo.AbstractAlgebra.Generic.MPoly{FieldElem}},PZ::SimpleSparsePolynomialZonotope,nb_iter::Int64,field::Field)
    i=0
    liste=[PZ]
    PZ_interm=PZ
    plots=[]
    while i<nb_iter
        println(i)
        PZ_interm=poly_apply_on_SSPZ(PZ,Polynomes,field)
        PZ=PZ_interm
        i+=1
        push!(liste,PZ)
    end
    return liste
end


@timeit to function plot_multiple(liste_PZ::Vector{Nemo.AbstractAlgebra.Generic.MPoly{FieldElem}},field::Field,filename::String,xlim=nothing,ylim=nothing,nbpoints=300000)
    i=1
    for PZ in liste_PZ
        println("coucou")
        nb_vars=size(expmat(PZ))[1]
        nb=size(genmat(PZ))[1]
        polynomes=get_polynomials_from_SSPZ(PZ,field)
        step=1/nbpoints
        liste=collect(rand(-1.0:step:1.0,nb_vars) for i in 1:nbpoints)
        points = [evaluate_polynomials_on_vector(polynomes,v) for v in liste]
        #pl=plot(first.(points),last.(points),legend=false)
        x=collect(points[i][1] for i in 1:length(points))
        y=collect(points[i][2] for i in 1:length(points))
        if i==1
            if xlim!==nothing
                scatter(x,y,xlimits=xlim,ylimits=ylim,legend=false)
            else 
                scatter(x,y,legend=false)
            end
        else
            if xlim!==nothing
                scatter!(x,y,xlimits=xlim,ylimits=ylim,legend=false)
            else 
                scatter!(x,y,legend=false)
            end
        end
        i=i+1
    end
    savefig(filename)
end

@timeit to function test_monomial_in_common(p::Nemo.AbstractAlgebra.Generic.MPoly{FieldElem},q::Nemo.AbstractAlgebra.Generic.MPoly{FieldElem})
    exponents=collect(exponent_vector(p,j) for j in 1:length(p))
    for e in exponents
        if coeff(p,e)!=0 && coeff(q,e)!=0
            println(e)
            return true
        end
    end
    return false
end

@timeit to function common_monomial_for_SPZ(PZ1::SimpleSparsePolynomialZonotope,PZ2::SimpleSparsePolynomialZonotope,field::Field)
    #va pas marcher pour le moment
    Polynomes1=get_polynomials_from_SSPZ(PZ1,field)
    Polynomes2=get_polynomials_from_SSPZ(PZ2,field)
    for i in eachindex(Polynomes1)
        if test_monomial_in_common(Polynomes1[i],Polynomes2[i])
            return true
        end
    end
    return false
end

@timeit to function even_exponent(expo::Vector{Int64})
    for i in expo
        if i%2!=0
            return false
        end
    end
    return true
end

@timeit to function reduce_order_SSPZ(PZ::SimpleSparsePolynomialZonotope,ordermax::Int64,toreduce::Int64,maxdegree::Int64,field::Field)
    """on conserve les monomes dont les coefficients sont grands ou ceux dont les exposants sont plus petits que maxdegree"""
    
    nb_vars=size(expmat(PZ))[1]
    Anneau,(x)=PolynomialRing(field,nb_vars+size(genmat(PZ))[1])
    Polynomes=get_polynomials_from_SSPZ(PZ,field)
    l=length(Polynomes)
    liste_exp=[]
    for i in 1:l
        for j in 1:length(Polynomes[i])
            push!(liste_exp,exponent_vector(Polynomes[i],j))
        end
    end
    maxi=0
    degreetoobig=[]
    for e in liste_exp
        m=sum(e)
        if m>maxdegree
            push!(degreetoobig,e)
        end
    end

    if size(expmat(PZ))[2]<=toreduce && length(degreetoobig)==0
        return PZ   
    end
    #println("reduction")
   
    big_expo=Array{Any}(undef,ordermax)
    big_coeffs=zeros(Float64,ordermax)
    cpt=1
    c=0
    for e in liste_exp
        c=maximum([(Float64(coeff(Polynomes[i],e))) for i in 1:length(Polynomes)])
        if cpt<ordermax
            big_expo[cpt]=e
            big_coeffs[cpt]=c
        elseif cpt==ordermax # on rajoute le dernier coeff/exposant et on trie pour mettre le + pletit coeff en dernier
            big_expo[cpt]=e
            big_coeffs[cpt]=c
            (min,ind)=findmin(big_coeffs)
            temp_e=big_expo[ordermax]
            temp_c=big_coeffs[ordermax]
            big_expo[ordermax]=big_expo[ind]
            big_coeffs[ordermax]=big_coeffs[ind]
            big_expo[ind]=temp_e
            big_coeffs[ind]=temp_c
        else
            if c>big_coeffs[ordermax]
                big_expo[ordermax]=e
                big_coeffs[ordermax]=c
                (min,ind)=findmin(big_coeffs)
                temp_e=big_expo[ordermax]
                temp_c=big_coeffs[ordermax]
                big_expo[ordermax]=big_expo[ind]
                big_coeffs[ordermax]=big_coeffs[ind]
                big_expo[ind]=temp_e
                big_coeffs[ind]=temp_c
            end
        end
        cpt=cpt+1
    end
    final_polys=[Anneau(0) for i in 1:l]
    zero=zeros(Int64,l)
    for i in 1:l
        p=Polynomes[i]
        sum=0
        for e in liste_exp
            if e in big_expo && !(e in degreetoobig)
                setcoeff!(final_polys[i],vcat(e,zero),coeff(p,e))
            elseif even_exponent(e)
                #println("exposant ",e)
                #probleme ici a cause du repassage meme quand on a un coeff qui vaut 0
                c=Float64(coeff(p,e))
                if c!=0
                    final_polys[i]=final_polys[i]+1/2
                    sum=sum+1/2*(Float64(coeff(p,e)))
                end
            else 
                sum=sum+Float64(coeff(p,e))
            end
        end
        final_polys[i]+=gens(Anneau)[nb_vars+i]*(sum)
    end
    return get_SSPZ_from_polynomials(final_polys)
end

@timeit to function affiche_polynomes(pol::Nemo.AbstractAlgebra.Generic.MPoly{FieldElem})
    #exponents=collect(exponent_vector(pol,j) for j in 1:length(pol))
    str=""
    println(str)
    for i in 1:length(pol)
        c=coeff(pol,i)
        m=monomial(pol,i)
        if sign(Float64(c))>= 0
            #println(string(Float64(c)),Float64(c))
            str=str*" + "*string(Float64(c))*string(m)
        else
            str=str*string(Float64(c))*string(m)
        end
    end
    println(str)
end
function affiche_liste(list_poly::Vector{Nemo.AbstractAlgebra.Generic.MPoly{FieldElem}})
    for p in list_poly
        affiche_polynomes(p)
    end
end

@timeit to function main()
    R=RealField()
    S,(x,y)=PolynomialRing(R,["x","y"])

    p1=x^3 -0.5*x^2+0.5
    p2=y^3 -0.5*y^2+0.5
    #p2=x*y^2
    p4=3/5*p1+4/5*p2
    p5=(-4/5)*p1+3/5*p2
    

    p6=(3/5*x +4/5*y)^3 -0.5*(3/5*x +4/5*y)^2+0.5
    p7=((-4/5)*x+3/5*y)^3 -0.5*((-4/5)*x+3/5*y)^2+0.5


    P1=get_SSPZ_from_polynomials([1/20*x + 0.4 ,1/20*y+0.40])
    #affichematrice(expmat(P1))

    start_time = now()
    fin=iterate_polynomials_over_PZ([p1,p2],P1,4,0,R,1500000000,2000000000,12000000000,1.1,false)
    end_time = now()
    elapsed = end_time - start_time
    println("temps des iterations:", elapsed)
    fin=reverse(fin)
    println("calcul passé")
   
    
    #savefig(pl,"Documents/essai6iterees")
    #plot([fin[1],fin[2],fin[3]],nsdiv=30)
    #end_time = now()
    #elapsed = end_time - start_time
    #println("temps de plot:", elapsed)
    #plot_sampling(fin[34],R,"Documents/34ieme")
    #@show(get_polynomials_from_SSPZ(fin[35],R))
    #plot_multiple(fin,R,"Documents/autresysteme_joinbary_3_8iter")
    #plot_sampling(fin[2],R,"Documents/enfinpolynomial")
  
end



main()
#to

R=RealField()

S,(x,y)=PolynomialRing(R,["x","y"])
typeof(S)

str="documents/cd"
typeof(str)
p1=x^3 -0.5*x^2+0.5
p2=y^3 -0.5*y^2+0.5
typeof(p1)
typeof(l)
    #p2=x*y^2
p4=3/5*p1+4/5*p2
p5=(-4/5)*p1+3/5*p2

p6=(3/5*x +4/5*y)^3 -0.5*(3/5*x +4/5*y)^2+0.5
p7=((-4/5)*x+3/5*y)^3 -0.5*((-4/5)*x+3/5*y)^2+0.5


P1=get_SSPZ_from_polynomials([p6,p7])
P1
reduce_order(P1,2)
E=expmat(P1)
typeof(E)==Matrix{Int64}
fin=iterate_polynomials_over_PZ([p1,p2],P1,2,0,R,25000,300,5000,1.1,false)
@show(get_polynomials_from_SSPZ(fin[1],R))
@show(get_polynomials_from_SSPZ(fin[2],R))

affiche_liste(get_polynomials_from_SSPZ(P1,R))
affiche_liste(get_polynomials_from_SSPZ(fin[1],R))
affiche_liste(get_polynomials_from_SSPZ(fin[2],R))
affiche_liste(get_polynomials_from_SSPZ(fin[3],R))

plot_sampling(fin[length(fin)],R,"Documents/classique_joinbary_2_iterations_resultatbis")



#Q=RationalField(16)
S,(x,y,z,t)=PolynomialRing(R,["x","y","z","t"])

P=get_SSPZ_from_polynomials([14*x^37+y^42+7*x+12*x*y*z+3*t^2,x*y+z^3])

@show(get_polynomials_from_SSPZ(P,R))

P2=reduce_order_SSPZ(P,4,3,4,R)
@show(get_polynomials_from_SSPZ(P2,R))
plot_sampling(P2,R,"Documents/pitie")
plot_sampling(P,R,"Documents/pitiedepart")
pl=plot(P,nsdiv=20)
savefig(pl,"Documents/pitiedepa")

p=x*y^2+t
q=x*y
list=[]
for i in 1:length(p)
    push!(list,exponent_vector(p,i))
end
e=list[1]
sum(e)

Ann=parent(p)
Ann
A,(a,b)=PolynomialRing(R,["a","b"])
pol=a^2+a*b
pol([p,q])
poly=evaluate(pol,[p,q])
parent(pol)
poly3=pol(p,q)
parent(poly3)
