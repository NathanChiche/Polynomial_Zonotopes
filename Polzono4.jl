using Nemo
using LazySets
using Plots
using Random 
using Symbolics
using Dates
using BenchmarkTools
#include("plot_polyzono.jl")

using TimerOutputs
#include("test.jl")
#using StatsPlots

#changement pour github

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


@timeit to function SPZ_to_SimpleSPZ(PZ::SparsePolynomialZonotope)
    dim,nb_indep=size(PZ.GI)
    
    z=zeros(Float64,dim,nb_indep)
    if z==PZ.GI
        return SimpleSparsePolynomialZonotope(PZ.c,PZ.G,PZ.E)
    end
    println("Mauvaise partie")
    z=zeros(Int64,nb_indep)
    Gnew=hcat(PZ.G,PZ.GI)
    n=size(PZ.E)[1]
    Enew=Array{Int64}(undef,nb_indep+n,0)
    
    v=zeros(Int64,nb_indep+n)
    for i in 1:size(PZ.E)[2]
        #println(size(Enew)[1])
        #println(size(PZ.E)[1])
        Enew=hcat(Enew,vcat(PZ.E[:,i],z))
    end

    for j in 1:nb_indep
        v_temp=v
        v_temp[n+j]=1
        Enew=hcat(Enew,v_temp)
    end
    return SimpleSparsePolynomialZonotope(PZ.c,Gnew,Enew)
end

@timeit to function remove_zero_lines(M::Matrix{Int64})
    z=zeros(size(M)[2])
    i=1
    while i<=size(M)[1]
        if M[i,:]==z
            if i==1
                M=M[2:end,:]
            elseif i==size(M)[1]
                M=M[1:end-1,:]
            else
                M=vcat(M[1:i-1,:],M[i+1:end,:])
            end
        else
            i=i+1
        end
    end
    return M
end

@timeit to function SimpleSPZ_to_SPZ(PZ::SimpleSparsePolynomialZonotope{Float64, Vector{Float64}, Matrix{Float64}, Matrix{Int64}})
    E=expmat(PZ)
    G=genmat(PZ)
    n=size(E)[1]
    m=size(G)[1]
    GI=Array{Float64}(undef,m,0)
    GD=Array{Float64}(undef,m,0)
    Enew=Array{Int64}(undef,n,0)
    len=size(E)[2]
    c=PZ.c
    zero=zeros(Int64,len)
    for i in 1:len
        expo=E[:,i]
        coeffs=G[:,i]
        simple=simple_exponent(expo)
       
        if expo==zero
            c=c+coeffs
        elseif simple>0 && zeros_except_index(E[simple,:],i)
            GI=hcat(GI,coeffs)
        else
            GD=hcat(GD,coeffs)
            Enew=hcat(Enew,expo)
        end
    end 
    #println(size(GD)[2],size(Enew)[2])
    #GD=remove_zero_lines(GD)
    Enew=remove_zero_lines(Enew)
    return SparsePolynomialZonotope(c,GD,GI,Enew)
end

function SimpletoSPZ(SSPZ)
    nb_vars=size(SSPZ.E)[1]
    id=collect(1:nb_vars)
    zer=zeros(Float64,length(SSPZ.c),1)
    return SparsePolynomialZonotope(SSPZ.c,SSPZ.G,zer,SSPZ.E,id)
end

#=R=RealField()
S,(x,y)=PolynomialRing(R,["x","y"])
lineaire1=-0.32*x+0.32*y
lineaire2= -0.42*x -0.92*y
Plineaire=get_SSPZ_from_polynomials([x,y])
image=poly_apply_on_SSPZ(Plineaire,[lineaire1,lineaire2],R)=#



@timeit to function simple_exponent(exponent::Vector{Int64})
    cp=0
    index=0
    for j in 1:length(exponent)
        e=exponent[j]
        cp=cp+e
        if e==1
            index=j
        end
        if cp>1
            return 0
        end
    end
    return index
end

@timeit to function zeros_except_index(line::Vector{Int64},index::Int64)
    for i in 1:length(line)
        if i!= index && line[i]!=0
            return false
        end
    end
    return true
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

@timeit to function get_SSPZ_from_polynomials(Polynomes)#il faut faire qqchose pour éviter les polynomes de dimesion 1
    """récupérer la forme PZ=<c,G,E> à partir de  PZ={(P1(x1,...xp),...,Pn(x1,...xp)) """

    liste_exp=Vector{Int64}[]
    n=length(Polynomes)
    #deg_list=Array{Int}(undef,n)
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
    #println("nb monomes dans getSSPZ from poly: ",nb_monomes)
    G=zeros(n,nb_monomes)
    #E=zeros([T={Int64}],nb_vars,nb_monomes)
    E=Array{Int64}(undef,nb_vars,nb_monomes)
    c=zeros(n)
    tmp=0
    for i in 1:nb_monomes
        for k in 1:nb_vars
            E[k,i]=liste_exp[i][k] # remplissage de la matrice des exposants
        end
        
        for j in 1:n
            exp=liste_exp[i]
            a=Float64(coeff(Polynomes[j],exp))
            if a==0
                tmp=tmp+1
            end  
            G[j,i]+=a #remplissage de la matrice génératrice 
        end
    end
    #println("nombre de coeffs nuls sur exposants: ",tmp)
    for i in 1:n
        c[i]=Float64(coeff(Polynomes[i],zero)) #remplissage ici du vecteur de centre 
    end

    #println("taille de E apres modifs dans getSSPZ: ",size(E)[2])
    return SimpleSparsePolynomialZonotope(c,G,E)
end

methods(get_SSPZ_from_polynomials)

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

@timeit to function poly_apply_on_SSPZ(PZ::SimpleSparsePolynomialZonotope,list_poly,field::Field)
    nb_vars=size(expmat(PZ))[1]#nombre de variables est le nombre de lignes de la matrice des exposants
    composit=[]
    #anneau,(x)=PolynomialRing(field,nb_vars)
    Poly_fromPZ=get_polynomials_from_SSPZ(PZ,field)
    #anneau=parent(Poly_fromPZ[1])
    
    for p in list_poly
        res=compose(p,Poly_fromPZ)
        push!(composit,res)
    end
    println("a l'interieur de la compo: ",length(composit[1]))
    #res=get_SSPZ_from_polynomials(composit)
    #println("apres changement de représentation: ",size(res.E)[2])
    return get_SSPZ_from_polynomials(composit)
end


@timeit to function affichematrice(M::Matrix{Any})
    for k in 1:size(M)[1]
        println(M[k,:])
    end
    println(" ")
end


@timeit to function inf_and_supbis(polynomial)
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

@timeit to function inf_and_sup(polynomial)
    exponents=collect(exponent_vector(polynomial,i) for i in 1:length(polynomial))
    l=length(exponents[1])
    zero=zeros(l)
    sup=Float64(evaluate(polynomial,zero))
    inf=Float64(evaluate(polynomial,zero))
    for e in exponents
        if e!=zero
            if even_exponent(e)
                sup=sup+abs(Float64(coeff(polynomial,e)))
            else
                sup=sup+abs(Float64(coeff(polynomial,e)))
                inf=inf-abs(Float64(coeff(polynomial,e)))
            end
        end
    end
    return [inf,sup]
end

p=x^2*y^4+y
el=collect(exponent_vector(p,i) for i in 1:length(p))
even_exponent(el[2])

@timeit to function mid_polynomials(p,q)
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

@timeit to function copy_poly(pol,Anneau)
    """on met le polynome pol dans un anneau multivarié Anneau ssi il est plus grand que parent(pol)"""
    n=length(gens(parent(pol)))
    m=length(gens(Anneau))
    #println("Anneau: ",Anneau)

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
    nb_vars=max(size(expmat(PZ1))[1],size(expmat(PZ2))[1])
    #println("nombre de variables avant de join: ",nb_vars)
    Polynomes1=get_polynomials_from_SSPZ(PZ1,field)
    #println("nombre de var: ", length(gens(parent(Polynomes1[1]))))
    Polynomes2=get_polynomials_from_SSPZ(PZ2,field)
    #println("nombre de var: ", length(gens(parent(Polynomes2[1]))))
    Anneau,(x)=PolynomialRing(field,nb_vars+1)
    Polynomes=[Anneau(0) for i in 1:nb]
    for i in 1:nb
        Polynomes[i]=(1/2*gens(Anneau)[nb_vars+1] + 1/2)*copy_poly(Polynomes1[i],Anneau) + copy_poly(Polynomes2[i],Anneau)*(1/2 - 1/2*gens(Anneau)[nb_vars+1])
    end
    return get_SSPZ_from_polynomials(Polynomes)
end

@timeit to function barycentric_join(SSPZ1,SSPZ2,field)
    G1=SSPZ1.G
    E1=SSPZ1.E
    G2=SSPZ2.G
    E2=SSPZ2.E
    c1=SSPZ1.c
    c2=SSPZ2.c
    n1=size(E1)[2]
    n2=size(E2)[2]
    center=ones(Int64,1,2)
    e1=zeros(Int64,1,n1)
    e_1=ones(Int64,1,n1)
    e2=zeros(Int64,1,n2)
    e_2=ones(Int64,1,n2)
    vecnew=hcat(center,e1,e_1,e2,e_2)
    nlines=size(E1)[1]
    t1=zeros(Int64,nlines,2)
    #Redundant=SimpleSparsePolynomialZonotope(0.5*c1+0.5*c2,hcat(0.5*c1,-0.5*c2,0.5*G1,0.5*G1,0.5*G2,-0.5*G2),vcat(hcat(t1,E1,E1,E2,E2),vecnew))

    return remove_redundant_generators(SimpleSparsePolynomialZonotope(0.5*c1+0.5*c2,hcat(0.5*c1,-0.5*c2,0.5*G1,0.5*G1,0.5*G2,-0.5*G2),vcat(hcat(t1,E1,E1,E2,E2),vecnew)))
end

PE=get_SSPZ_from_polynomials([x^2+y+3+x^7*y^4,x*y+x*x*y])
PR=get_SSPZ_from_polynomials([y+7*y^2,x*y*7*9*x+x])

PJ=barycentric_join(PE,PR,R)

PJ2=barycentre_union(PE,PR,R)


@timeit to function scale_SSPZ(factor::Float64,PZ::SimpleSparsePolynomialZonotope)
    G=factor*genmat(PZ)
    return SimpleSparsePolynomialZonotope(LazySets.center(PZ),G,expmat(PZ))
end


@timeit to function iterate_polynomials_over_PZ(Polynomes,PZ::SimpleSparsePolynomialZonotope,nb_iter::Int64,borne_union::Int64,field::Field,choice;max_order::Int64,toreduce::Int64=200,maxdegree::Int64=50,scale_factor::Float64=1.1,power::Int64=1)
    """il faudrait quand même trouver un moyen efficace de tester l'inclusion entre polynomial zonotopes"""
    i=0
    nb_reduc=0
    liste=[PZ]
    if choice=="zono"
        f=union_pol_zono
    elseif choice=="bary"
        f=barycentre_union
        #f=barycentric_join
    else
        print("Le join n'est pas disponible")
        return 0
    end
    while i<nb_iter
        #=if i%2==0
            f=barycentre_union
        else 
            f=union_pol_zono
        end =#
        
        println(i)
        println("coucou")
        println("nb variables PZ: ",size(PZ.E)[1])
        fPZ=liste[end]
        PZ_previous=fPZ
        for p in 1:power 
            fPZ=poly_apply_on_SSPZ(fPZ,Polynomes,field)
            
        end
        #println("nb variables PZ_interm: ",size(fPZ.E)[1])
        println("nombre de termes avant la réduction: ",size(fPZ.E)[2])
        polytest1=get_polynomials_from_SSPZ(fPZ,field)
        PZ_reduc=SimpletoSPZ(fPZ)
        if Float64(LazySets.order(PZ_reduc))>max_order
            println("voici l'ordre du SSPZ: ",LazySets.order(PZ_reduc))
            PZ_reduc=reduce_order(PZ_reduc,max_order)
            println("ordre après reduction: ",LazySets.order(PZ_reduc))
            nb_reduc=nb_reduc+1
        end
        fPZ=SPZ_to_SimpleSPZ(PZ_reduc)
        polytest2=get_polynomials_from_SSPZ(fPZ,field)

        if nb_reduc==0 && polytest1!=polytest2
            affiche_liste(polytest1)
            affiche_liste(polytest2)
            println("LOUPE")
            @show(polytest1)
            @show(polytest2)
            return liste
        end

        #plot_sampling(PZ_interm,field,filename*string(i)*".png")
        if i>= borne_union
            PZ=f(PZ_previous,fPZ,field)
        else
            PZ=fPZ
        end
        println("number of terms/monomials after join ",size(expmat(PZ))[2])
        #println("nb of terms bis p1: ",length(get_polynomials_from_SSPZ(PZ,field)[1]))
        #println("nb of terms bis p2: ",length(get_polynomials_from_SSPZ(PZ,field)[2]))
        #=if inclusion_test(get_polynomials_from_SSPZ(PZ,field),get_polynomials_from_SSPZ(PZ_interm,field))
            println("on a trouvé notre invariant")
            return PZ_interm
        end=#
       

        i+=1
        push!(liste,PZ)
    end
    println("nombre de réductions, ",nb_reduc)
    return liste
end

@timeit to function evaluate_polynomials_on_vector(polynomials,vector::Vector{Float64})
    return [Float64(evaluate(polynomials[i],vector)) for i in 1:length(polynomials)]
end

@timeit to function plot_sampling(PZ::SimpleSparsePolynomialZonotope,field::Field,filename::String;nbpoints=300000,xlim=nothing,ylim=nothing)
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

@timeit to function strict_iterates(Polynomes,PZ::SimpleSparsePolynomialZonotope,nb_iter::Int64,field::Field)
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


@timeit to function plot_multiple(liste_PZ,field::Field,filename::String;nbpoints=300000,xlim=nothing,ylim=nothing)
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

@timeit to function test_monomial_in_common(p,q)
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

@timeit to function affiche_polynomes(pol)
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

function affiche_liste(list_poly)
    for p in list_poly
        affiche_polynomes(p)
    end
end

@timeit to function main()
    R=RealField()
    S,(x,y)=PolynomialRing(R,["x","y"])

    p1=x^3 -0.5*x^2+0.5
    p2=y^3 -0.5*y^2+0.5

    p6=(3/5*x +4/5*y)^3 -0.5*(3/5*x +4/5*y)^2+0.5
    p7=((-4/5)*x+3/5*y)^3 -0.5*((-4/5)*x+3/5*y)^2+0.5

    chatal1= x+y
    chatal2=-0.5952 + x^2

    p8=1/4*x^2+1/2 -y
    p9=2*y-y^2 + x

    henon1=1 - 1.4*x^2 + y
    henon2=0.3*x

    lineaire1=-0.32*x+0.32*y
    lineaire2= -0.42*x -0.92*y

    step=0.25
    vanpol1=step*y + x
    vanpol2=step*(1-x^2)*y - step*x + y

    brusselator1=step*(x)
    brusselator2=step*(y)

    lokta1=step*(x*(1.5 - y)) + x
    lokta2=step*(-y*(3 - x)) + y

    parillo1= step*(-x+x*y) + x
    parillo2=step*(-y) + y

    PHenon=get_SSPZ_from_polynomials([1/5*x ,1/5*y])
    Pparillo=get_SSPZ_from_polynomials([1/10*x+12 ,1/10*y+2])
    PChatal=get_SSPZ_from_polynomials([1/6*x+1/6 ,1/6*y+1/6])
    Plineaire=get_SSPZ_from_polynomials([x-2,y+1])
    P1=get_SSPZ_from_polynomials([x + 3  ,y+4])
    PVanPol=get_SSPZ_from_polynomials([0.15*x + 1.4  ,0.05*y+2.30])
    Plokta=get_SSPZ_from_polynomials([1/5*x + 5  ,1/5*y+2])
    #Pquad= get_SSPZ_from_polynomials([1/2*x ,3/4*y])
    #affichematrice(expmat(P1))

    start_time = now()
    fin=iterate_polynomials_over_PZ([chatal1,chatal2],PChatal,4,2,R,"zono",max_order=25000000)
    end_time = now()
    elapsed = end_time - start_time
    println("temps des iterations:", elapsed)
    #fin=reverse(fin)
    fini=fin[end]
    println("nombre de variables à la fin: ",size((fini).E)[1])
    println("nombre de monomes à la fin: ",size(fini.E)[2])
    somme=sum(fini.E,dims=1)
    som=vec(somme)
    println("degré maximal à la fin: ",maximum(som))
    #affiche_liste(get_polynomials_from_SSPZ(fini,R))

    plot_multiple(fin,R,"Documents/julia/plots_julia/chatal^1_4iter_joinzonoraffine2_30000pts_max_order=40000000bis",nbpoints=30000)
    #plot_sampling(fini,R,"Documents/julia/plots_julia/Parillo^1_4iter_joinbary_100000pts_",nbpoints=100000)
    #plot_multiple(fin,R,"Documents/julia/plots_julia/lineaire^1_4iter_joinbary1_40000pts_max_order=25_x-2_y+1",nbpoints=40000)
    
    #derniere=poly_apply_on_SSPZ(fini,[lineaire1,lineaire2],R)

    #plot_multiple([derniere,fini],R,"Documents/julia/plots_julia/lineaire^1_5iter_joinbary0_Inclusion?_1000000pt",nbpoints=1000000)
    end_time2 = now()
    elapsed = end_time2 - end_time
    println("temps de plot:", elapsed)

    
    

    return fini
    #return derniere
   
    
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



r=main()
r.G
r.E
r.G[:,1:10]
r.G[:,11:20]
r.G[:,21:23]
r.G[1,:]

to

R=RealField()
S,(x,y)=PolynomialRing(R,["x","y"])
lineaire1=-0.32*x+0.32*y
lineaire2= -0.42*x -0.92*y
chatal1= x+y
chatal2=-0.5952 + x^2
#Plineaire=get_SSPZ_from_polynomials([x,y])
image=poly_apply_on_SSPZ(r,[lineaire1,lineaire2],R)
image.G[:,1:10]
image.G[:,11:23]
image.G[:,1:10]
image.c
chatal5=poly_apply_on_SSPZ(r,[chatal1,chatal2],R)
plot_multiple([r,image],R,"Documents/julia/plots_julia/lineaire^1_3iter_joinbary1_70000pts_max_order=25_x-2_y+1_inclusion?",nbpoints=70000)

G=image.G
0.5*G

GZ=hcat(G,G,G)
image2=SimpleSPZ_to_SPZ(image)
plot(image2,nsdiv=30)

image3=SPZ_to_SimpleSPZ(image2)
plot([image,image3],nsdiv=45)
image2

rad=rand(SparsePolynomialZonotope)

@show(reducto.G)




println("coucou")
#=r.G

to
agb=get_polynomials_from_SSPZ(r,R)
affiche_liste(agb)

R=RealField()

function testconversion(n::Int64)
    bool=true
    for i in 1:n
        ran=rand(SimpleSparsePolynomialZonotope)
        sparse=SimpleSPZ_to_SPZ(ran)
        simple=SPZ_to_SimpleSPZ(sparse)
        if ran!=simple
            println(1)
            @show(ran)
            @show(simple)
            affiche_liste(get_polynomials_from_SSPZ(ran,R))
            affiche_liste(get_polynomials_from_SSPZ(simple,R))
            
            return [ran,simple]
        end
    end
    for i in 1:n
        rand=rand(SparsePolynomialZonotope)
        simple=SPZ_to_SimpleSPZ(rand)
        sparse=SimpleSPZ_to_SPZ(simple)
        if rand!=sparse
            println(2)
            @show(rand)
            @show(sparse)
            return false
        end
    end
    return bool
end

li=testconversion(100)
plot_multiple(li,R,"Documents/julia/plots_julia/testconversion")

rando=rand(SimpleSparsePolynomialZonotope,num_generators=2)
rando2=SimpleSPZ_to_SPZ(rando)
rando3=SPZ_to_SimpleSPZ(rando2)
rando
rando==rando3

rspz=SimpleSPZ_to_SPZ(r)
LazySets.order(rspz)
res=reduce_order(rspz,20)
LazySets.order(res)
res.E

a=rand

plot()
plot(r)
plot(r,nsdiv=20)
R=RealField()
plot_sampling(r,R,"Documents/julia/plots_julia/testrotat_3iter_joinzono",100000)
plot(r)
r.E
r.G
Z = Zonotope([0.0, 0.0], [1 2.0 ; 1.0 4 ])
Z2=Zonotope([6.0, 0.0], [1 2.0 ; -3.0 1 ])
Z3=Zonotope([3.0, 0.0], [1 2.0 3 0; 0 1.0 0 4 ])
plot([Z,Z2,Z3])
plot(Z3)
#to

R=RealField()

S,(x,y)=PolynomialRing(R,["x","y"])
#P3D=get_SSPZ_from_polynomials([x,y,z])
#plot_sampling(P3D,R,"Documents/julia/plots_julia/testtroisdmension")
APRES=get_SSPZ_from_polynomials([x-x*y,y-x])
APRES2=get_SSPZ_from_polynomials([x-x*y,y-x])
APRES==APRES2
typeof(APRES)
random1=rand(SimpleSparsePolynomialZonotope)
random1.E
APRE=SimpleSPZ_to_SPZ(APRES)
LazySets.order(APRE)
AVANT=get_SSPZ_from_polynomials([x,y])
plot([AVANT,APRES],nsdiv=40)

S,(x,y,s,t)=PolynomialRing(R,["x","y","s","t"])
Join=get_SSPZ_from_polynomials([x-x*y+s,y-x*y+2*t])
plot([Join,AVANT,APRES],nsdiv=20)

typeof(S)

s=1/4*x^2 + 1/2 - y
t=2*y - y^2 + x
pointx=1/3
pointy=1/4
for cpt in 1:n
    tpointx=s(pointx,pointy)
    pointy=t(pointx,pointy)
    pointx=tpointx
end


p1=x^3 -0.5*x^2+0.5
p2=y^3 -0.5*y^2+0.5


function testdegre(n::Int64,pol,qol,q1,q2)
    for i in 1:n
        #println(q1)
        q1t=compose(pol,[q1,q2])
        q2=compose(qol,[q1,q2])
        q1=q1t
        #q1=pol(q1,q2)
        #q2=qol(q1,q2)
        print(i)
        println(": nombre de termes de p ->",length(q1))
        
    end
    return q1
end





P1=get_SSPZ_from_polynomials([p6,p7])
P2=get_SSPZ_from_polynomials([p1,p2])

SPZ1=get_SSPZ_from_polynomials([x,y])
SPZ2=get_SSPZ_from_polynomials([x+3,-y+3])
plot([SPZ1,SPZ2])
SPZ3=barycentre_union(SPZ1,SPZ2,R)
plot_sampling(SPZ3,R,"Documents/julia/plots_julia/pasconvexhull")=#






