
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

function even_exponent(expo::Vector{Int64})
    for i in expo
        if i%2!=0
            return false
        end
    end
    return true
end

function inf_and_supbis(polynomial)
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

function inf_and_sup(polynomial)
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

function mid_polynomials(p,q)
    return Float64((max(inf_and_sup(p)[2],inf_and_sup(q)[2])+min(inf_and_sup(p)[1],inf_and_sup(q)[1]))/2)
end

function union_pol_zono(PZ1,PZ2,field)
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

function copy_poly(pol,Anneau)
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


function barycentre_union(PZ1::SimpleSparsePolynomialZonotope,PZ2::SimpleSparsePolynomialZonotope,field::Field)
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

function barycentric_join(SSPZ1,SSPZ2,field)
    G1=SSPZ1.G
    E1=SSPZ1.E
    G2=SSPZ2.G
    E2=SSPZ2.E
    c1=SSPZ1.c
    c2=SSPZ2.c
    m1,n1=size(E1)
    m2,n2=size(E2)
    m=max(m1,m2)
    center=ones(Int64,1,2)
    e1=zeros(Int64,1,n1)
    e_1=ones(Int64,1,n1)
    e2=zeros(Int64,1,n2)
    e_2=ones(Int64,1,n2)
    vecnew=hcat(center,e1,e_1,e2,e_2)
    if m2>m1
        adjustment=zeros(Int64,m2-m1,n1)
        E1=vcat(E1,adjustment)
    elseif m1>m2
        adjustment=zeros(Int64,m1-m2,n2)
        E2=vcat(E2,adjustment)
    end
    t1=zeros(Int64,m,2)
    Redundant=SimpleSparsePolynomialZonotope(0.5*c1+0.5*c2,hcat(0.5*c1,-0.5*c2,0.5*G1,0.5*G1,0.5*G2,-0.5*G2),vcat(hcat(t1,E1,E1,E2,E2),vecnew))
    res=remove_redundant_generators(Redundant)
    return remove_redundant_generators(res)#on remove deux fois car sinon le premier remove peut donner des colonnes nulles dans G
end



#=using LazySets
using Nemo

R=RealField()
S,(x,y)=PolynomialRing(R,["x","y"])
include("polynomap.jl")

chatal1= x^2*y^2
chatal2=x^2

P=get_SSPZ_from_polynomials([x,y]);P=poly_apply_on_SSPZ(P,[chatal1,chatal2],R)
P=poly_apply_on_SSPZ(P,[chatal1,chatal2],R)
P.E
P3.E
P4=union_pol_zono(P2,P3,R)
P4.G
P4.E
P4.c=#