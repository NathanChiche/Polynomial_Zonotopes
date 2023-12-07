
function get_polynomials_from_SSPZ(PZ::SimpleSparsePolynomialZonotope,field::Field)#checked

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

function get_SSPZ_from_polynomials(Polynomes)#il faut faire qqchose pour éviter les polynomes de dimesion 1
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

function SimpletoSPZ(SSPZ)
    nb_vars=size(SSPZ.E)[1]
    id=collect(1:nb_vars)
    zer=zeros(Float64,length(SSPZ.c),1)
    return SparsePolynomialZonotope(SSPZ.c,SSPZ.G,zer,SSPZ.E,id)
end

function SPZ_to_SimpleSPZ(PZ::SparsePolynomialZonotope)
    dim,nb_indep=size(PZ.GI)
    @assert dim==length(PZ.c) "erreur delogique avec la dimension"
    z=zeros(Float64,dim,nb_indep)
    if nb_indep==0 || z==PZ.GI
        return SimpleSparsePolynomialZonotope(PZ.c,PZ.G,PZ.E)
    end
    println("Mauvaise partie CONVERSION GI NON VIDE")
    Gnew=hcat(PZ.G,PZ.GI)
    n,m=size(PZ.E)
    Z=zeros(Int64,nb_indep,m)
    Enew=vcat(PZ.E,Z)
    
    v=zeros(Int64,nb_indep+n)
    v_temp=zeros(Int64,nb_indep+n)

    for j in 1:nb_indep
        v_temp=copy(v)
        v_temp[n+j]=1
        Enew=hcat(Enew,v_temp)
    end
    return SimpleSparsePolynomialZonotope(PZ.c,Gnew,Enew)
end

#=z=zeros(Float64,2)
G=[1.0 2 4; 2 2 2.0]
E=[1 0 2; 3 2 1; 4 4 4;1 1 1]
GI=[1 2.0; 2 2.0]
P=SparsePolynomialZonotope(z,G,GI,E)
SP=SPZ_to_SimpleSPZ(P)
SP.c
SP.G
P.G
P.GI
SP.E

a=2
b=3
b=a
b=b+1
b
a=#