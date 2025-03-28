
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

function inf_and_sup(polynomial)
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

function inf_and_supbis(polynomial)
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
    Anneau,(x)=polynomial_ring(field,nb_vars+nb) # on sait qu'on va ajouter autant de variables qu'il y a de formes polynomiales dans le vecteur
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


function zonotopic_joinbisold(PZ1,PZ2,solver)
    """warning, it is required that PZ1 and PZ2 are defined over the same variables"""

    if solver=="NaturalEnclosure"
        sol=NaturalEnclosure()
    elseif solver=="BranchAndBoundEnclosure"
        sol=BranchAndBoundEnclosure()
    elseif solver=="TaylorModel"
        sol=TaylorModelsEnclosure()
    elseif solver=="SumOfSquares"
        sol=SumOfSquaresEnclosure()
    end

    c1=PZ1.c 
    c2=PZ2.c 
    functionsPZ1=[]
    functionsPZ2=[]
    n1,n2=size(PZ1.G)
    m1,m2=size(PZ2.G)
    e1,e2=size(PZ1.E)
    f1,f2=size(PZ2.E)
    nbvar=max(e1,f1)
    #Gnew=zeros(Float64,n1,n2+n1)
    #Enew=zeros(Int64,n1+e1,e2+n1)
    Gnew=Matrix{Float64}(undef,n1,0) 
    Enew=Matrix{Int64}(undef,n1+e1,0) 
    #println(typeof(Gnew))
    for i in 1:n2
        for j in 1:m2
            if PZ1.E[:,i]==PZ2.E[:,j]
                #Enew[:,i]=vcat(PZ1.E[:,i],zeros(n1))#add the dimension 
                #Gnew[:,i]=[arg_min(PZ1.G[k,i],PZ2.G[k,j]) for k in 1:n1]
                Enew=hcat(Enew,vcat(PZ1.E[:,i],zeros(Int64,n1)))
                Gnew=hcat(Gnew,[Float64(arg_min(PZ1.G[k,i],PZ2.G[k,j])) for k in 1:n1])
            end
        end
    end
    #println(typeof(Gnew))
    ranges1=[]
    domain=IntervalBox(-1..1, e1)
    for i in 1:n1
        a1(x)=sum(PZ1.G[i,j]*prod(x[k]^PZ1.E[k,j] for k in 1:e1) for j in 1:n2)+PZ1.c[i]
        push!(ranges1,enclose(a1,domain,sol))
        #println("enclose de la premiere partie :",enclose(a1,domain,BranchAndBoundEnclosure()))
    end
    println("ranges1: ",ranges1)
    domain=IntervalBox(-1..1, f1)
    ranges2=[]
    for i in 1:m1
        a2(x)=sum(PZ2.G[i,j]*prod(x[k]^PZ2.E[k,j] for k in 1:f1) for j in 1:m2)+PZ2.c[i]
        push!(ranges2,enclose(a2,domain,sol))
        #println("enclose de la deuxieme partie :",enclose(a2,domain,BranchAndBoundEnclosure()))
    end
    #println("ranges2: ",ranges2)
    center=zeros(n1)
    #ln(typeof(Gnew))
    for i in 1:n1
        mini=min(ranges1[i].lo,ranges2[i].lo)
        maxi=max(ranges1[i].hi,ranges2[i].hi)
        mid=Float64(1/2*(mini+maxi))
        center[i]=mid
        z=zeros(Int64,n1+e1)
        z[i+e1]=1
        Enew=hcat(Enew,z)
        z=zeros(Float64,n1)
        z[i]=Float64(maxi-mid-sum(abs(Gnew[i,l]) for l in 1:size(Gnew)[2]))
        #println(typeof(z))
        Gnew=hcat(Gnew,z)
    end

    return SimpleSparsePolynomialZonotope(center,Gnew,Enew)
end

function zonotopic_joinbis(PZ1,PZ2,solver,tolerance,maxdept)#Checked
    """warning, it is required that PZ1 and PZ2 are defined over the same variables"""
    if solver=="NaturalEnclosure"
        sol=NaturalEnclosure()
    elseif solver=="BranchAndBoundEnclosure"
        sol=BranchAndBoundEnclosure(tol=tolerance,maxdepth=maxdept)
    elseif solver=="TaylorModel"
        sol=TaylorModelsEnclosure()
    elseif solver=="SumOfSquares"
        sol=SumOfSquaresEnclosure(; backend=backend)
    end

    n1,n2=size(PZ1.G)
    m1,m2=size(PZ2.G)
    e1,e2=size(PZ1.E)
    f1,f2=size(PZ2.E)
    if (e1!=f1)
        throw(DomainError(e1, "nombre de variables différent"))
    end

    Gnew=Matrix{Float64}(undef,n1,0) 
    Enew=Matrix{Int64}(undef,n1+e1,0) 
    GP1=Matrix{Float64}(undef,n1,n2)
    GP2=Matrix{Float64}(undef,m1,m2)
    interm=zeros(Float64,n1)
    test=1
    for i in 1:n2
        for j in 1:m2
            if PZ1.E[:,i]==PZ2.E[:,j]
                test=0
                #Enew[:,i]=vcat(PZ1.E[:,i],zeros(n1))#add the dimension 
                #Gnew[:,i]=[arg_min(PZ1.G[k,i],PZ2.G[k,j]) for k in 1:n1]
                Enew=hcat(Enew,vcat(PZ1.E[:,i],zeros(Int64,n1)))
                Gnew=hcat(Gnew,[Float64(arg_min(PZ1.G[k,i],PZ2.G[k,j])) for k in 1:n1])
            end
        end
    end
    nh=size(Enew)[2]
    GP1=fill_matrix_after_argmin(PZ1.G,PZ1.E,Gnew,Enew,n2,nh,e1)
    GP2=fill_matrix_after_argmin(PZ2.G,PZ2.E,Gnew,Enew,m2,nh,f1)

    ranges1=[]
    domain=IntervalBox(-1..1, e1)
    #@polyvar y[1:e1]
    for i in 1:n1
        a1(x)=sum(GP1[i,j]*prod(x[k]^PZ1.E[k,j] for k in 1:e1) for j in 1:n2)+PZ1.c[i]
        #println("premier reste")
        #ez=sum(GP1[i,j]*prod(y[k]^PZ1.E[k,j] for k in 1:e1) for j in 1:n2)+PZ1.c[i]
        #@show(ez)
        push!(ranges1,enclose(a1,domain,sol))
    end
   # @show(ranges1)

    domain=IntervalBox(-1..1, f1)
    ranges2=[]
    for i in 1:m1
        a2(x)=sum(GP2[i,j]*prod(x[k]^PZ2.E[k,j] for k in 1:f1) for j in 1:m2)+PZ2.c[i]
        #println("deuxieme reste")
        #ez=sum(GP1[i,j]*prod(y[k]^PZ1.E[k,j] for k in 1:e1) for j in 1:n2)+PZ2.c[i]
        #@show(ez)
        push!(ranges2,enclose(a2,domain,sol))
    end
    #@show(ranges2)
    center=zeros(n1)
    
    for i in 1:n1
        mini=min(ranges1[i].lo,ranges2[i].lo)
        maxi=max(ranges1[i].hi,ranges2[i].hi)
        mid=Float64(1/2*(mini+maxi))
        #@show(mini,maxi,mid)
        #println(mini," ", maxi," ", mid)
        center[i]=mid
        z=zeros(Int64,n1+e1)
        z[i+e1]=1
        Enew=hcat(Enew,z)
        z=zeros(Float64,n1)
        z[i]=Float64(maxi-mini)/2
        #println(z[i])
        Gnew=hcat(Gnew,z)
    end

    return SimpleSparsePolynomialZonotope(center,Gnew,Enew)

end

function zonotopic_joinold(PZ1,PZ2,solver)
    """warning, it is required that PZ1 and PZ2 are defined over the same variables"""
    if solver!="bernstein"
        return zonotopic_joinbis(PZ1,PZ2,solver)
    end
    c1=PZ1.c 
    c2=PZ2.c 
    functionsPZ1=[]
    functionsPZ2=[]
    n1,n2=size(PZ1.G)
    m1,m2=size(PZ2.G)
    e1,e2=size(PZ1.E)
    f1,f2=size(PZ2.E)
    nbvar=max(e1,f1)
    #Gnew=zeros(Float64,n1,n2+n1)
    #Enew=zeros(Int64,n1+e1,e2+n1)
    Gnew=Matrix{Float64}(undef,n1,0) 
    Enew=Matrix{Int64}(undef,n1+e1,0) 
    #println(typeof(Gnew))
    for i in 1:n2
        for j in 1:m2
            if PZ1.E[:,i]==PZ2.E[:,j]
                #Enew[:,i]=vcat(PZ1.E[:,i],zeros(n1))#add the dimension 
                #Gnew[:,i]=[arg_min(PZ1.G[k,i],PZ2.G[k,j]) for k in 1:n1]
                Enew=hcat(Enew,vcat(PZ1.E[:,i],zeros(Int64,n1)))
                Gnew=hcat(Gnew,[Float64(arg_min(PZ1.G[k,i],PZ2.G[k,j])) for k in 1:n1])
            end
        end
    end

    domain=IntervalBox(-1..1, e1)
    ranges1=ranges_from_Bernsteincoeff(PZ1.G,PZ1.E,domain)
    for i in 1:n1 #WE HAVE TO ADD THE CENTER TO THE COMPUTED RANGES
        ranges1[i]=ranges1[i] .+ PZ1.c[i]

    end
    #println("ranges1: ",ranges1)
    
    domain=IntervalBox(-1..1, f1)
    ranges2=ranges_from_Bernsteincoeff(PZ2.G,PZ2.E,domain)
    for i in 1:n1
        ranges2[i]=ranges2[i] .+ PZ2.c[i]
    end
    #println("ranges2: ",ranges2)

    center=zeros(n1)
    
    for i in 1:n1
        mini=min(ranges1[i][1],ranges2[i][1])
        maxi=max(ranges1[i][2],ranges2[i][2])
        mid=Float64(1/2*(mini+maxi))
        println(mini," ", maxi," ", mid)
        center[i]=mid
        z=zeros(Int64,n1+e1)
        z[i+e1]=1
        Enew=hcat(Enew,z)
        z=zeros(Float64,n1)
        z[i]=Float64(maxi-mid-sum(abs(Gnew[i,l]) for l in 1:size(Gnew)[2]))
        println(z[i])
        Gnew=hcat(Gnew,z)
    end

    return SimpleSparsePolynomialZonotope(center,Gnew,Enew)
end

function zonotopic_join(PZ1,PZ2,solver,tolerance,maxdep)
    """warning, it is required that PZ1 and PZ2 are defined over the same variables"""
    if solver!="bernstein"
        return zonotopic_joinbis(PZ1,PZ2,solver,tolerance,maxdep)
    end
    c1=PZ1.c 
    c2=PZ2.c 
    n1,n2=size(PZ1.G)
    m1,m2=size(PZ2.G)
    e1,e2=size(PZ1.E)
    f1,f2=size(PZ2.E)
    if (e1!=f1)
        throw(DomainError(e1, "nombre de variables différent"))
    end

    Gnew=Matrix{Float64}(undef,n1,0) 
    Enew=Matrix{Int64}(undef,n1+e1,0) 

    interm=zeros(Float64,n1)
    test=1
    for i in 1:n2
        for j in 1:m2
            if PZ1.E[:,i]==PZ2.E[:,j]
                test=0
                #Enew[:,i]=vcat(PZ1.E[:,i],zeros(n1))#add the dimension 
                #Gnew[:,i]=[arg_min(PZ1.G[k,i],PZ2.G[k,j]) for k in 1:n1]
                Enew=hcat(Enew,vcat(PZ1.E[:,i],zeros(Int64,n1)))
                Gnew=hcat(Gnew,[Float64(arg_min(PZ1.G[k,i],PZ2.G[k,j])) for k in 1:n1])
            end
        end
    end
    nh=size(Enew)[2]
    GP1=fill_matrix_after_argmin(PZ1.G,PZ1.E,Gnew,Enew,n2,nh,e1)
    GP2=fill_matrix_after_argmin(PZ2.G,PZ2.E,Gnew,Enew,m2,nh,f1)

    domain=IntervalBox(-1..1, e1)
    ranges1=ranges_from_Bernsteincoeff(GP1,PZ1.E,domain)[1]
    for i in 1:n1 #WE HAVE TO ADD THE CENTER TO THE COMPUTED RANGES
        ranges1[i]=ranges1[i] .+ PZ1.c[i]
    end
    @show(ranges1)
    domain=IntervalBox(-1..1, f1)
    ranges2=ranges_from_Bernsteincoeff(GP2,PZ2.E,domain)[1]
    for i in 1:n1
        ranges2[i]=ranges2[i] .+ PZ2.c[i]
    end
    @show(ranges2)

    center=zeros(n1)
    
    for i in 1:n1
        mini=min(ranges1[i][1],ranges2[i][1])
        maxi=max(ranges1[i][2],ranges2[i][2])
        mid=Float64(1/2*(mini+maxi))
        @show(mini,maxi,mid)
        #println(mini," ", maxi," ", mid)
        center[i]=mid
        z=zeros(Int64,n1+e1)
        z[i+e1]=1
        Enew=hcat(Enew,z)
        z=zeros(Float64,n1)
        z[i]=Float64(maxi-mini)/2
        println(z[i])
        Gnew=hcat(Gnew,z)
    end

    return SimpleSparsePolynomialZonotope(center,Gnew,Enew)

end


function bernstein_zonotopic_join(PZ1,PZ2,field)
    """warning, it is required that PZ1 and PZ2 are defined over the same variables"""
    c1=PZ1.c 
    dim=length(c1)
    c2=PZ2.c 
    n1,n2=size(PZ1.G)
    m1,m2=size(PZ2.G)
    e1,e2=size(PZ1.E)
    f1,f2=size(PZ2.E)
    if (e1!=f1)
        throw(DomainError(e1, "nombre de variables différent"))
    end

    maxdegree=map(max,maxdeg(PZ1.E),maxdeg(PZ2.E))
    domain=IntervalBox(-1..1, e1)
    rangesbis1,coeffsbis1=ranges_from_Bernsteincoeff(PZ1.G,PZ1.E,domain)
    @show(rangesbis1)
    ranges1,coeffs1=ranges_from_Bernsteincoeff(PZ1.G,PZ1.E,domain,maxdeg=maxdegree)
    domain=IntervalBox(-1..1, f1)
    ranges2,coeffs2=ranges_from_Bernsteincoeff(PZ2.G,PZ2.E,domain,maxdeg=maxdegree)
    @show(ranges1,ranges2)
    for i in 1:dim #WE HAVE TO ADD THE CENTER TO THE COMPUTED RANGES 
        ranges1[i]=ranges1[i] .+ c1[i]
        ranges2[i]=ranges2[i] .+ c2[i]
    end
    @show(ranges1,ranges2)
    h=[]
    nbcoeffs=length(coeffs1[1])
    Anneau,(x)=polynomial_ring(field,e1+dim)
    @show(coeffs1)
    @show(coeffs2)
    
    res=[]
    for i in 1:dim
        
        A=min(ranges1[i][1],ranges2[i][1])
        B=max(ranges1[i][2],ranges2[i][2])

        println("coucou")
        coeffs1[i]=coeffs1[i] .+ c1[i]
        coeffs2[i]=coeffs2[i] .+ c2[i]
        argm=map(arg_min,coeffs1[i],coeffs2[i])
        @show(argm)
        push!(h,zeros(Float64,nbcoeffs))
        s=coeffs1[i]-argm
        t=coeffs2[i]-argm
        ms,Ms=minimum(s),maximum(s)
        mt,Mt=minimum(t),maximum(t)
        #on calcule notre polynome commun h
        for j in 1:length(argm) #on calcule notre polynome commun h
            if argm[j]>0
                h[i][j]=min(argm[j],B-Ms,B-Mt)
                #@show(Bshift,Ms,Mt)
            elseif argm[j]<0
                h[i][j]=max(argm[j],A-ms,A-mt)
                #@show(Ashift,ms,mt,Ashift-ms)
            else
                h[i][j]=0
            end
        
        end
        #on convertit les coeffs h[i] en un vrai polynome polyh
        @show(coeffs1[i],h[i],s)
        @show(coeffs2[i],h[i],t)
        println("polynomialfrombernsteincoeffs")
        hbis=polynomial_from_bernstein_coeffs(Anneau,maxdegree,domain,h[i])
        println("polynomialfrombernsteincoeffs ended")
        #hbis=copy_poly(h[i],Anneau)
        #on recalcule les polynomes f-h et g-h
        for j in 1:length(coeffs1[i])
            s[j]=coeffs1[i][j]-h[i][j]
        end
        for j in 1:length(coeffs2[i])
            t[j]=coeffs2[i][j]-h[i][j]
        end
        #s=s.-shift
        #t=t.-shift
        #t=coeffs2[i]-h .+ PZ2.c[i]
        mini=min(minimum(s),minimum(t))
        maxi=max(maximum(s),maximum(t))  
        #@show(mini,maxi)
        #on recalcule les polynomes f-h et g-h
        mid=Float64(1/2*(mini+maxi))
        #@show(mini,maxi,mid)
        #println(mini," ", maxi," ", mid)
        temp=mid+hbis+(Float64(maxi-mini)/2)*gens(Anneau)[e1+i]
        #@show(mid,hbis,(maxi-mini)/2,temp)
        push!(res,temp)
    end
    println("fin du join, reste à convertir")
    Ph=get_SSPZ_from_polynomials(res)
    return Ph

end


function fill_matrix_after_argmin(G,E,Gh,Eh,ngen,ngenh,nbvar)
    Gn=deepcopy(G)
    for i in 1:ngen
        for j in 1:ngenh
            if E[:,i]==Eh[1:nbvar,j]
                Gn[:,i]=Gn[:,i]-Gh[:,j]
            end
        end
    end
    return Gn
end

M=[ 1 2 3.0 3 3; 2 2 1 3 3;3 1 2 1 1; 1 1 1 2 2]
M[1:end,end]



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


function barycentre_union(PZ1::SimpleSparsePolynomialZonotope,PZ2::SimpleSparsePolynomialZonotope,field::Nemo.Field)
    """union barycentrique entre PZ1 et PZ2"""
    nb=size(genmat(PZ1))[1]
    nb_vars=max(size(expmat(PZ1))[1],size(expmat(PZ2))[1])
    #println("nombre de variables avant de join: ",nb_vars)
    Polynomes1=get_polynomials_from_SSPZ(PZ1,field)
    #println("nombre de var: ", length(gens(parent(Polynomes1[1]))))
    Polynomes2=get_polynomials_from_SSPZ(PZ2,field)
    #println("nombre de var: ", length(gens(parent(Polynomes2[1]))))
    Anneau,(x)=polynomial_ring(field,nb_vars+1)
    Polynomes=[Anneau(0) for i in 1:nb]
    for i in 1:nb
        Polynomes[i]=(1/2*gens(Anneau)[nb_vars+1] + 1/2)*copy_poly(Polynomes1[i],Anneau) + copy_poly(Polynomes2[i],Anneau)*(1/2 - 1/2*gens(Anneau)[nb_vars+1])
    end
    return get_SSPZ_from_polynomials(Polynomes)
end

function barycentre_union_simplifiee(PZ1::SimpleSparsePolynomialZonotope,PZ2::SimpleSparsePolynomialZonotope,field::Nemo.Field)
    """union barycentrique entre PZ1 et PZ2"""
    dim=size(genmat(PZ1))[1]
    nb_vars=max(size(expmat(PZ1))[1],size(expmat(PZ2))[1])
    #println("nombre de variables avant de join: ",nb_vars)
    Polynomes1=get_polynomials_from_SSPZ(PZ1,field)
    #println("nombre de var: ", length(gens(parent(Polynomes1[1]))))
    Polynomes2=get_polynomials_from_SSPZ(PZ2,field)
    #println("nombre de var: ", length(gens(parent(Polynomes2[1]))))
    Anneau,(x)=polynomial_ring(field,nb_vars+2)
    Polynomes=[Anneau(0) for i in 1:dim]
    for i in 1:dim
        Polynomes[i]=(1/2*gens(Anneau)[nb_vars+i] + 1/2)*copy_poly(Polynomes1[i],Anneau) + copy_poly(Polynomes2[i],Anneau)*(1/2 - 1/2*gens(Anneau)[nb_vars+i])
    end
    return get_SSPZ_from_polynomials(Polynomes)
end

function barycentric_join(SSPZ1,SSPZ2)
    c1=SSPZ1.c
    c2=SSPZ2.c
    m1,n1=size(SSPZ1.E)
    m2,n2=size(SSPZ2.E)
    m=max(m1,m2)
    center=ones(Int64,1,2)
    e1=zeros(Int64,1,n1)
    e_1=ones(Int64,1,n1)
    e2=zeros(Int64,1,n2)
    e_2=ones(Int64,1,n2)
    vecnew=hcat(center,e1,e_1,e2,e_2)
    t1=zeros(Int64,m,2)
    if m2>m1
        adjustment=zeros(Int64,m2-m1,n1)
        #E1=vcat(E1,adjustment)
        return remove_useless_terms!(
            SimpleSparsePolynomialZonotope(0.5*c1+0.5*c2,
            hcat(0.5*c1,-0.5*c2,0.5*SSPZ1.G,0.5*SSPZ1.G,0.5*SSPZ2.G,-0.5*SSPZ2.G),
            vcat(hcat(t1,vcat(SSPZ1.E,adjustment),vcat(SSPZ1.E,adjustment),SSPZ2.E,SSPZ2.E),vecnew))
        )
    elseif m1>m2
        adjustment=zeros(Int64,m1-m2,n2)
        #E2=vcat(E2,adjustment)
        return remove_useless_terms!(
            SimpleSparsePolynomialZonotope(0.5*c1+0.5*c2,
            hcat(0.5*c1,-0.5*c2,0.5*SSPZ1.G,0.5*SSPZ1.G,0.5*SSPZ2.G,-0.5*SSPZ2.G),
            vcat(hcat(t1,SSPZ1.E,SSPZ1.E,vcat(SSPZ2.E,adjustment),vcat(SSPZ2.E,adjustment)),vecnew))
        )
    end
    
    return remove_useless_terms!(
        SimpleSparsePolynomialZonotope(0.5*c1+0.5*c2,
        hcat(0.5*c1,-0.5*c2,0.5*SSPZ1.G,0.5*SSPZ1.G,0.5*SSPZ2.G,-0.5*SSPZ2.G),
        vcat(hcat(t1,SSPZ1.E,SSPZ1.E,SSPZ2.E,SSPZ2.E),vecnew))
    )

    #Redundant=SimpleSparsePolynomialZonotope(0.5*c1+0.5*c2,hcat(0.5*c1,-0.5*c2,0.5*G1,0.5*G1,0.5*G2,-0.5*G2),vcat(hcat(t1,E1,E1,E2,E2),vecnew))
    #res=remove_redundant_generators(Redundant) remove 2 fois pour enlever les colonnes nulles

end

function remove_useless_terms!(SSPZ)
    dim=size(SSPZ.G)[2]
    nb_vars,nb_gens=size(SSPZ.E)
    cnew = copy(SSPZ.c)
    visited_exps = Dict{Vector{Int},Int}()
    toremove=Int[]
    for i in 1:nb_gens
        gi=SSPZ.G[:,i]
        ei=@view SSPZ.E[:,i]
        if iszero(gi)
            push!(toremove,i)
        elseif iszero(ei)
            cnew+=gi
            push!(toremove,i)
        elseif haskey(visited_exps, ei) #repeated exponent
            #@show(SSPZ.G[:,visited_exps[ei]])
            SSPZ.G[:,visited_exps[ei]]+=gi
            #@show(SSPZ.G[:,visited_exps[ei]])
            #@show(gi)
            push!(toremove,i)
        else
            visited_exps[ei]=i
        end
    end
    for i in 1:nb_gens
        if iszero(@view SSPZ.G[:,i])
            push!(toremove,i)
        end
    end
    #@show(SSPZ.G)
    return SimpleSparsePolynomialZonotope(cnew,
    SSPZ.G[1:end,setdiff(1:end,toremove)],
    SSPZ.E[1:end,setdiff(1:end,toremove)])
end






