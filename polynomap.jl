function compose(p,list_poly)
    return Nemo.evaluate(p,list_poly)
end

function applyfilter(PZ::SimpleSparsePolynomialZonotope,A::Matrix{Float64},B::Vector{Float64},m::Float64,M::Float64) #VERIFIEE
    nbvars=size(PZ.E)[1] 
    ngen=size(PZ.G)[2]
    cnew=[PZ.c[2], A[2,1]*PZ.c[1]+A[2,2]*PZ.c[2]+(B[1]+B[2]+B[3])*(m+M)/2]
    exponew=zeros(Int64,nbvars+1)
    exponew[nbvars+1]=1
    genew=[0.0 , B[3]*(M-m)/2]
    vectores=zeros(Float64,ngen)
    vectores[end],vectores[end-1]=B[2]*(M-m)/2,B[1]*(M-m)/2
    Enew=hcat(vcat(PZ.E,zeros(Int64,1,ngen)),exponew)

    Gnew=hcat(vcat(reshape(PZ.G[2,:],1,:),reshape(A[2,1]*PZ.G[1,:]+A[2,2]*PZ.G[2,:]+vectores,1,:)),genew)
    return SimpleSparsePolynomialZonotope(cnew,Gnew,Enew)

end

function powerfilter(PZ::SimpleSparsePolynomialZonotope,A::Matrix{Float64},B::Vector{Float64},m::Float64,M::Float64,power::Int64) #VERIFIEE SUR EXEMPLE ERIC SYLVIE
    for p in 1:power
        PZ=applyfilter(PZ,A,B,m,M)
    end
    return PZ
end
st=SimpleSparsePolynomialZonotope([0.0, 0.0],[0.0 0; 0 0],[1 0; 0 1])
HU=applyfilter(st,[0 1; -0.7 1.4],[1.1, -1.3, 0.7],0.0,1.0) 
@time sparse_poly_zono_to_dynamic(HU)[1]
HU.G
HU.c
HU.E
@time sparse_poly_zono_to_dynamic(HU)[1]
PHU=applyfilter(HU,[0 1; -0.7 1.4],[1.1, -1.3, 0.7],0.0,1.0)
PHU1=powerfilter(st,[0 1; -0.7 1.4],[1.1, -1.3, 0.7],0.0,1.0,98)
PHU1.c
PHU1.G
@time sparse_poly_zono_to_dynamic(PHU)[1]

function align_exponent_matrices(PZ1, PZ2)
    e1, _ = size(PZ1.E)
    f1, _ = size(PZ2.E)
    if e1 < f1
        padding = zeros(Int64, f1 - e1, size(PZ1.E, 2))
        PZ1 = SimpleSparsePolynomialZonotope(PZ1.c, PZ1.G, vcat(PZ1.E, padding))
    elseif f1 < e1
        padding = zeros(Int64, e1 - f1, size(PZ2.E, 2))
        PZ2 = SimpleSparsePolynomialZonotope(PZ2.c, PZ2.G, vcat(PZ2.E, padding))
    end
    return PZ1, PZ2
end

function inclusion_after_filter(FPZ,PZ,join,tolerance,depth)

    if join=="barycentric"
        a1(v)=sum(FPZ.G[2,j]*prod(v[k]^FPZ.E[k,j] for k in 1:2) for j in 1:size(FPZ.G)[2])+FPZ.c[2]
        b1(v)=sum(PZ.G[2,j]*prod(v[k]^PZ.E[k,j] for k in 1:2) for j in 1:size(PZ.G)[2])+PZ.c[2]
        domain=IntervalBox(-1..1, 2)
        a=enclose(a1,domain,BranchAndBoundEnclosure(tol=tolerance,maxdepth=depth))
        b=enclose(b1,domain,BranchAndBoundEnclosure(tol=tolerance,maxdepth=depth))
        @show(a,b)
        if a.lo>=b.lo +2*tolerance && a.hi+2*tolerance<=b.hi
            return 1
        end
    else
        a=interval(FPZ.c[2] - sum(abs(FPZ.G[2,j]) for j in 1:size(FPZ.E)[2]), FPZ.c[2] + sum(abs(FPZ.G[2,j]) for j in 1:size(FPZ.E)[2]))
        b=interval(PZ.c[2] - sum(abs(PZ.G[2,j]) for j in 1:size(PZ.E)[2]), PZ.c[2] + sum(abs(PZ.G[2,j]) for j in 1:size(PZ.E)[2]))
        #@show(sparse_poly_zono_to_dynamic(FPZ)[1][2])
        @show(a,b)
        #@show(sparse_poly_zono_to_dynamic(PZ)[1][2])
        if a.lo>=b.lo  && a.hi<=b.hi
            return 1
        end
    end
    return "false"
end

function iterate_filter_over_PZ(PZ::SimpleSparsePolynomialZonotope,nbiter::Int64,A::Matrix{Float64},B::Vector{Float64},m::Float64,M::Float64,power::Int64,join::String; solver::String="NaturalEnclosure",tolerance::Float64=1e-3,depth::Int64=10)
    i=1
    while i<=nbiter
        if i>1
            FPZ=applyfilter(PZ,A,B,m,M)
            @show(sparse_poly_zono_to_dynamic(PZ)[1])
            println("inclusion test",i)
            if inclusion_after_filter(FPZ,PZ,join,tolerance,depth)==1
                println("INCLUSION!")
                return PZ
            end
        end
        FPZ=powerfilter(PZ,A,B,m,M,power)
        if join=="barycentric"
            PZ=barycentric_join(PZ,FPZ)
        else
            PZ,FPZ=align_exponent_matrices(PZ, FPZ)
            PZ=zonotopic_join(PZ,FPZ,solver,tolerance,depth)
            PZ=remove_redundant_generators(PZ)
        end
        i=i+1
    end
    return sparse_poly_zono_to_dynamic(PZ)[1],PZ
end

"""SAVE THE RESULT FOR I=5 and zonotopic join: 24 iterations"""

@time iterate_filter_over_PZ(st,2,[0 1; -0.7 1.4],[1.1, -1.3, 0.7],0.0,1.0,5,"zono",solver="NaturalEnclosure",tolerance=1e-2,depth=10)
"""SAVE THE RESULT FOR I=16"""

remove_redundant_generators(h).E
sparse_poly_zono_to_dynamic(remove_redundant_generators(h))[1][1]
e[1]
@time iterate_filter_over_PZ(st,15,[0 1; -0.7 1.4],[1.1, -1.3, 0.7],0.0,1.0,16,"zono",solver="NaturalEnclosure",tolerance=1e-2,depth=10)


function poly_apply_on_SSPZ(PZ::SimpleSparsePolynomialZonotope,list_poly,field::Nemo.Field)
    nb_vars=size(expmat(PZ))[1]#nombre de variables est le nombre de lignes de la matrice des exposants
    composit=[]
    #anneau,(x)=polynomial_ring(field,nb_vars)
    Poly_fromPZ=get_polynomials_from_SSPZ(PZ,field)
    #anneau=parent(Poly_fromPZ[1])
    #println("On va entrer dans la composition de l'application")
    for p in list_poly
        res=compose(p,Poly_fromPZ)
        push!(composit,res)
    end
    #println("Onest sorti de la composition de l'application")
    #println("a l'interieur de la compo: ",length(composit[1]))
    #res=get_SSPZ_from_polynomials(composit)
    #println("apres changement de représentation: ",size(res.E)[2])
    return get_SSPZ_from_polynomials(composit)
end

function poly_apply_on_SSPZdynamic(PZ,list_poly)
    #println("Bonne transformation polynomiale")
    polyPZ=sparse_poly_zono_to_dynamic(PZ)[1]
    #@show(polyPZ)
    Transformed=polynomial_transformation(polyPZ,list_poly)
    #@show(Transformed)
    return dynamic_to_sparse_poly_zono(Transformed)
end

function iterate_polynomials_over_PZ(Polynomes,PZ::SimpleSparsePolynomialZonotope,nb_iter::Int64,borne_union::Int64,field::Nemo.Field,choice;max_order::Int64,toreduce::Int64=200,maxdegree::Int64=50,scale_factor::Float64=1.1,power::Int64=1,inclusiontest::Int64=1,solver="bernstein",nbinitialvariables::Int64=2,linearsystem::Bool=false,tolerance::Float64=1e-3,maxdepth::Int64=10)
    """il faudrait quand même trouver un moyen efficace de tester l'inclusion entre polynomial zonotopes"""
    #println("on entre dans l'itération")
    i=0
    nb_reduc=0
    liste=[PZ]
    
    while i<nb_iter
        println("itération numéro",i)
        fPZ=liste[end]
        PZ_previous=fPZ
        
        for p in 1:power 
            println("avant appli")
            fPZ=poly_apply_on_SSPZ(fPZ,Polynomes,field)#pas de problème d'aliasing entre les arrays ici
            println("apres appli")
        end
        #@show(get_polynomials_from_SSPZ(fPZ,R))
        if i>=borne_union+1 && inclusiontest==1
            inclusion=inclusion_test(fPZ,PZ,1.5,nbinitialvariables,linearsystem)
            @show(inclusion)
            if inclusion!="false"
                println("INCLUSION!")
                #push!(liste,fPZ)
                return liste
            end
        end

        #=if nb_reduc==0 && polytest1!=polytest2
            println("LOUPE")
            return liste
        end=#
        if i>= borne_union
            println("entre join")
            if choice=="bernstein"
                println("ON FAIT BIEN LE JOIN DE BERNSTEIN")
                PZ=bernstein_zonotopic_join(PZ_previous,fPZ,field)
            elseif choice=="zono"
                
                PZ=zonotopic_join(PZ_previous,fPZ,solver,tolerance,maxdepth)
                PZ=remove_unused_variables(PZ)
            else
                PZ=barycentric_join(PZ_previous,fPZ)
                #PZ=barycentre_union_simplifiee(PZ_previous,fPZ,field)
                #PZ=remove_useless_terms!(PZ) PAS BESOIN PUISQUE CEST DEJA DANS LE JOIN
            end
            println("sort join")
        else
            PZ=fPZ
        end

        #println("number of terms/monomials after join ",size(expmat(PZ))[2])
        if size(PZ.G)[2]>=max_order
            PZ=Simple_reduce_order(PZ,max_order)
            nb_reduc+=1
        end

        i+=1
        push!(liste,PZ)
    end
    println("nombre de réductions, ",nb_reduc)
    return liste
end

function iterate_oned(Polynomes,PZ::SimpleSparsePolynomialZonotope,nb_iter::Int64,borne_union::Int64,field::Nemo.Field,choice;max_order::Int64,toreduce::Int64=200,maxdegree::Int64=50,power::Int64=1,inclusiontest::Int64=1,solver="NaturalEnclosure",tolerance::Float64=1e-2,maxdepth::Int64=10)
    i=0
    nb_reduc=0
    liste=[PZ]
    
    while i<nb_iter
        println("itération numéro",i)
        fPZ=liste[end]
        PZ_previous=fPZ
        for p in 1:power 
            println("avant appli")
            fPZ=poly_apply_on_SSPZ(fPZ,Polynomes,field)
            println("apres appli")
        end
        if i>=borne_union+1 && inclusiontest==1
            inclusion=inclusion_after_filter(fPZ,PZ,choice,tolerance,maxdepth)
            @show(inclusion)
            if inclusion!="false"
                println("INCLUSION!")
                return liste
            end
        end
        if i>= borne_union
            println("entre join")
            if choice=="zono"
                PZ=zonotopic_join(PZ_previous,fPZ,solver,tolerance,maxdepth)
                PZ=remove_unused_variables(PZ)
            else
                PZ=barycentric_join(PZ_previous,fPZ)
            end
            println("sort join")
        else
            PZ=fPZ
        end
        if size(PZ.G)[2]>=max_order
            PZ=Simple_reduce_order(PZ,max_order)
            nb_reduc+=1
        end
        i+=1
        push!(liste,PZ)
    end
    println("nombre de réductions, ",nb_reduc)
    return liste
end