function compose(p,list_poly)
    return evaluate(p,list_poly)
end

function poly_apply_on_SSPZ(PZ::SimpleSparsePolynomialZonotope,list_poly,field::Field)
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

function iterate_polynomials_over_PZ(Polynomes,PZ::SimpleSparsePolynomialZonotope,nb_iter::Int64,borne_union::Int64,field::Field,choice;max_order::Int64,toreduce::Int64=200,maxdegree::Int64=50,scale_factor::Float64=1.1,power::Int64=1)
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