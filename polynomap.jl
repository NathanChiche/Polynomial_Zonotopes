function compose(p,list_poly)
    return evaluate(p,list_poly)
end

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

function iterate_polynomials_over_PZ(Polynomes,PZ::SimpleSparsePolynomialZonotope,nb_iter::Int64,borne_union::Int64,field::Nemo.Field,choice;max_order::Int64,toreduce::Int64=200,maxdegree::Int64=50,scale_factor::Float64=1.1,power::Int64=1,inclusiontest::Int64=1,solver="bernstein")
    """il faudrait quand même trouver un moyen efficace de tester l'inclusion entre polynomial zonotopes"""
    #println("on entre dans l'itération")
    i=0
    nb_reduc=0
    liste=[PZ]
    
    while i<nb_iter
        println("itération numéro",i)
        #println("nb variables PZ: ",size(PZ.E)[1])
        fPZ=liste[end]
        PZ_previous=fPZ
        
        for p in 1:power 
            fPZ=poly_apply_on_SSPZ(fPZ,Polynomes,field)#pas de problème d'aliasing entre les arrays ici
        end
        if i>=borne_union+1 && inclusiontest==1
            inclusion=inclusion_test(fPZ,PZ,1.5)
            @show(inclusion)
            if inclusion!="false"
                println("INCLUSION!")
                return fPZ
            end
        end
        #println("number of terms after the composition: ",size(fPZ.E)[2])
        #println("nb variables PZ_interm: ",size(fPZ.E)[1])
        #println("nombre de termes avant la réduction/après fonctionnelle: ",size(fPZ.E)[2])
        #polytest1=get_polynomials_from_SSPZ(fPZ,field)
        
        #polytest2=get_polynomials_from_SSPZ(fPZ,field)

        #=if nb_reduc==0 && polytest1!=polytest2
            println("LOUPE")
            return liste
        end=#

        if i>= borne_union
            if choice=="bernstein"
                println("ON FAIT BIEN LE JOIN DE BERNSTEIN")
                PZ=bernstein_zonotopic_join(PZ_previous,fPZ,field)
            elseif choice=="zono"
                PZ=zonotopic_join(PZ_previous,fPZ,solver)
                PZ=remove_unused_variables(PZ)
            else
                #PZ=barycentric_join(PZ_previous,fPZ)
                PZ=barycentre_union_simplifiee(PZ_previous,fPZ,field)
                #PZ=remove_useless_terms!(PZ) PAS BESOIN PUISQUE CEST DEJA DANS LE JOIN
            end
        else
            PZ=fPZ
        end
        #println("number of terms/monomials after join ",size(expmat(PZ))[2])
        if size(PZ.G)[2]>=max_order
            PZ=Simple_reduce_order(PZ,max_order)
            nb_reduc+=1
        end
        #println("number of terms after reduction: ",size(PZ.E)[2])
        i+=1
        push!(liste,PZ)
    end
    #println("nombre de réductions, ",nb_reduc)
    return liste
end