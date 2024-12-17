using Nemo
using LazySets


include("testinclusion_robrange.jl")
include("inclusiongeom_robustrange.jl")
include("conversions.jl")
include("joins.jl")
include("plotsample.jl")
include("polynomap.jl")
include("matlabmatrix.jl")
include("reduction.jl")



println("_____________________________________________")

function inclusion_lineaire_fonctionnelle()
    R=RealField()
    S,(x)=PolynomialRing(R,4)
    PZ=get_SSPZ_from_polynomials([x[1]+x[3],x[2]+x[4]])
    FPZ=get_SSPZ_from_polynomials([(x[1]+x[3]),(x[2]+x[4])])


    return inclusion_test(FPZ,PZ,1.1)
    #FPZ=poly_apply_on_SSPZ(PZ,[x[1]^2,x[2]^2],R)
    
end
#inclusion_lineaire_fonctionnelle()
println("_____________________________________________")

function inclusion_lineaire_fonctionnelle2()
    R=RealField()
    S,(x,y,s,t)=PolynomialRing(R,["x","y","s","t"])

    PZ=get_SSPZ_from_polynomials([x+t,y+s])
    FPZ=get_SSPZ_from_polynomials([x^2+x^3+s+t,y^2+s])

    return inclusion_test(FPZ,PZ,1.1)
    #FPZ=poly_apply_on_SSPZ(PZ,[x[1]^2,x[2]^2],R)
    
end

function inclusion_lineaire_geometrique()
    R=RealField()
    S,(x)=PolynomialRing(R,4)
    PZ=get_SSPZ_from_polynomials([x[1]+x[3],x[2]+x[4]])
    FPZ=get_SSPZ_from_polynomials([0.8*(x[1]+x[3]),0.9*(x[2]+x[4])])


    return geometrical_inclusion(FPZ,PZ,1.1)
    #FPZ=poly_apply_on_SSPZ(PZ,[x[1]^2,x[2]^2],R)
    
end
#inclusion_lineaire_geometrique()

println("_____________________________________________")

function inclusion_basique()
    R=RealField()
    S,(x)=PolynomialRing(R,4)
    PZ=get_SSPZ_from_polynomials([x[1]+x[3],x[2]+x[4]])
    FPZ=get_SSPZ_from_polynomials([0.1*(x[1]+x[3])^2,0.1*(x[2]+x[4])^2])

    PZ1=get_SSPZ_from_polynomials([x[1],x[2]+2])
    PZ2=get_SSPZ_from_polynomials([x[1]^2,x[2]*x[1]])
    BSI=barycentre_union_simplifiee(PZ1,PZ2,R)
    @show(get_polynomials_from_SSPZ(BSI,R))
    #return inclusion_test(FPZ,PZ,0.9)
    #FPZ=poly_apply_on_SSPZ(PZ,[x[1]^2,x[2]^2],R)
    
end

function inclusion_geometrique()
    R=RealField()
    S,(x)=PolynomialRing(R,2)
    PZ1=get_SSPZ_from_polynomials([x[1]^2+x[1]*x[2],x[1]*x[2]])
    PZ2=get_SSPZ_from_polynomials([x[1],x[2]])
end
#inclusion_basique()
#h