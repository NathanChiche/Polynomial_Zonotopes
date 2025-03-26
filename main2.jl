using Nemo
using LazySets
using Plots
using Random
using Dates
using TimerOutputs
using Profile
using PProf
using ProfileView
using RangeEnclosures
using IntervalArithmetic
using ForwardDiff
using IntervalOptimisation
using LinearAlgebra
using DynamicPolynomials
using BernsteinExpansions
using TaylorModels
#using SumOfSquares
#using AffineArithmetic
#using ProfileVega


include("conversions.jl")
include("joins.jl")
include("plotsample.jl")
include("polynomap.jl")
include("matlabmatrix.jl")
include("reduction.jl")
#include("bernstein.jl")
include("testinclusion_robrange.jl")

R=RealField()
A=SimpleSparsePolynomialZonotope([0.0, 0.0],[1 1.0 0; 0 1 1],[1 2 1; 0 0 1])
get_polynomials_from_SSPZ(A,R)
get_polynomials_from_SSPZ(B,R)
B=SimpleSparsePolynomialZonotope([1.0, 0.0],[0 1.0 0;1 3 2],[1 2 1; 0 0 1])
#plot_multiple([A,B],R,"Documents/julia/plots_julia/ajoin",nbpoints=80000)
C1=zonotopic_joinbisold(A,B,"BranchAndBoundEnclosure")
#plot_multiple([C1,A,B],R,"Documents/julia/plots_julia/C1joinbisold",nbpoints=800000)
C2=zonotopic_joinbis(A,B,"BranchAndBoundEnclosure")
#plot_multiple([C2,A,B],R,"Documents/julia/plots_julia/C2joinbis",nbpoints=800000)
get_polynomials_from_SSPZ(C2,R)

dei(x)=2*x[1]^2+x[1]*x[2]+x[1]
dozk=IntervalBox(-1..1,2)
der=enclose(dei,dozk,BranchAndBoundEnclosure()).lo
der2=enclose(dei,dozk,BranchAndBoundEnclosure()).hi
(der2-der)/2


R=RealField()
Annea,(x,y)=polynomial_ring(R,["x","y"])

function lineaireconverge()
    R=RealField()
    S,(x,y)=polynomial_ring(R,["x","y"])
    #f1=0.2*x - 0.1*y
    #f2=0.05*x + 0.3*y
    f1=0.4*x+0.1*y
    f2=0.2*y+0.2*x
    Pstart=get_SSPZ_from_polynomials([x+1,y+1])
    fin=iterate_polynomials_over_PZ([f1,f2],Pstart,4,1,R,"zono",max_order=200000,inclusiontest=0,solver="NaturalEnclosure")
    #popfirst!(fin)
    #fin=reverse(fin)
    #plot_multiple(fin,R,"Documents/julia/plots_julia/lineaireconverge1_zono_4iter_join81_test",nbpoints=15000)
    #fini=fin[end]
    #next=poly_apply_on_SSPZ(fini,[f1,f2],R)
    return fin

end
@time fin=lineaireconverge()
fin[1].G
@time lin,next=lineaireconverge()
lin.E
lin.G
plot(lin)
plot_multiplefinnext([lin,next],R,"Documents/julia/plots_julia/lineaire1_zono_4iter_join1_test",nbpoints=15000)

function lineaireconverge2()
    R=RealField()
    S,(x,y)=polynomial_ring(R,["x","y"])
    #f1=0.2*x - 0.1*y
    #f2=0.05*x + 0.3*y
    f1=0.4*x+0.2*y
    f2=0.2*y+0.2*x
    Pstart=get_SSPZ_from_polynomials([x+0.18,y-1.8])
    fin=iterate_polynomials_over_PZ([f1,f2],Pstart,30,1,R,"zono",max_order=200000,inclusiontest=1,solver="NaturalEnclosure",power=2,linearsystem=true)
    #popfirst!(fin)
    #fin=reverse(fin)
    #plot_multiple(fin,R,"Documents/julia/plots_julia/lineaireconverge1_zono_4iter_join81_test",nbpoints=15000)
    #fini=fin[end]
    #next=poly_apply_on_SSPZ(fini,[f1,f2],R)
    return fin
end
@time res=lineaireconverge2()
res[7].E

function lineaireconvergebary2()
    R=RealField()
    S,(x,y)=polynomial_ring(R,["x","y"])
    #f1=0.2*x - 0.1*y
    #f2=0.05*x + 0.3*y
    f1=0.4*x+0.2*y
    f2=0.2*y+0.2*x
    Pstart=get_SSPZ_from_polynomials([x,y])
    fin=iterate_polynomials_over_PZ([f1,f2],Pstart,7,1,R,"bary",max_order=200000,inclusiontest=0,solver="NaturalEnclosure",power=2,linearsystem=false)
    popfirst!(fin)
    fin=reverse(fin)
    for i in 1:length(fin)
        @show(size(fin[i].E)[1])
    end
    #plot_multiple(fin,R,"Documents/julia/plots_julia/lineaireconverge2_bary_4iter_join81_testreversed",nbpoints=10000)
    #fini=fin[end]
    next=poly_apply_on_SSPZ(fin[1],[f1,f2],R)
    MatlabMatrix(fin[1],"Documents/julia/Traduct_Matlab/finlineaire2bary.txt","finlinbary")
    MatlabMatrix(next,"Documents/julia/Traduct_Matlab/nextlineaire2bary.txt","nextlinbary")
    #plot_multiple([fin[1],next],R,"Documents/julia/plots_julia/finnextlineaireconverge2_bary_4iter_join1_testreversed",nbpoints=20000)
    return fin,next
end
@time res=lineaireconvergebary2()





function lineaireconvergegp()
    R=RealField()
    S,(x,y)=polynomial_ring(R,["x","y"])
    f1=-0.5488*x + 0.2156*y
    f2=  -0.15092*x - -0.24696*y
    Pstart=get_SSPZ_from_polynomials([x,y])
    fin=iterate_polynomials_over_PZ([f1,f2],Pstart,4,1,R,"zono",max_order=200000,power=1,solver="NaturalEnclosure")
    popfirst!(fin)

    #fpp=[fin[length(fin)],fin[length(fin)-1]]
    #plot_multiple(fin,R,"Documents/julia/plots_julia/lineairegp_zono_3iter_nojoin",nbpoints=15000)
    return fin
end
@time res=lineaireconvergegp()
res[end].G


function henon()
    R=RealField()
    S,(x,y)=polynomial_ring(R,["x","y"])
    a=0.2
    b=0.3
    h1=y+1-a*x^2
    h2=b*y
    Pstart=get_SSPZ_from_polynomials([x; y])
    fin=iterate_polynomials_over_PZ([h1,h2],Pstart,18,1,R,"zono",max_order=200000,power=1,inclusiontest=1,solver="BranchAndBoundEnclosure")
    #popfirst!(fin)
    #fini=fin[end]
    #next=poly_apply_on_SSPZ(fini,[h1,h2],R)
    #return fini,next
    #plot_multiple(fpp,R,"Documents/julia/plots_julia/fpp_henon_zono_4iter_join_2",nbpoints=40000)
    #@show(get_polynomials_from_SSPZ(fin[1],R))
    return fin
end
@time reshenon1=henon()
plot_multiplefinnext(reshenon1,R,"Documents/julia/plots_julia/henon1_zono_5iter_join2",nbpoints=400000)



function polynomialnonconverge()
    R=RealField()
    S,(x,y)=polynomial_ring(R,["x","y"])

    p1=1/4*(x+x^2)
    p2=1/4*(y+x)
    Pstart=get_SSPZ_from_polynomials([x; y])
    fin=iterate_polynomials_over_PZ([p1,p2],Pstart,5,1,R,"zono",max_order=20000000,power=1,inclusiontest=0,solver="BranchAndBoundEnclosure")
    popfirst!(fin)
    fini=fin[end]
    next=poly_apply_on_SSPZ(fini,[p1,p2],R)
    #return fin
    return fini,next
end
@time respolynome1=polynomialnonconverge()
#plot_multiple(respolynome1,R,"Documents/julia/plots_julia/polynome1_zono_5iter_join7decal",nbpoints=200000)
plot_multiplefinnext(respolynome1,R,"Documents/julia/plots_julia/polynome1_bary_5iter_join2decal",nbpoints=400000)

@VSCodeServer.profview res=polynomialnonconverge()
#@VSCodeServer.profview_allocs g=goal()






function exemplejoinzono()
    R=RealField()
    S,(x,y)=polynomial_ring(R,["x","y"])

    A=get_SSPZ_from_polynomials([x+y^2;x*y+y^3])
    B=get_SSPZ_from_polynomials([x+x^2+0.5*y^2;x*y+1])
    J=zonotopic_join(A,B,"BranchAndBoundEnclosure")
    @show(get_polynomials_from_SSPZ(J,R))
    #plot_multiple([J,A,B],R,"Documents/julia/plots_julia/testjoinzonopapier",nbpoints=60000)
    #MatlabMatrix(A,"Documents/julia/Traduct_Matlab/poisson1.txt","pois1")
    #MatlabMatrix(B,"Documents/julia/Traduct_Matlab/poisson2.txt","pois2")
    #MatlabMatrix(J,"Documents/julia/Traduct_Matlab/poissonj.txt","poisj")
    return A,B,J
end

exemplejoinzono()

