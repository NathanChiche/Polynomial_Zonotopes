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

function identite(nbiter,dx,dy,powerf)
    R=RealField()
    S,(x,y)=polynomial_ring(R,["x","y"])
    #f1=0.2*x - 0.1*y
    #f2=0.05*x + 0.3*y
    f1=x
    f2=y
    Pstart=get_SSPZ_from_polynomials([x+dx,y+dy])
    fin=iterate_polynomials_over_PZ([f1,f2],Pstart,nbiter,1,R,"zono",max_order=200000,inclusiontest=1,solver="NaturalEnclosure",power=powerf,linearsystem=true)
    #popfirst!(fin)
    #fin=reverse(fin)
    #plot_multiple(fin,R,"Documents/julia/plots_julia/lineaireconverge1_zono_4iter_join81_test",nbpoints=15000)
    #fini=fin[end]
    #next=poly_apply_on_SSPZ(fini,[f1,f2],R)
    #return fin
    println("finito")

end
@time identite(20,0,0,1)

function lineaireconverge(nbiter,dx,dy,powerf)
    R=RealField()
    #S,(x,y)=polynomial_ring(R,["x","y"])
    #f1=0.2*x - 0.1*y
    #f2=0.05*x + 0.3*y
    @polyvar x1 x2
    f1=0.3*x1+0.1*x2
    f2=0.4*x2
    #Pstart=get_SSPZ_from_polynomials([x+dx,y+dy])
    Pstart=dynamic_to_sparse_poly_zono([x1+dx; x2+dx])
    fin=iterate_polynomials_over_PZ([f1,f2],Pstart,nbiter,1,R,"zono",max_order=200000,inclusiontest=1,solver="NaturalEnclosure",power=powerf,linearsystem=true)
    #popfirst!(fin)
    #fin=reverse(fin)
    #plot_multiple(fin,R,"Documents/julia/plots_julia/lineaireconverge1_zono_4iter_join81_test",nbpoints=15000)
    #fini=fin[end]
    #next=poly_apply_on_SSPZ(fini,[f1,f2],R)
    #return fin
    println("finito")

end
@time lineaireconverge(4,0,0,1)
24/100
fin[1].G
@time lin,next=lineaireconverge()
lin.E
lin.G
plot(lin)
plot_multiplefinnext([lin,next],R,"Documents/julia/plots_julia/lineaire1_zono_4iter_join1_test",nbpoints=15000)

function lineaireconverge2(nbiter)
    R=RealField()
    S,(x,y)=polynomial_ring(R,["x","y"])
    #f1=0.2*x - 0.1*y
    #f2=0.05*x + 0.3*y
    f1=0.4*x+0.2*y
    f2=0.2*y+0.2*x
    Pstart=get_SSPZ_from_polynomials([x+0.17,y-1.8])
    fin=iterate_polynomials_over_PZ([f1,f2],Pstart,nbiter,1,R,"zono",max_order=200000,inclusiontest=1,solver="NaturalEnclosure",power=2,linearsystem=true)
    #popfirst!(fin)
    #fin=reverse(fin)
    #plot_multiple(fin,R,"Documents/julia/plots_julia/lineaireconverge1_zono_4iter_join81_test",nbpoints=15000)
    #fini=fin[end]
    #next=poly_apply_on_SSPZ(fini,[f1,f2],R)
    return fin
end
@time res=lineaireconverge2(6)
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

1/2^4

function inversesquareroot(nbiter,dx,dy,powerf)
    R=RealField()
    S,(x,y)=polynomial_ring(R,["x","y"])
    p1=x+x*(0.5*(1-y*x^2)+0.375*(1-y*x^2)^2)
    p2=y
    Pstart=get_SSPZ_from_polynomials([0*x+dx; 2*y+dy])
    fin=iterate_polynomials_over_PZ([p1,p2],Pstart,nbiter,7,R,"zono",max_order=20000000,power=powerf,inclusiontest=0,solver="BranchAndBoundEnclosure",tolerance=1e-4,maxdepth=12)
    popfirst!(fin)
    fini=fin[end]
    next=poly_apply_on_SSPZ(fini,[p1,p2],R)
    #return fin
    return fini,next
end
respolynome1=inversesquareroot(4,1/2^4,18,1)
rangesfrompolynomialzonotope(respolynome1[1],10e-3,10)
#plot_multiple(respolynome1,R,"Documents/julia/plots_julia/polynome1_zono_5iter_join7decal",nbpoints=200000)
plot_multiplefinnext(respolynome1,R,"Documents/julia/plots_julia/polynome1_bary_5iter_join2decal",nbpoints=400000)

#@VSCodeServer.profview res=polynomialnonconverge()
#@VSCodeServer.profview_allocs g=goal()

1/0.273539
1/0.111227
1/0.223602
1/0.249969

function rangesfrompolynomialzonotope(PZ,tolerance,maxde)
    d=length(PZ.c)
    nbvars=size(PZ.E)[1]
    ngener=size(PZ.G)[2]
    domain=IntervalBox(-1..1, nbvars)
    for i in 1:d
        a1(v)=sum(PZ.G[i,j]*prod(v[k]^PZ.E[k,j] for k in 1:nbvars) for j in 1:ngener)+PZ.c[i]
        println(enclose(a1,domain,BranchAndBoundEnclosure(tol=tolerance,maxdepth=maxde)))
        #println(enclose(a1,domain,MooreSkelboeEnclosure()))
    end
end

function exemplejoinzono()
    R=RealField()
    S,(x,y)=polynomial_ring(R,["x","y"])

    A=get_SSPZ_from_polynomials([x+y^2;x*y+y^3])
    Testrange=get_SSPZ_from_polynomials([x+x^2;x*y+y^2+y^3])
    @show(rangesfrompolynomialzonotope(Testrange,1e-4,14))
    B=get_SSPZ_from_polynomials([x+x^2+0.5*y^2;x*y+1])
    J=zonotopic_join(A,B,"BranchAndBoundEnclosure")
    #@show(get_polynomials_from_SSPZ(J,R))
    #plot_multiple([J,A,B],R,"Documents/julia/plots_julia/testjoinzonopapier",nbpoints=60000)
    #MatlabMatrix(A,"Documents/julia/Traduct_Matlab/poisson1.txt","pois1")
    #MatlabMatrix(B,"Documents/julia/Traduct_Matlab/poisson2.txt","pois2")
    #MatlabMatrix(J,"Documents/julia/Traduct_Matlab/poissonj.txt","poisj")
    return A,B,J
end

exemplejoinzono()

