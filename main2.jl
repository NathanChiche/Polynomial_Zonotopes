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
#using TaylorModels
#using SumOfSquares
#using AffineArithmetic
#using ProfileVega


include("conversions.jl")
include("joins.jl")
include("plotsample.jl")
include("polynomap.jl")
include("matlabmatrix.jl")
include("reduction.jl")
include("bernstein.jl")
include("testinclusion_robrange.jl")

R=RealField()
Annea,(x,y)=polynomial_ring(R,["x","y"])

function lineairenonconverge()
    R=RealField()
    S,(x,y)=polynomial_ring(R,["x","y"])
    f1=-0.032*x -0.3968*y
    f2=0.5208*x + 0.712*y
    Pstart=get_SSPZ_from_polynomials([3*x,3*y])
    fin=iterate_polynomials_over_PZ([f1,f2],Pstart,3,1,R,"bernstein",max_order=200000,inclusiontest=1,solver="bernstein")
    popfirst!(fin)
    fini=fin[end]
    #@show(get_polynomials_from_SSPZ(fin[1],R))
    #@show(get_polynomials_from_SSPZ(fin[2],R)) 
    next=poly_apply_on_SSPZ(fini,[f1,f2],R)
    #fin=reverse(fin)
    #fpp=[fin[length(fin)],fin[length(fin)-1]]
    #plot_multiple([fini,next],R,"Documents/julia/plots_julia/lineaire1_zono_2iter_join1decal",nbpoints=3000000)
    return fini,next

end
@time lin,next=lineairenonconverge()
 plot(lin,next)


#plot_multiple(lin,R,"Documents/julia/plots_julia/lineaire1_zono_5iter_join7decal",nbpoints=800000)
plot([lin,next],nsdiv=18)



function lineaireconverge()
    R=RealField()
    S,(x,y)=polynomial_ring(R,["x","y"])
    f1=-0.32*x + 0.32*y
    f2=-0.42*x - 0.92*y
    Pstart=get_SSPZ_from_polynomials([2*x-5,1.5*y+1])
    fin=iterate_polynomials_over_PZ([f1,f2],Pstart,6,7,R,"zono",max_order=200000,power=1,solver="bernstein")
    popfirst!(fin)

    #fpp=[fin[length(fin)],fin[length(fin)-1]]
    plot_multiple(fin,R,"Documents/julia/plots_julia/lineaire1_zono_6iter_nojoin",nbpoints=40000)

end

function henon()
    R=RealField()
    S,(x,y)=polynomial_ring(R,["x","y"])
    a=0.2
    b=1
    h1=y+1-a*x^2
    h2=b*y
    Pstart=get_SSPZ_from_polynomials([x; y])
    fin=iterate_polynomials_over_PZ([h1,h2],Pstart,4,2,R,"zono",max_order=2000,power=1,solver="bernstein")
    
    popfirst!(fin)
    fpp=[fin[length(fin)],fin[length(fin)-1]]
    plot_multiple(fpp,R,"Documents/julia/plots_julia/fpp_henon_zono_4iter_join_2",nbpoints=40000)
    #@show(get_polynomials_from_SSPZ(fin[1],R))
    return fin
end



function polynomialnonconverge()
    R=RealField()
    S,(x,y)=polynomial_ring(R,["x","y"])

    p1=1/4*(x+x^2)
    p2=1/4*(y+x)
    Pstart=get_SSPZ_from_polynomials([x; y])
    fin=iterate_polynomials_over_PZ([p1,p2],Pstart,8,1,R,"zono",max_order=2000000000,power=1,solver="bernstein")
    
    return fin
end

res=polynomialnonconverge()

@VSCodeServer.profview_allocs res=polynomialnonconverge()
#@VSCodeServer.profview_allocs g=goal()


function exemplejoinzono()
    R=RealField()
    S,(x,y)=polynomial_ring(R,["x","y"])

    p1=1/4*(x+x^2)
    p2=1/4*(y+x)
    Pstart=get_SSPZ_from_polynomials([x; y])
    fPstart=poly_apply_on_SSPZ(Pstart,[p1,p2],R)
    S,(x,y,s,t)=polynomial_ring(R,["x","y","s","t"])
    J=get_SSPZ_from_polynomials([x/16 + 5*x^2/64 + 217/704 + 289/704*s, y/16 + x/8 + 10/16*t])

    return fPstart,FPZ,J
end
fPstart,FPZ,J=exemplejoinzono()

R=RealField()
get_polynomials_from_SSPZ(J,R)
plot([J,fPstart,FPZ],nsdiv=20)

function testjoinzonopolynomial()
    R=RealField()
    S,(x,y)=polynomial_ring(R,["x","y"])

    p1=1/4*(x+x^2)
    p2=1/4*(y+x)
    Pstart=get_SSPZ_from_polynomials([x; y])
    fPstart=poly_apply_on_SSPZ(Pstart,[p1,p2],R)
    FPZ=poly_apply_on_SSPZ(fPstart,[p1,p2],R)
    J=zonotopic_join(FPZ,fPstart,"BranchAndBoundEnclosure")
    println(get_polynomials_from_SSPZ(J,R))
    return fPstart,FPZ,J
end

fPstart,FPZ,J=testjoinzonopolynomial()
plot([J,fPstart,FPZ],nsdiv=20)