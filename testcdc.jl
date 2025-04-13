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
include("conversions.jl")
include("joins.jl")
include("plotsample.jl")
include("polynomap.jl")
include("matlabmatrix.jl")
include("reduction.jl")
#include("bernstein.jl")
include("testinclusion_robrange.jl")

#@polyvar x1,x2

function rangesfrompolynomialzonotope(PZ::SimpleSparsePolynomialZonotope,tolerance::Float64,maxde::Int64)
    d=length(PZ.c)
    nbvars=size(PZ.E)[1]
    ngener=size(PZ.G)[2]
    domain=IntervalBox(-1..1, nbvars)
    intervallesurapprox=Array{IntervalArithmetic.Interval{Float64}}(undef,d)
    for i in 1:d
        a1(v)=sum(PZ.G[i,j]*prod(v[k]^PZ.E[k,j] for k in 1:nbvars) for j in 1:ngener)+PZ.c[i]
        intervallesurapprox[i]=enclose(a1,domain,BranchAndBoundEnclosure(tol=tolerance,maxdepth=maxde))
        #println(enclose(a1,domain,MooreSkelboeEnclosure()))
    end
    return intervallesurapprox
end

rangesfrompolynomialzonotope(get_SSPZ_from_polynomials([x,y]),1e-3,10)

R=RealField()
S,(x,y)=polynomial_ring(R,["x","y"])

function testinclusion(p1,p2,dx,dy,a1,a2,nbiter,withoutjoin,join,maxo,inclusiontes,solve,powerf,linearsys)
    Pstart=get_SSPZ_from_polynomials([a1*x + dx; a2*y + dy])
    fin=iterate_polynomials_over_PZ([p1,p2],Pstart,nbiter,withoutjoin,R,join,max_order=maxo,inclusiontest=inclusiontes,solver=solve,power=powerf,linearsystem=linearsys)
    fini=fin[end]
    next=poly_apply_on_SSPZ(fini,[p1,p2],R)
    return fin,next
end

function testlineaire(nbiter,withoutjoin,join,inclusiotes)
    p1=0.25*x+0.2*y +7/20
    p2=0.2*x+0.2*y 
    dx=0
    dy=0
    zono,nextzono=testinclusion(p1,p2,dx,dy,1,1,nbiter,withoutjoin,"zono",200000,1,"NaturalEnclosure",1,true)
    bary,nextbary=testinclusion(p1,p2,dx,dy,1,1,nbiter,withoutjoin,"bary",200000,0,"NaturalEnclosure",1,true)
    println()
    return zono,bary

end