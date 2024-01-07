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
#using TaylorModels
using SumOfSquares
#using AffineArithmetic
#using ProfileVega

Profile.init(delay=0.01)
include("conversions.jl")
include("joins.jl")
include("plotsample.jl")
include("polynomap.jl")
include("matlabmatrix.jl")

const to = TimerOutput();
list=Int[]



function main()
    R=RealField()
    S,(x,y)=PolynomialRing(R,["x","y"])

    p1=x^3 -0.5*x^2+0.5
    p2=y^3 -0.5*y^2+0.5

    p6=(3/5*x +4/5*y)^3 -0.5*(3/5*x +4/5*y)^2+0.5
    p7=((-4/5)*x+3/5*y)^3 -0.5*((-4/5)*x+3/5*y)^2+0.5

    chatal1= x+y
    chatal2=-0.5952 + x^2

    p8=1/4*x^2+1/2 -y
    p9=2*y-y^2 + x

    henon1=1 - 1.4*x^2 + y
    henon2=0.3*x

    lineaire1=-0.32*x+0.32*y
    lineaire2= -0.42*x -0.92*y

    step=0.25
    vanpol1=step*y + x
    vanpol2=step*(1-x^2)*y - step*x + y

    brusselator1=step*(x)
    brusselator2=step*(y)

    lokta1=step*(x*(1.5 - y)) + x
    lokta2=step*(-y*(3 - x)) + y

    parillo1= step*(-x+x*y) + x
    parillo2=step*(-y) + y

    PHenon=get_SSPZ_from_polynomials([1/5*x ,1/5*y])
    Pparillo=get_SSPZ_from_polynomials([1/10*x+12 ,1/10*y+2])
    PChatal=get_SSPZ_from_polynomials([1/6*x+1/6 ,1/6*y+1/6])
    Plineaire=get_SSPZ_from_polynomials([x-2,y+1])
    P1=get_SSPZ_from_polynomials([x + 3  ,y+4])
    PVanPol=get_SSPZ_from_polynomials([0.15*x + 1.4  ,0.05*y+2.30])
    Plokta=get_SSPZ_from_polynomials([1/5*x + 5  ,1/5*y+2])


    start_time = now()
    fin=iterate_polynomials_over_PZ([chatal1,chatal2],PChatal,6,1,R,"bary",max_order=5000000,power=1)

    end_time = now()
    elapsed = end_time - start_time
    println("temps des iterations:", elapsed)
    fin=reverse(fin)
    fini=fin[end]
    println("nombre de variables à la fin: ",size((fini).E)[1])
    println("nombre de monomes à la fin: ",size(fini.E)[2])
    somme=sum(fini.E,dims=1)
    som=vec(somme)
    println("degré maximal à la fin: ",maximum(som))
    #affiche_liste(get_polynomials_from_SSPZ(fini,R))

    #plot_multiple(fin,R,"Documents/julia/plots_julia/lineaire^1_4iter_joinzonoBRABOUN1_100000pts_",nbpoints=100000)
    #plot_sampling(fini,R,"Documents/julia/plots_julia/lineaire^1_4iter_joinzonoBRABOUN1_100000pts_",nbpoints=100000)
    
    #derniere=poly_apply_on_SSPZ(fini,[lineaire1,lineaire2],R)

    #plot_multiple([derniere,fini],R,"Documents/julia/plots_julia/lineaire^1_5iter_joinbary0_Inclusion?_1000000pt",nbpoints=1000000)
    end_time2 = now()
    elapsed = end_time2 - end_time
    println("temps de plot:", elapsed)

    #MatlabMatrix(fini,"Documents/julia/Traduct_Matlab/Henon_6iter_zono1_doublereduc50.txt")
    

    return fini
    #return derniere
end
r=main()
r.E



ProfileView.@profview main()
#Profile.clear()
#@profile main()
#Profile.print()
#pprof()
@allocations main()

R=RealField()
S,(x,y)=PolynomialRing(R,["x","y"])
using BernsteinExpansions
using IntervalArithmetic
dom = IntervalBox(-1..1, 2)

chatal1= x+y
enclose(chatal1,dom)
range(chatal1,dom)
chatal2=-0.5952 + x^2
point=[-2.0,2.0]
point=[evaluate(chatal1,point),evaluate(chatal2,point)]
evaluate(
    chatal1,point
)
r=main()
PZ_reduc=SimpletoSPZ(r)
RED=reduce_order(PZ_reduc,50)
rb=SPZ_to_SimpleSPZ(RED)
plot_sampling(rb,R,"Documents/julia/plots_julia/Henonreduit50_6iter_joinbary1_100000pts_",nbpoints=100000)
plot_sampling(r,R,"Documents/julia/plots_julia/Henon_6iter_joinbary1_100000pts_",nbpoints=100000)
r.E
r.G
G
MatlabMatrix(r,"Documents/julia/Traduct_Matlab/Chatala^2_3iter_joinbary0_noreduc.txt")




R=RealField()
S,(x,y,s,t)=PolynomialRing(R,["x","y","s","t"])

Trian2=get_SSPZ_from_polynomials([x*y,x])

Trian3=get_SSPZ_from_polynomials([-x*y+3,-x+3])
Trian4=get_SSPZ_from_polynomials([-s*t+3,-s+3])
Trian5=barycentric_join(Trian2,Trian4,3)
Trian6=barycentric_join(Trian2,Trian3,3)
plot_multiple([Trian2,Trian3],R,"Documents/julia/plots_julia/exemplepourjoin2",nbpoints=100000)
plot_multiple([Trian5,Trian6],R,"Documents/julia/plots_julia/exemplepourjoinVSconvexhull2",nbpoints=3000000)

Trian3.G

@view Trian3.G[1,:]==@view Trian3.G[2,:]
A,(x,y)=PolynomialRing(R,["x","y"])
Test1=get_SSPZ_from_polynomials([x+y,-0.5+x^2])
Test2=get_SSPZ_from_polynomials([x+y+x^2-0.5,(x+y)^2])
get_polynomials_from_SSPZ(union_pol_zono(Test1,Test2,R),R)
get_polynomials_from_SSPZ(zonotopic_join(Test1,Test2),R)
REST=zonotopic_join(Test1,Test2)
get_polynomials_from_SSPZ(REST,R)
REST.G
n1,n2=size(Test1.G)
e1,e2=size(Test1.E)
a1(x)=sum(Test1.G[2,j]*prod(x[k]^Test1.E[k,j] for k in 1:e1) for j in 1:n2)+Test1.c[2]
a1([1,1])