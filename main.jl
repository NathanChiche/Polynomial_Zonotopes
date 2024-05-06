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



Profile.init(delay=0.01)
include("conversions.jl")
include("joins.jl")
include("plotsample.jl")
include("polynomap.jl")
include("matlabmatrix.jl")
include("reduction.jl")
include("bernstein.jl")

const to = TimerOutput();

function main()
    R=RealField()
    S,(x,y)=PolynomialRing(R,["x","y"])

    p1=x^3 -0.5*x^2+0.5
    p2=y^3 -0.5*y^2+0.5

    p6=(3/5*x +4/5*y)^3 -0.5*(3/5*x +4/5*y)^2+0.5
    p7=((-4/5)*x+3/5*y)^3 -0.5*((-4/5)*x+3/5*y)^2+0.5

    p8=1/5*(3/5*x +4/5*y)^2
    p9=1/5*((-4/5)*x+3/5*y)^2

    g1=3/5*(0.4*x^2) - 4/5*(0.6*y^2)
    g2= -4/5*(0.4*x^2) + 3/5*(0.6*y^2)

    chatal1= x+y
    chatal2=-0.5952 + x^2
    chatalb1= 1/(1.2)*(x+y)
    chatalb2=1/(1.2)*(-0.5952 + x^2)


    #chatalb1= 3/5*x+ -4/5*y+y^2
    #chatalb2= 3/5*x^2+ 4/5*y

    #p8=1/4*x^2+1/2 -y
    #p9=2*y-y^2 + x

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

    basique1=1/7*(x^2+x*y+x)
    basique2=1/6*(y+x*y)

    basique1bis=1/10*((3/5*x +4/5*y)^2+(3/5*x +4/5*y)*((-4/5)*x+3/5*y))
    basique2bis=1/12*(((-4/5)*x+3/5*y)+(3/5*x +4/5*y)*((-4/5)*x+3/5*y))

    mix_linch1=lineaire1 + 0.2*chatal1
    mix_linch2=lineaire2 + 0.2*chatal2


    

    PHenon=get_SSPZ_from_polynomials([1/5*x ,1/5*y])
    Pparillo=get_SSPZ_from_polynomials([1/10*x+12 ,1/10*y+2])
    PChatal=get_SSPZ_from_polynomials([1/6*x+1/6 ,1/6*y+1/6])
    PChatalb=get_SSPZ_from_polynomials([1/6*x+1/3 ,1/6*y+1/3])
    Plineaire=get_SSPZ_from_polynomials([x-2,y+1])
    P1=get_SSPZ_from_polynomials([1/20*x + 0.4 ,1/20*y+0.4])
    P2=get_SSPZ_from_polynomials([1/5*x + 0.5 ,1/5*y+0.5])
    P89=get_SSPZ_from_polynomials([1/5*x + 1.5 ,1/5*y+1])
    PVanPol=get_SSPZ_from_polynomials([0.15*x + 1.4  ,0.05*y+2.30])
    Plokta=get_SSPZ_from_polynomials([1/5*x + 5  ,1/5*y+2])
    Pbasique=get_SSPZ_from_polynomials([1/3*x + 1  ,1/3*y+1])


    start_time = now()
    fin=iterate_polynomials_over_PZ([lineaire1,lineaire2],Plineaire,8,1,R,"bary",max_order=20,power=1,solver="bernstein")
   
    end_time = now()
    elapsed = end_time - start_time
    println("temps des iterations:", elapsed)
    #fin=reverse(fin)
    fini=fin[end]
    fin=reverse(fin)
    println("nombre de variables à la fin: ",size((fini).E)[1])
    println("nombre de monomes à la fin: ",size(fini.E)[2])
    somme=sum(fini.E,dims=1)
    som=vec(somme)
    println("degré maximal à la fin: ",maximum(som))#affiche_liste(get_polynomials_from_SSPZ(fini,R))

    #plot_multiple(fin,R,"Documents/julia/plots_julia/mixlincha_3it_nojoin_ordre=500",nbpoints=40000)
    #plot_sampling(fini,R,"Documents/julia/plots_julia/chatal^1_4iter_joinbarymatriciel1_maxorder=50_100000pts_",nbpoints=100000)
    

    #plot_multiple([derniere,fini],R,"Documents/julia/plots_julia/lineaire^1_5iter_joinbary0_Inclusion?_1000000pt",nbpoints=1000000)
    

    MatlabMatrix(fini,"Documents/julia/Traduct_Matlab/fini_lineaire_8iter_bary1_reduc20.txt","flin8_1") 
    next=poly_apply_on_SSPZ(fini,[lineaire1,lineaire2],R)
    #next=poly_apply_on_SSPZ(next,[basique1,basique2],R)
    MatlabMatrix(next,"Documents/julia/Traduct_Matlab/next_lineaire_8iter_bary1_reduc20.txt","nlin8_1") 

    

    return fini,next
    #return derniere
    return fini
end

@time r=main()



@polyvar x[1:4]
x[1]
r[1].G
R=RealField()
get_polynomials_from_SSPZ(r[1],R)
plot_multiple(r,R,"Documents/julia/plots_julia/basiquebisjoinzono5iter_2_reduc50",nbpoints=40000)
chatalb1= 3/5*x+ -4/5*y+y^2
    chatalb2= 3/5*x^2+ 4/5*y
rb=poly_apply_on_SSPZ(r,[chatalb1,chatalb2],R)


rbj=zonotopic_join(r,rb,"bernstein")
r.c
rb2=zonotopic_join(r,rb,"BranchAndBoundEnclosure")
rb3=zonotopic_join(r,rb,"NaturalEnclosure")

plot_multiple([rbj,r,rb],R,"Documents/julia/plots_julia/testsjoinpluscentre")
plot_multiple([rb2,r,rb],R,"Documents/julia/plots_julia/testsjoinbranch")

box=SimpleSparsePolynomialZonotope([0, 0.],[1 0;0 1.],[1 0;0 1])
box2=SimpleSparsePolynomialZonotope([5, 3.],[2 0;0 1.],[1 0;0 1])
T
box4=bernstein_zonotopic_join(box,box2,"bernstein",R)
box5=union_pol_zono(box,box2,R)
get_polynomials_from_SSPZ(box5,R)
box4
typeof(box4[1])
get_SSPZ_from_polynomials(box4)

box3=zonotopic_join(box,box2,"bernstein")
box4=zonotopic_join(box,box2,"BranchAndBoundEnclosure")
box3.G
plot_multiple([box3,box2,box],R,"Documents/julia/plots_julia/testjoincasbasedecal")




S1=SimpleSparsePolynomialZonotope([0 , 0.],[1 1 0 0 0; 0 0 1 1 1.],[2 0 3 0 1;1 1 2 2 0])
plot_sampling(S1,R,"Documents/julia/plots_julia/S1")
S2=SimpleSparsePolynomialZonotope([0 , 0.],[1 0.5 0 0 ; 0 0 1 2],[1 0 3 1; 0 1 2 0 ])
S2=remove_useless_terms!(S2)
S2=SimpleSparsePolynomialZonotope([0 , 0.],[1 0.5 0 ; 2 0 1],[1 0 3; 0 1 2 ])
plot_sampling(S2,R,"Documents/julia/plots_julia/S2")



get_polynomials_from_SSPZ(S2,R)

S3=zonotopic_join(S1,S2,"bernstein")
get_polynomials_from_SSPZ(S3,R)
plot_multiple([S4,S1,S2],R,"Documents/julia/plots_julia/S4S2S1")
S4=barycentric_join(S1,S2)
plot_sampling(S4,R,"Documents/julia/plots_julia/S4")


rred=Simple_reduce_order(r,50)
rred.E
plot_multiple([rred,r],R,"Documents/julia/plots_julia/test_reductionoversimple_")

r.E
rb=remove_unused_variables(r)
rb.E


ProfileView.@profview main()
#Profile.clear()
#@profile main()
#Profile.print()
#pprof()
@allocations main()


R=RealField()
S,(x,y)=PolynomialRing(R,["x","y"])
g=gens(S)
ze=S(0)
T,(x,y,z,s)=PolynomialRing(R,["x","y","z","s"])
zer=T(ze)
g[1]*g[2]
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


P1=SparsePolynomialZonotope([0 , 0.0],[3 2.5 2.0 2 2 2; 2 1.0 0 1.3 1 1],[0.0 0.0 ; 0.0 0],[2 1 1 1 1 2; 4 1 0 3 3 3])
SP1=SPZ_to_SimpleSPZ(P1)
x
SP1.E
list=[Monomial(x, SP1.E[:,i]) for i in 1:2]

SP1=SimpleSparsePolynomialZonotope([0 , 0.0],[3 2.5 2.0 2 2 2 2 1 0 0.2; 2 1.0 0 1.3 1 1 3 1 1 1],[2 1 1 1 1 2 1 2 0 1; 4 1 0 3 3 3 4 2 1 7])
size(SP1.G)
SP1=remove_useless_terms!(SP1)
size(SP1.G)
RSP1=Simple_reduce_order(SP1,3)
RSP2=Simple_reduce_order(SP1,3,nor=1)
RSP1==RSP2
@show(RSP1.G)
@show(RSP2.G)
@show(RSP1.E)
@show(RSP2.E)
plot([RSP2,RSP1],nsdiv=20)
plot_multiple([RSP1,RSP2,SP1],R,"Documents/julia/plots_julia/deuxreductions",nbpoints=100000)
G=SP1.G
n=norm
norms = [n(g) for g in eachcol(G)]
norms2=[norm_combination(g) for g in eachcol(G)]

@show(SP1.G)
@show(RSP1.G)
@show(SP1.E)
@show(RSP1.E)
get_polynomials_from_SSPZ(SP1,R)
get_polynomials_from_SSPZ(RSP1,R)



plot
RP1=reduce_order(P1,2)
plot(RP1)
plot(P1,nsdiv=20,parti)
plot(SP1)

include("deuxieme.jl")

LazySets.center(SP1)
ZO=Zonotope([0 , 0.0],[1 2.0 2 2 2; 1.0 0 1 1 1])
LazySets.center(ZO)
ZO.generators
ZO.center
reduce_order(P1,2)
typeof(LazySets.GIR05())
a=LazySets.GIR05()
subtypes(AbstractReductionMethod)
function looping(Z,method=LazySets.GIR05())
    for i in 1:10000
        reduce_order(Z,1,method)
    end
    return Z
end

ProfileView.@profview looping(ZO)

R=RealField()
S,(x,y,s,t)=PolynomialRing(R,["x","y","s","t"])

Trian2=get_SSPZ_from_polynomials([x*y,x])
overapproximate(Trian2, Zonotope)
norm(Trian2.G[:,1])

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
