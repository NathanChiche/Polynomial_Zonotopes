using LazySets
using Nemo 
using IntervalArithmetic
using IntervalContractors

include("conversions.jl")
include("joins.jl")
include("plotsample.jl")
include("polynomap.jl")
include("matlabmatrix.jl")
include("reduction.jl")
include("bernstein.jl")


"""TESTS SUR LES REMOVE REDUNDANT GENERATORS ET LA REDUCTION"""
SP1=SimpleSparsePolynomialZonotope([0 , 0.0],[ 0 3 2.5 2.0 2 2 2 2 1 0 0.2 3; 2 2 1.0 0 1.3 1 1 3 1 1 1 3],[2 2 1 1 1 1 2 1 2 0 1 3; 5 2 1 0 3 3 3 4 2 1 1 3])
SP3=deepcopy(SP1)
SP3.E
SP3==SP1
SP2=remove_useless_terms!(SP1)
plot([SP3,SP1],nsdiv=28)
plot(SP1,nsdiv=28)
plot([SP2,SP3],nsdiv=28)
SP2==SP3

SP2.E
SP3.G==SP1.G

size(SP1.G)
RSP1=Simple_reduce_order(SP2,5)
RSP2=Simple_reduce_order(SP2,4,nor=1)
RSP2.G
RSP1.G
RSP1==RSP2
plot([RSP2,SP2],nsdiv=25)
plot([RSP1,SP2],nsdiv=25)
SP2.G
SP2.E
RSP1.G
RSP1.E

plot_multiple([RSP2,SP2],R,"Documents/julia/plots_julia/comparaisonreductionnor=1_ord=5")
plot_multiple([RSP1,SP2],R,"Documents/julia/plots_julia/comparaisonreductionnor=2_ord=3")


R=RealField()
S,(x,y)=PolynomialRing(R,["x","y"])

"""CET EXEMPLE DE JOIN ESR TRES JOLI"""

P1=get_SSPZ_from_polynomials([x^2+2*y^2+1 ; x*y])
P2=get_SSPZ_from_polynomials([x^3+y^2;2+x+x*y])
P3=barycentric_join(P1,P2)
get_polynomials_from_SSPZ(P1,R)
get_polynomials_from_SSPZ(P2,R)
@show(get_polynomials_from_SSPZ(P3,R))

P4=zonotopic_join(P1,P2,"bernstein")
@show(get_polynomials_from_SSPZ(P4,R))
plot_multiple([P4,P1,P2],R,"Documents/julia/plots_julia/deuxiemetestjoinzono")
plot(P4,nsdiv=15)

P5=get_SSPZ_from_polynomials([0.4*x^2+2*y^2+1 ; x*y])
P6=get_SSPZ_from_polynomials([0.7*x^3+1.3*y^2; 2+0.8*x+x*y])
P7=zonotopic_join(P5,P6,"bernstein")

plot_multiple([P7,P5,P6],R,"Documents/julia/plots_julia/quatrtestjoinzono")
plot([P5,P6])
@show(P6.G,P6.E)
plot(P6)

PZ=SimpleSparsePolynomialZonotope([0.0;0],[0.7 0.0 0.0 0.0; 0.0 0.0 0.0 0.8],[3 0 1 1; 0 2 1 0; 0 0 0 0; 0 0 0 0])
get_polynomials_from_SSPZ(PZ,R)
plot(PZ,nsdiv=16)

function affiche_polynomes(pol)
    #exponents=collect(exponent_vector(pol,j) for j in 1:length(pol))
    str=""
    println(str)
    for i in 1:length(pol)
        c=coeff(pol,i)
        m=Nemo.monomial(pol,i)
        if sign(Float64(c))>= 0
            #println(string(Float64(c)),Float64(c))
            str=str*" + "*string(Float64(c))*string(m)
        else
            str=str*string(Float64(c))*string(m)
        end
    end
    println(str)
end

function affiche_liste(list_poly)
    for p in list_poly
        affiche_polynomes(p)
    end
end



plot_multiple([P3,P2,P1],R,"Documents/julia/plots_julia/testjoinbisbary")
X=get_SSPZ_from_polynomials([-x^2+y;x*y])
Y=get_SSPZ_from_polynomials([-1/5*x^2+y;x*y])
Z=zonotopic_join(X,Y,"bernstein")

plot_multiple([Z,Y,X],R,"Documents/julia/plots_julia/premiertest")
affiche_liste(get_polynomials_from_SSPZ(Z,R))

S,(x,y,s,t)=PolynomialRing(R,["x","y","s","t"])
X=get_SSPZ_from_polynomials([x*y+y^2,y^2])
Y=get_SSPZ_from_polynomials([x*y+2*y^2+2,y^2+2])
Z=get_SSPZ_from_polynomials([2+x*y+y^2+s,y^2+1.5+0.5*t])


plot_multiple([Z,X,Y],R,"Documents/julia/plots_julia/rangetest")

"""TEST POUR INCLUSION ET METHODES INTERVALLES"""

Inc=get_SSPZ_from_polynomials([x^2+x*y,x*y])
plot_sampling(Inc,R,"Documents/julia/plots_julia/imagepourinclusioncontracteurs")

Box::IntervalBox
typeof(Box)
sawtooth(Box)
IntervalContractors.constant_contractor 

"""TEST POUR JOLI JOIN ZONO"""
AZ=get_SSPZ_from_polynomials([ x + x*y + x^2 - 1; 2*y + x^2 + 0.5*x*y + 1 ])
QZ=get_SSPZ_from_polynomials([0.3*x^4 + x^2 + y^3  + x*y + x + 1; 0.5*y^4 + 0.5*x^2 + 2*y + x^2*y + 1.5*x*y + y ])
SZ=get_SSPZ_from_polynomials([(x^2+y^3+x*y + x+1)^2*0.5*y^4 + 2*y + x^2*y+ 1.5*x^3*y + (0.5*y^4 + 2*y + x^2*y+ 1.5x^3*y)*(x^2+y^3+x*y + x+1) + 0.3*(x^2+y^3+x*y + x+1) ;
 0.6*(0.5*y^4 + 2*y + x^2*y+ 1.5x^3*y ) + 1.5*(0.5*y^4 + 2*y + x^2*y+ 1.5x^3*y)* (x^2+y^3+x*y + x+1) + x^2+y^3+x*y + x+1 + 2])
SZ=get_SSPZ_from_polynomials([(x^2+y^3+x*y + x+1)^2*0.5*y^4 + 2*y + x^2*y+ 1.5*x^3*y + (0.5*y^4 + 2*y + x^2*y+ 1.5x^3*y)*(x^2+y^3+x*y + x+1) + 0.3*(x^2+y^3+x*y + x+1) ;
0.6*(0.5*y^4 + 2*y + x^2*y+ 1.5x^3*y ) + 1.5*(0.5*y^4 + 2*y + x^2*y+ 1.5x^3*y)* (x^2+y^3+x*y + x+1) + x^2+y^3+x*y + x+1 + 2])
SZ.E
QZ.E
ZZ=zonotopic_join(AZ,QZ,"bernstein")
BZ=barycentric_join(AZ,QZ)
ZZ.E

po=x^2 + x*y + x + x*y
qo=y^2 + x + y

T1=get_SSPZ_from_polynomials([ x^2 + x*y + x + x*y; y^2 + x + y ])
T2=poly_apply_on_SSPZ(T1,[po,qo],R)
T2.G
T2.E
T1.E
@show(get_polynomials_from_SSPZ(T2,R))

T3=zonotopic_join(T1,T2,"bernstein")
T3.E
T3.G
plot_multiple([T3,T1,T2],R,"Documents/julia/plots_julia/joinzonocool7",nbpoints=50000)
plot_multiple([ZZ,AZ,QZ],R,"Documents/julia/plots_julia/joinzonocool6",nbpoints=300000)
plot_multiple([BZ,QZ,AZ],R,"Documents/julia/plots_julia/joinbarycool3")

poq= 0.2*x*y + x +1
qoq=0.15*x*y + x + y +2
Box=get_SSPZ_from_polynomials([x , y])
R1=poly_apply_on_SSPZ(Box,[poq,qoq],R)
R2=poly_apply_on_SSPZ(R1,[poq,qoq],R)
@show(get_polynomials_from_SSPZ(R2,R))
R3=zonotopic_join(R1,R2,"bernstein")
R2.G
R2.E
plot_multiple([R3,R2,R1],R,"Documents/julia/plots_julia/joinzonocool9",nbpoints=500000)