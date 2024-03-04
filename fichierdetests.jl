using LazySets
using Nemo 

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
S,(x,y,s,t)=PolynomialRing(R,["x","y","s","t"])

"""CET EXEMPLE DE JOIN ESR TRES JOLI"""

P1=get_SSPZ_from_polynomials([x^2+2*y^2+1;x*y])
P2=get_SSPZ_from_polynomials([x^3+y^2;2+x+x*y])
P3=barycentric_join(P1,P2)
get_polynomials_from_SSPZ(P1,R)
get_polynomials_from_SSPZ(P2,R)
@show(get_polynomials_from_SSPZ(P3,R))
affiche_liste(get_polynomials_from_SSPZ(P1,R))
affiche_liste(get_polynomials_from_SSPZ(P2,R))
affiche_liste(get_polynomials_from_SSPZ(P3,R))

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


plot_multiple([P3,P2,P1],R,"Documents/julia/plots_julia/testjoinbis")

X=get_SSPZ_from_polynomials([x*y+y^2,y^2])
Y=get_SSPZ_from_polynomials([x*y+2*y^2+2,y^2+2])
Z=get_SSPZ_from_polynomials([2+x*y+y^2+s,y^2+1.5+0.5*t])


plot_multiple([Z,X,Y],R,"Documents/julia/plots_julia/rangetest")