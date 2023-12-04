using LazySets
using Nemo
using Plots
include("conversions.jl")

function test()
    R=RealField()
    S,(x,y)=PolynomialRing(R,["x","y"])

    p1=y^3 -0.5*x^2+0.5
    p2=y^3 -0.5*y*x^2+0.5
    P=get_SSPZ_from_polynomials([p1,p2])
    P2=get_polynomials_from_SSPZ(P,R)
    return P
end

z=test()
plot(z)