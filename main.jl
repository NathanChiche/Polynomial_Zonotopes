using Nemo
using LazySets
using Plots
using Random
using Dates

include("conversions.jl")
include("joins.jl")
include("plotsample.jl")
include("polynomap.jl")

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
    fin=iterate_polynomials_over_PZ([chatal1,chatal2],PChatal,4,2,R,"zono",max_order=25000000)
    end_time = now()
    elapsed = end_time - start_time
    println("temps des iterations:", elapsed)
    #fin=reverse(fin)
    fini=fin[end]
    println("nombre de variables à la fin: ",size((fini).E)[1])
    println("nombre de monomes à la fin: ",size(fini.E)[2])
    somme=sum(fini.E,dims=1)
    som=vec(somme)
    println("degré maximal à la fin: ",maximum(som))
    #affiche_liste(get_polynomials_from_SSPZ(fini,R))

    plot_multiple(fin,R,"Documents/julia/plots_julia/chatal^1_4iter_joinzonoraffine2_30000pts_max_order=40000000bis",nbpoints=30000)
    #plot_sampling(fini,R,"Documents/julia/plots_julia/Parillo^1_4iter_joinbary_100000pts_",nbpoints=100000)
    #plot_multiple(fin,R,"Documents/julia/plots_julia/lineaire^1_4iter_joinbary1_40000pts_max_order=25_x-2_y+1",nbpoints=40000)
    
    #derniere=poly_apply_on_SSPZ(fini,[lineaire1,lineaire2],R)

    #plot_multiple([derniere,fini],R,"Documents/julia/plots_julia/lineaire^1_5iter_joinbary0_Inclusion?_1000000pt",nbpoints=1000000)
    end_time2 = now()
    elapsed = end_time2 - end_time
    println("temps de plot:", elapsed)

    
    

    return fini
    #return derniere

end
r=main()