function evaluate_polynomials_on_vector(polynomials,vector::Vector{Float64})
    return [Float64(evaluate(polynomials[i],vector)) for i in 1:length(polynomials)]
end

function plot_sampling(PZ::SimpleSparsePolynomialZonotope,field::Field,filename::String;nbpoints=300000,xlim=nothing,ylim=nothing)
    """enregistre dans filename le tracé de PZ avec nbpoints différents"""
    nb_vars=size(expmat(PZ))[1]
    nb=size(genmat(PZ))[1]
    polynomes=get_polynomials_from_SSPZ(PZ,field)
    step=1/nbpoints
    liste=collect(rand(-1.0:step:1.0,nb_vars) for i in 1:nbpoints)
    points = [evaluate_polynomials_on_vector(polynomes,v) for v in liste]
    #pl=plot(first.(points),last.(points),legend=false)
    x=collect(points[i][1] for i in 1:length(points))
    y=collect(points[i][2] for i in 1:length(points))
    pl=scatter(x,y,#=xlimits=xlim,ylimits=ylim,=#legend=false)
    savefig(pl,filename)
    return pl
end    

function plot_multiple(liste_PZ,field::Field,filename::String;nbpoints=300000,xlim=nothing,ylim=nothing)
    i=1
    for PZ in liste_PZ
        println("coucou")
        nb_vars=size(expmat(PZ))[1]
        nb=size(genmat(PZ))[1]
        polynomes=get_polynomials_from_SSPZ(PZ,field)
        step=1/nbpoints
        liste=collect(rand(-1.0:step:1.0,nb_vars) for i in 1:nbpoints)
        points = [evaluate_polynomials_on_vector(polynomes,v) for v in liste]
        #pl=plot(first.(points),last.(points),legend=false)
        x=collect(points[i][1] for i in 1:length(points))
        y=collect(points[i][2] for i in 1:length(points))
        if i==1
            if xlim!==nothing
                scatter(x,y,xlimits=xlim,ylimits=ylim,legend=false)
            else 
                scatter(x,y,legend=false)
            end
        else
            if xlim!==nothing
                scatter!(x,y,xlimits=xlim,ylimits=ylim,legend=false)
            else 
                scatter!(x,y,legend=false)
            end
        end
        i=i+1
    end
    savefig(filename)
end


#=pas=1/10
pointspas=collect(rand(-1.0:pas:1.0,5) for i in 1:100)
@show(pointspas)

pointspas2=collect(rand(-1.0:1.0,3) for i in 1:10)
@show(pointspas2)

@time pointspas=collect(rand(-1.0:1/10000:1.0,5) for i in 1:100000)
@time samples = [rand(3) .* 2 .- 1 for _ in 1:100000]=#