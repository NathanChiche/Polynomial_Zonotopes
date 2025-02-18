using Symbolics
using Nemo
using IntervalArithmetic
using StaticArrays
using Combinatorics
using ForwardDiff

include("polzono_to_symbo.jl")


function O(range_Jf, i)
    return @interval(-(abs(range_Jf[i]).hi),abs(range_Jf[i]).hi)
end

function I(range_Jf,i)
    return @interval(-(abs(range_Jf[i]).lo),abs(range_Jf[i]).lo)
end

function outer_approx(range_Dg,center,quantifiers)
    min=center
    max=center
    p=length(quantifiers)
    for i in 0:p-1
        if quantifiers[p-i]=="exists"

            temp=O(range_Dg,p-i)
            @show(temp)
            min=min+temp.lo
            max=max+temp.hi
        else
            temp=I(range_Dg,p-i)
            min=min+temp.hi
            max=max+temp.lo
            if max<min
                #println("resultat de sur-approx vide")
                return "false"
            end
        end  
    end
    #println(min,typeof(min))
    #println(max,typeof(max))
    return interval(min,max)
end

function inner_approx(range_Dg,center,quantifiers)
    min=center
    max=center
    p=length(quantifiers)
    for i in 0:p-1
        if quantifiers[p-i]=="exists"
            temp=I(range_Dg,p-i)
            #@show(p-i,temp)
            min=min+temp.lo
            max=max+temp.hi
        else
            temp=O(range_Dg,p-i)
            #@show(p-i,temp)
            min=min+temp.hi
            max=max+temp.lo
            if max<min
                #println("resultat de sous-approx vide")
                return "false"
            end
        end  
    end
    return interval(min,max)
end
#on doit faire un tableau qui liste dans quel ordre on intervient

function centers(Polynomes1,Polynomes2,n1,n2)
    @assert length(Polynomes1)==length(Polynomes2)
    z1=zeros(n1)
    z2=zeros(n2)
    centers=[]
    for i in 1:length(Polynomes1)
        push!(centers,Float64(evaluate(Polynomes1[i],z1)-evaluate(Polynomes2[i],z2)))
    end
    return centers
end

function call_multiple_outers_bis(rangelist,centers,quantifiers)
    @show(rangelist)
    #println("centers :",centers)
    intervals=[]
    #m=length(quantifiers)
    for j in 1:length(rangelist)
        i=outer_approx(rangelist[j],centers[j],quantifiers)
        if !(0 in i)
            #println("0 n'est même pas dans l'outer approx")
            return "false"
        end
        push!(intervals,i)
        #push!(ranges_Df,range_Dg)
    end
    return intervals
end



function call_multiple_inners_bis(rangelist,centers,quantifiers)
    #println("on entre dans la sous approx multidimensionnelle")
    #println("centers :",centers)
    intervals=[]
    #m=length(x)
    for j in 1:length(rangelist)
        i=inner_approx(rangelist[j],centers[j],quantifiers)
        if !(0 in i) || i=="false"
            #println("0 n'est pas dans l'inner pour quant: ",quantifiers)
            return "false"
        end
        push!(intervals,i)
        #push!(ranges_Df,range_Dg)
    end
    return intervals
end

function widthofintervalbox(intervalbox,nbvars)
    width=maximum([intervalbox[i].hi-intervalbox[i].lo for i in 1:nbvars])
    return width
end


function compute_derivatives_approx(fun,dim,intervalbox)
    ranges=[]
    for a in 1:dim
        push!(ranges,ForwardDiff.gradient(fun[a],intervalbox))
    end
    return ranges
end

function compute_ranges(listgradient,dim,intervalbox)
    ranges=[]
    for a in 1:dim
        @show(listgradient[a](intervalbox))
        push!(ranges,listgradient[a](intervalbox))
    end
    println("range apres calcul sur boite 2d")
    @show(ranges)
    return ranges
end


function zeros_in_intervalvectors(intervals,dim)
    for o in 1:dim
        if !(0 in intervals[o])
            return 0
        end
    end
    return 1
end

function intervals_centers(intervals)
    m=zeros(Float64,length(intervals))
    for i in 1:length(intervals)
        m[i]=(intervals[i].hi+intervals[i].lo)/2
    end
    return m
end


function bisect_at_component(interva, component::Int)
    # Find the midpoint of the chosen dimension
    a, b = interva[component].lo, interva[component].hi
    m = (a + b) / 2

    # Create two new intervals: one for each half of the bisection
    left = copy(interva)
    right = copy(interva)

    # Split the component at the midpoint
    left[component] = a..m
    right[component] = m..b

    return left, right
end

function paverobust(f,listgradie,dim,nbvars,epsilon,quantifiers)
    intervalbox=[interval(-1..1) for c in 1:nbvars]
    listintervalbox=[intervalbox]
    width=2.0
    tour=1
    #@show(dim,nbvars)
    while width>epsilon
        l=length(listintervalbox)
        for i in 1:l
            interval_to_check=popfirst!(listintervalbox)
            #@show(interval_to_check)
            #ranges=compute_ranges(listgradie,dim,interval_to_check)
            ranges=compute_derivatives_approx(f,dim,interval_to_check)
            #@show(ranges)
            center=[f[b](intervals_centers(interval_to_check)) for b in 1:dim]
            #@show(center)

            if call_multiple_outers_bis(ranges,center,quantifiers)==false
                println("pas dans la surapprox")
                return "false"
            else 
                vectinner=call_multiple_inners_bis(ranges,center,quantifiers)
                #@show(vectinner)
                if vectinner=="false"
                    left,right=bisect_at_component(interval_to_check,tour)
                    push!(listintervalbox,left,right)
                end
            end
            if length(listintervalbox)==0
                return 1
            end
        end
        tour=(tour)%(nbvars-dim)+1
        #@show(tour)
        if tour==1
            width=width/2
        end

        #@show(width)
    end
    return "false"
end 


function inclusion_test(FPZ,PZ,epsilon)
    dim=length(PZ.c)
    listfunc,listgrad,quantifiers,nbvars=functionalinclusion_polynomial_zonotopes_to_function(FPZ,PZ)
    inte=interval(-1..1)
    intervalbox=[inte for i in 1:nbvars]
    
    @show(listgrad[1](intervalbox))
    @show(listgrad[2](intervalbox))
    return paverobust(listfunc,listgrad,dim,nbvars,epsilon,quantifiers)

end

function test_primal() #ATTENTION: plutôt gros problème sur la tête des ranges jacobiennes qui doivent être des intervalles
    R=RealField()
    Anneau,(x,y)=PolynomialRing(R,["x","y"])
    pol1=[2*x^2+2*x]
    pol2=[3*x^2+5*x]
    inclusion_test(pol1,pol2)
end


@variables x, y, z, a, b, c
u = a*(1 - 1/(1+(z/c)^2) *  exp(-2*(x^2 + y^2)/(b^2 * (1+(z/c)^2))))
Symbolics.gradient(u, [x,y,z])

