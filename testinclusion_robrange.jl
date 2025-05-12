using Symbolics
using Nemo
using IntervalArithmetic
#using StaticArrays
using Combinatorics
using ForwardDiff

include("polzono_to_symbo.jl")

function outer_approx(absrange_Dg,center,quantifiers,interval_to_check)
    min=center
    max=center
    p=length(quantifiers)
    for i in 0:p-1
        if quantifiers[p-i]=="exists"

            temp=(absrange_Dg[p-i].hi)*(interval_to_check[p-i].hi-interval_to_check[p-i].lo)/2
            #@show(p-i,temp)
            min=min-temp
            max=max+temp
        else
            temp=(absrange_Dg[p-i].lo)*(interval_to_check[p-i].hi-interval_to_check[p-i].lo)/2
            #@show(p-i,temp)
            min=min+temp
            max=max-temp
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

function inner_approx(absrange_Dg,center,quantifiers,interval_to_check)
    min=center
    max=center
    p=length(quantifiers)
    #@show(interval_to_check)
    for i in 0:p-1
        if quantifiers[p-i]=="exists"
            #println("exists")
            #@show(absrange_Dg[p-i].lo)
            temp=(absrange_Dg[p-i].lo)*(interval_to_check[p-i].hi-interval_to_check[p-i].lo)/2
            #@show(p-i,temp)
            min=min-temp
            max=max+temp
        else
            #println("forall")
            #@show(absrange_Dg[p-i].hi)
            temp=(absrange_Dg[p-i].hi)*(interval_to_check[p-i].hi-interval_to_check[p-i].lo)/2
            #@show(p-i,temp)
            min=min+temp
            max=max-temp
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

function call_multiple_outers_bis(rangelist,centers,quantifiers,interval_to_check)
    #@show(rangelist)
    #println("centers :",centers)
    intervals=[]
    #m=length(quantifiers)
    for j in 1:length(rangelist)
        i=outer_approx(rangelist[j],centers[j],quantifiers,interval_to_check)
        if !(0 in i)
            #println("0 n'est même pas dans l'outer approx")
            return "false"
        end
        push!(intervals,i)
        #push!(ranges_Df,range_Dg)
    end
    return intervals
end



function call_multiple_inners_bis(rangelist,centers,quantifiers,interval_to_check)
    #println("on entre dans la sous approx multidimensionnelle")
    #println("centers :",centers)
    intervals=[]
    #m=length(x)
    for j in 1:length(rangelist)
        i=inner_approx(rangelist[j],centers[j],quantifiers,interval_to_check)
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


function compute_ranges(gradient,intervalbox,quantifiers,intervalboxcenters)
    n=length(quantifiers)
    ranges=Array{IntervalArithmetic.Interval{Float64}}(undef,n)
    centered_interval=copy(intervalbox)
    for k in 1:n
        if quantifiers[k]=="exists"
            if typeof(intervalboxcenters[k])!= IntervalArithmetic.Interval
                centered_interval[k]=intervalboxcenters[k]..intervalboxcenters[k]
            else
                centered_interval[k]=intervalboxcenters[k]
            end
        end
    end
    #NOW centered_interval equals our intervalbox where all existentially quantified variables are "reduced to their center"
    for k in 1:n
        if quantifiers[k]=="exists"
            if typeof(abs.(gradient(intervalbox)[k]))!= IntervalArithmetic.Interval{Float64} && typeof(abs.(gradient(centered_interval)[k]))!= IntervalArithmetic.Interval{Int64}
                ranges[k]=abs.(gradient(intervalbox)[k])..abs.(gradient(intervalbox)[k])
            else
                ranges[k]=abs.(gradient(intervalbox)[k])
            end
        else    
            if typeof(abs.(gradient(centered_interval)[k]))!= IntervalArithmetic.Interval{Float64} && typeof(abs.(gradient(centered_interval)[k]))!= IntervalArithmetic.Interval{Int64}
                #@show(abs.(gradient(centered_interval)[k]))
                #@show(typeof(abs.(gradient(centered_interval)[k])))
                ranges[k]=abs.(gradient(centered_interval)[k])..abs.(gradient(centered_interval)[k])
            else 
                ranges[k]=abs.(gradient(centered_interval)[k]) # we apply strictly the theorem 2 form goubault putot 2020 which states we can evaluate in the center of the existentially quantified variables
            end
        end
    end
    return ranges
end

IntervalArithmetic.Interval{Float64}==IntervalArithmetic.Interval

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

function linearfunctionalinclusion(f,listgradie,dim,nbvars,epsilon,quantifiers,nbinitialvariables)
    interval_to_check=[interval(-1..1) for c in 1:nbvars]
    intervcenters=intervals_centers(interval_to_check)
    rangelist=[compute_ranges(listgradie[k],interval_to_check,quantifiers,intervcenters) for k in 1:dim]
    center=[f[b](intervals_centers(interval_to_check)) for b in 1:dim]
    if call_multiple_outers_bis(rangelist,center,quantifiers,interval_to_check)==false
        println("pas dans la surapprox")
        return "false"
    else 
        return call_multiple_inners_bis(rangelist,center,quantifiers,interval_to_check)
    end
end

function paverobust(f,listgradie,dim,nbvars,epsilon,quantifiers,nbinitialvariables,linearsystem)
    println("entre pavage")
    if linearsystem
        return linearfunctionalinclusion(f,listgradie,dim,nbvars,epsilon,quantifiers,nbinitialvariables)
    end
    intervalbox=[interval(-1..1) for c in 1:nbvars]
    listintervalbox=[intervalbox]
    width=2.0
    tour=1
    #@show(dim,nbvars)
    while width>epsilon
        l=length(listintervalbox)

        #@show(l)
        for i in 1:l
            println(i)
            interval_to_check=popfirst!(listintervalbox)
            #@show(interval_to_check)
            #ranges=compute_ranges(listgradie,dim,interval_to_check)
            intervcenters=intervals_centers(interval_to_check)
            rangelist=[compute_ranges(listgradie[k],interval_to_check,quantifiers,intervcenters) for k in 1:dim]
            #@show(rangelist)
            center=[f[b](intervals_centers(interval_to_check)) for b in 1:dim]
            #@show(center)

            if call_multiple_outers_bis(rangelist,center,quantifiers,interval_to_check)==false
                println("pas dans la surapprox")
                return "false"
            else 
                vectinner=call_multiple_inners_bis(rangelist,center,quantifiers,interval_to_check)
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
        @show(tour)
        tour=(tour)%(nbinitialvariables)+1
        if tour==1
            width=width/2
        end

        #@show(width)
    end
    println("sort pavage")
    return "false"
end 


function inclusion_test(FPZ,PZ,epsilon,nbinitialvariables,linearsystem)
    dim=length(PZ.c)
    println("entre traduction ")
    listfunc,listgrad,quantifiers,nbvars=functionalinclusion_polynomial_zonotopes_to_function(FPZ,PZ)
    println("sort traduction")
    inte=interval(-1..1)
    #intervalbox=[inte for i in 1:nbvars]
    
    #@show(listgrad[1](intervalbox))
    #@show(listgrad[2](intervalbox))
    return paverobust(listfunc,listgrad,dim,nbvars,epsilon,quantifiers,nbinitialvariables,linearsystem)

end


function test_primal() #ATTENTION: plutôt gros problème sur la tête des ranges jacobiennes qui doivent être des intervalles
    R=RealField()
    Anneau,(x,y)=polynomial_ring(R,["x","y"])
    pol1=[2*x^2+2*x]
    pol2=[3*x^2+5*x]
    inclusion_test(pol1,pol2)
end

