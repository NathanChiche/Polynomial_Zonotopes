using Symbolics
using Nemo
using IntervalArithmetic
using StaticArrays
using Combinatorics
using ForwardDiff

include("polzono_to_symbo.jl")


function gfind_range_derivatives(g,nb_vars,v)
    #g=x[1]^2/4 + (x[2]+1)*(x[3]+2) + (x[3]+3)^2
    #@variables v[1:nb_vars]
    inp = #=@SVector=# [@interval(-1.0,1.0) for i= 1:nb_vars] # CAREFUL, range is evaluated at global scope 
    Dg = Symbolics.jacobian([g], [v[i] for i=1:nb_vars])
    Dg_expr = build_function(Dg, [v[i] for i=1:nb_vars], expression = Val{false})
    my_Dg = eval(Dg_expr[1])
    #return my_Dg
    range_Dg = my_Dg(inp)
    return range_Dg
end


function make_quant2(n1,n2) # C'est nul
    start=["forall" for i in 1:n1+n2]
    res=[]
    index=[i for i in 1:n2]
    for i in 1:n2-1
        q1=start
        q2=vcat(["forall" for i in 1:n1],["exists" for i in 1:n2])
        comb=collect(combinations(index,i))
        #@show(comb)
        #deleteat!(comb, findall(x->x==index,comb))
        #println("apres delete")
        #@show(comb)
        for j in 1:i
            #println(comb[j])
            q1[n1+comb[j][1]]="exists"
            q2[n1+comb[j][1]]="forall"
        end
        push!(res,[q1,q2])
    end
    for k in eachindex(res)
        tmp=res[k]
        push!(res,[tmp[2],tmp[1]])
    end
    return res
end

function mysetdiff(y, x)
    res = Vector{eltype(y)}(undef, length(y) - length(x))
    i = 1
    @inbounds for el in y
        el ∈ x && continue
        res[i] = el
        i += 1
    end
    res
end
collect(3:5)





function gO(range_Jf, i)
    return @interval(-(abs(range_Jf[i]).hi),abs(range_Jf[i]).hi)
end

function gI(range_Jf,i)
    return @interval(-(abs(range_Jf[i]).lo),abs(range_Jf[i]).lo)
end

function gouter_approx(range_Dg,center,quantifiers,varsordering)
    min=center
    max=center
    p=length(quantifiers)
    for i in 0:p-1
        if quantifiers[varsordering[p-i]]=="exists"
            temp=gO(range_Dg,p-i)
            min=min+temp.lo
            max=max+temp.hi
        else
            temp=gI(range_Dg,p-i)
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

function ginner_approx(range_Dg,center,quantifiers,varsordering)
    min=center
    max=center
    p=length(quantifiers)
    for i in 0:p-1
        if quantifiers[varsordering[p-i]]=="exists"
            temp=gI(range_Dg,p-i)
            #@show(p-i,temp)
            min=min+temp.lo
            max=max+temp.hi
        else
            temp=gO(range_Dg,p-i)
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

function gcenters(Polynomes1,Polynomes2,n1,n2)
    @assert length(Polynomes1)==length(Polynomes2)
    z1=zeros(n1)
    z2=zeros(n2)
    centers=[]
    for i in 1:length(Polynomes1)
        push!(centers,Float64(evaluate(Polynomes1[i],z1)-evaluate(Polynomes2[i],z2)))
    end
    return centers
end

function gcall_multiple_outers_bis(rangelist,centers,quantifierslist,listofvarsordering)
    #println("centers :",centers)
    intervals=[]
    #m=length(quantifiers)
    for j in 1:length(rangelist)
        i=gouter_approx(rangelist[j],centers[j],quantifierslist[j],listofvarsordering[j])
        if !(0 in i)
            #println("0 n'est même pas dans l'outer approx")
            return "false"
        end
        push!(intervals,i)
        #push!(ranges_Df,range_Dg)
    end
    return intervals
end



function gcall_multiple_inners_bis(rangelist,centers,quantifierslist,listofvarsordering)
    #println("on entre dans la sous approx multidimensionnelle")
    #println("centers :",centers)
    intervals=[]
    #m=length(x)
    for j in 1:length(rangelist)
        i=ginner_approx(rangelist[j],centers[j],quantifierslist[j],listofvarsordering[j])
        if !(0 in i) || i=="false"
            #println("0 n'est pas dans l'inner pour quant: ",quantifiers)
            return "false"
        end
        push!(intervals,i)
        #push!(ranges_Df,range_Dg)
    end
    return intervals
end

function gwidthofintervalbox(intervalbox,nbvars)
    width=maximum([intervalbox[i].hi-intervalbox[i].lo for i in 1:nbvars])
    return width
end


function gcompute_derivatives_approx(fun,dim,intervalbox)
    ranges=[]
    for a in 1:dim
        push!(ranges,ForwardDiff.gradient(fun[a],intervalbox))
    end
    return ranges
end

function gcompute_ranges(listgradient,dim,intervalbox)
    ranges=[]
    for a in 1:dim
        push!(ranges,listgradient[a](intervalbox))
    end
    return ranges
end


function gzeros_in_intervalvectors(intervals,dim)
    for o in 1:dim
        if !(0 in intervals[o])
            return 0
        end
    end
    return 1
end

function gintervals_centers(intervals)
    m=zeros(Float64,length(intervals))
    for i in 1:length(intervals)
        m[i]=(intervals[i].hi+intervals[i].lo)/2
    end
    return m
end


function gbisect_at_component(interva, component::Int)
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

function gpaverobust(f,dim,n_vars1,n_vars2,epsilon,quantifiers,listofvarsordering)

    
    intervalbox=[interval(-1..1) for c in 1:n_vars1+n_vars2]
    listintervalbox=[intervalbox]
    width=2.0
    tour=1
    #@show(dim,nbvars)
    while width>epsilon
        l=length(listintervalbox)
        for i in 1:l
            #@show(i,l)
            #@show(listintervalbox)
            interval_to_check=popfirst!(listintervalbox)
            #@show(interval_to_check)
            #ranges=gcompute_ranges(listgradie,dim,interval_to_check)
            ranges=gcompute_derivatives_approx(f,dim,interval_to_check)
            @show(ranges)
            #@show(ranges)
            center=[f[b](gintervals_centers(interval_to_check)) for b in 1:dim]
            #@show(center)

            if gcall_multiple_outers_bis(ranges,center,quantifiers,listofvarsordering)==false
                println("pas dans la surapprox")
                return "false"
            else 
                vectinner=gcall_multiple_inners_bis(ranges,center,quantifiers,listofvarsordering)
                @show(vectinner)
                if vectinner=="false"
                    left,right=gbisect_at_component(interval_to_check,tour)
                    push!(listintervalbox,left,right)
                end
            end
            #@show(length(listintervalbox))
            if length(listintervalbox)==0
                return 1
            end
        end
        #@show(listintervalbox)
        tour=(tour)%(n_vars1)+1
        #@show(tour)
        if tour==1
            width=width/2
        end
        #@show(width)
    end
    return "false"
end 


function geometrical_inclusion(FPZ,PZ,epsilon)
    dim=length(PZ.c)
    listfunc,quantif,nv1,nv2=geometricalinclusion_polynomial_zonotopes_to_function(FPZ,PZ)
    for i in 1:length(quantif)
        listofquantifiers=collect(quantif[1][i][1] for i in 1:dim)
        listofvarsordering=collect(quantif[1][i][2] for i in 1:dim)
        if gpaverobust(listfunc,dim,nv1,nv2,epsilon,listofquantifiers,listofvarsordering)==1
            return true
        end
    end
    return false
end



