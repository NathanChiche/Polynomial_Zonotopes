using Symbolics
using Nemo
using IntervalArithmetic
using StaticArrays
using Combinatorics
using ForwardDiff


function find_range_derivatives(g,nb_vars,v)
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

function nemo_to_symbo(polynomes1,x,start_index)
    #nb_vars=length(gens(parent(p1)))*2
    #@variables x[1:nb_vars]
    nb_vars=length(x)
    n=length(polynomes1)
    f=[]
    for i in 1:n
        p=polynomes1[i]
        exponents1=collect(exponent_vector(p,i) for i in 1:length(p))
        terms=[]
        for e in exponents1
            term=1
            for j in eachindex(e)
                term=term*x[j+start_index]^e[j]
            end
            term=Float64(coeff(p,e))*term
            push!(terms,term)
        end
        funct=sum(terms)
        push!(f,funct)
    end
    return f
end

function make_quant2(n1,n2) # C'est nul
    start=["forall" for i in 1:n1+n2]
    res=[]
    index=[i for i in 1:n2]
    for i in 1:n2-1
        q1=start
        q2=vcat(["forall" for i in 1:n1],["exists" for i in 1:n2])
        comb=collect(combinations(index,i))
        deleteat!(comb, findall(x->x==index,comb))
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
make_quant2(1,2)


function O(range_Jf, i)
    return @interval(-(abs(range_Jf[i]).hi),abs(range_Jf[i]).hi)
end

function I(range_Jf,i)
    return @interval(-(abs(range_Jf[i]).lo),abs(range_Jf[i]).lo)
end

function outer_approx(range_Dg,center,quantifiers,variablesorder)
    min=center
    max=center
    p=length(quantifiers)
    for i in 0:p-1

        if quantifiers[variablesorder[p-i]]=="exists"
            temp=O(range_Dg,variablesorder[p-i])
            println(temp ,"exists")
            min=min+temp.lo
            max=max+temp.hi
        else
            temp=I(range_Dg,variablesorder[p-i])
            println(temp ,"forall")
            min=min+temp.hi
            max=max+temp.lo
            if max<min
                println("resultat de sur-approx vide")
                return false
            end
        end  
    end
    #println(min,typeof(min))
    #println(max,typeof(max))
    return @interval(min,max)
end

function inner_approx(range_Dg,center,quantifiers,variablesorder)
    min=center
    max=center
    p=length(quantifiers)
    for i in 0:p-1
        if quantifiers[variablesorder[p-i]]=="exists"
            temp=I(range_Dg,variablesorder[p-i])
            min=min+temp.lo
            max=max+temp.hi
        else
            temp=O(range_Dg,variablesorder[p-i])
            min=min+temp.hi
            max=max+temp.lo
            if max<min
                println("resultat de sous-approx vide")
                return false
            end
        end  
    end
    return @interval(min,max)
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

function call_multiple_outers_bis(rangelist,centers,quantifiers,variablesorder)
    println("centers :",centers)
    intervals=[]
    #m=length(quantifiers)
    for j in 1:length(rangelist)
        i=outer_approx(rangelist[j],centers[j],quantifiers,variablesorder)
        if !(0 in i)
            println("0 n'est même pas dans l'outer approx")
            return false
        end
        push!(intervals,i)
        #push!(ranges_Df,range_Dg)
    end
    return intervals
end

function call_multiple_inners_bis(rangelist,centers,quantifiers,variablesorder)
    println("centers :",centers)
    intervals=[]
    #m=length(x)
    for j in 1:length(rangelist)
        i=inner_approx(rangelist[j],centers[j],quantifiers,variablesorder)
        if !(0 in i) || i==false
            println("0 n'est pas dans l'inner pour quant: ",quantifiers)
            return false
        end
        push!(intervals,i)
        #push!(ranges_Df,range_Dg)
    end
    return intervals
end

function call_multiple_inners(nD_funct,centers,quantifiers,input,x)
    println("centers :",centers)
    intervals=[]
    m=length(x)
    for j in 1:length(nD_funct)
        println("fonction sur laquelle on teste: ",nD_funct[j])
        Dg = Symbolics.jacobian([nD_funct[j]], [x[i] for i=1:m])
        println("sa jacobienne: ",Dg)
        Dg_expr = build_function(Dg, [x[i] for i=1:m],expression=Val{false})
        my_Dg = eval(Dg_expr[1])
        range_Dg = my_Dg(input) 
        println("ici la range de la jaco",range_Dg)
        #println(range_Dg)
        i=inner_approx(range_Dg,centers[j],quantifiers[j])
        if !(0 in i) || i==false
            println("0 n'est pas dans l'inner pour quant: ",quantifiers)
            return false
        end
        push!(intervals,i)
        #push!(ranges_Df,range_Dg)
    end
    return intervals
end

function inclusion_test(Polynomes1,Polynomes2)
    dim=length(Polynomes1)
    n1=length(gens(parent(Polynomes1[1])))
    n2=length(gens(parent(Polynomes2[1])))
    quantifiers=vcat(["forall" for i in 1:n1],["exists" for i in 1:n2])
    nb_vars=n1+n2

    #input = @SVector [@interval(-1.0,1.0) for i= 1:nb_vars]   
    inte=@interval(-1.0,1.0)
    input=[inte for i in 1:nb_vars] 

    @variables x[1:nb_vars]
    f=nemo_to_symbo(Polynomes1,x,0)
    g=nemo_to_symbo(Polynomes2,x,n1)
    symb_fct=[]
    for j in 1:dim
        push!(symb_fct,f[j]-g[j])
    end
    #println(symb_fct)
    #built_functions=[]
    ranges=[]
    for l in 1:dim
        g=build_function(symb_fct[l],[x[i] for i=1:nb_vars],expression=Val{false})
        my_g=eval(g)
        #push!(built_functions,my_g)
        a=ForwardDiff.gradient(my_g,input)
        b=find_range_derivatives(symb_fct[l],nb_vars,x)
        println(a)
        println(b)
        println(a==b)
        push!(ranges,ForwardDiff.gradient(my_g,input))
    end
    #return ranges
    ctrs=centers(Polynomes1,Polynomes2,n1,n2)
    if call_multiple_outers_bis(ranges,ctrs,quantifiers) ==false
        println("on ne peut pas avoir d'inclusion à cause de la surapprox")
        return false
    end
    list_quantif=make_quant2(n1,n2)
    println("on passe aux inners")
    for quant in list_quantif
        #println(quant)
        if call_multiple_inners_bis(ranges,ctrs,quant)!=false
            println(quant,"inclusion")
            return true
        end
    end
    return false
    #génération de la liste des exists/forall
end

function test_primal() #ATTENTION: plutôt gros problème sur la tête des ranges jacobiennes qui doivent être des intervalles
    R=RealField()
    Anneau,(x,y)=PolynomialRing(R,["x","y"])
    pol1=[2*x^2+2*x]
    pol2=[3*x^2+5*x]
    #pol1=[2*x^2+2*x,2*y^3+2*y]
    #pol2=[3*x^2+5*x,3*(y^3)+3*y] #ici c'est censé etre un true dans les bons cas
    #pol1=[2*x*y^2+y+3,y^2*x^3+1]
    #pol2=[x^2+x,2*y^2*x+1]
    #println(centers(pol1,pol2,2,2)) celle-ci fonctionne pas mal
    inclusion_test(pol1,pol2)
end

test_primal()
f(x)=(1//4)*(x[1]^2) + (1 + x[2])*(2 + x[3]) + (3 + x[3])^2

function defineintricate()
    
    f(x)=(1//4)*(x[1]^2) + (1 + x[2])*(2 + x[3]) + (3 + x[3])^2
    return f
end

h=defineintricate()
ForwardDiff.gradient(h,input)

function paverobust(f,dim,vars,epsilon)
    inte=IntervalBox(-1..1,dim)
    listintervals=[inte]
    width=2
    tour=1

    while width>epsilon
        l=length(listintervals)
        interval_to_check=popfirst!(listintervals)
        for i in 1:l
            rangejf=ForwardDiff.gradient(f,interval_to_check)
            if outer_approx(rangejf,center,quantifiers)
                return 0
            else !(0 in inner_approx(rangejf,center,quantifiers))
                left,right=bisect_at_component(interval_to_check,tour)
                push!(listintervals,left,right)
            end
        end
        if length(listintervals)==0
            return 1
        end
        tour=(tour+1)%dim
        width=last(listintervals).hi-last(listintervals).lo
    end
    return 0
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

