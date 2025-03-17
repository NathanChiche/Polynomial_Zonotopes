using LazySets
using Symbolics
using Nemo
using Combinatorics

include("conversions.jl")

function quant_and_varsorders(n1,n2)
    indexbase=collect(1:n1)
    nt=n1+n2
    index=collect(n1+1:nt)
    base=["forall" for j in 1:n1]
    res=[]
    for i in 1:n2-1
        q1=vcat(base,["forall" for j in 1:i],["exists" for j in 1:n2-i])
        q2=vcat(base,["forall" for j in 1:n2-i],["exists" for j in 1:i])
        #@show(q1)
        comb=collect(combinations(index,i))
        #@show(comb)
        #diff=setdiff(index,comb)
        
        for c in comb
            #@show(c)
            diff=setdiff(index,c)
            #@show(diff)
            ordersq1=vcat(indexbase,c,diff)
            ordersq2=vcat(indexbase,diff,c)
            push!(res,[[q1,ordersq1],[q2,ordersq2]])
        end
    end
    return res
end



function functionalinclusion_polynomial_zonotopes_to_function(fpzonotope,pzonotope)
    """f(pZ) à gauche absolument"""
    fcenter = fpzonotope.c  # Center of the zonotope f(PZ)
    fgenerators = fpzonotope.G  # Polynomial generators of the zonotope f(PZ)
    fexponents = fpzonotope.E # Exponents of the zonotope f(PZ)

    # Step 2: Extract center and generators
    center = pzonotope.c  # Center of the zonotope
    generators = pzonotope.G  # Polynomial generators
    exponents = pzonotope.E # Exponents
    dim=size(generators)[1]
    # Step 3: Define symbolic variables corresponding to each variables in the zonotope
    n_generators = size(generators)[2]
    nf_generators=size(fgenerators)[2]
    n_vars=size(exponents)[1]
    
    #create variables
    @variables var[1:n_vars+dim]
    y=[var...]

    # Step 4: Construct symbolic polynomial expression
    # Initialize with the center
    symbolic_expr = []
    
    # Add terms for each generator
    for i in 1:dim
        temp_symbol_expr=fcenter[i]
        for j in 1:nf_generators
            poly_term = fgenerators[i,j] * prod(var[k]^fexponents[k,j] for k in 1:n_vars)
            temp_symbol_expr+=poly_term
        end
        push!(symbolic_expr,temp_symbol_expr)
    end

    #@show(n_vars-dim)
    for i in 1:dim
        symbolic_expr[i] =  symbolic_expr[i] - center[i]
        for j in 1:n_generators
            symbolic_expr[i]= symbolic_expr[i] - (generators[i,j] * prod(var[k]^exponents[k,j] for k in 1:n_vars-dim)*prod(var[k+dim]^exponents[k,j] for k in n_vars-dim+1:n_vars))
        end
    end

    #@show(symbolic_expr)

    quantifiers=vcat(["forall" for i in 1:n_vars],["exists" for i in 1:dim])
    listfunctions=[]
    listgradient=[]
    #listg=[]
    for i in 1:dim
        f=build_function(symbolic_expr[i],[var[j] for j=1:n_vars+dim],expression=Val{false})
        push!(listfunctions,f)
        gradfsymbo=Symbolics.gradient(symbolic_expr[i],y)
        #@show(gradfsymbo)
        #push!(listg,gradfsymbo)
        gradfun=build_function(gradfsymbo,[y[j] for j=1:n_vars+dim],expression=Val{false})[1]
        push!(listgradient,gradfun)

    end

    #@show(symbolic_expr)
    #@show(listg)
    #return listg
    return listfunctions,listgradient,quantifiers,n_vars+dim#=,listg=#
    return symbolic_expr
end

function geometricalinclusion_polynomial_zonotopes_to_function(fpzonotope,pzonotope)
    """f(pZ) à gauche absolument"""
    fcenter = fpzonotope.c  # Center of the zonotope f(PZ)
    fgenerators = fpzonotope.G  # Polynomial generators of the zonotope f(PZ)
    fexponents = fpzonotope.E # Exponents of the zonotope f(PZ)

    # Step 2: Extract center and generators
    center = pzonotope.c  # Center of the zonotope
    generators = pzonotope.G  # Polynomial generators
    exponents = pzonotope.E # Exponents
    dim=size(generators)[1]
    # Step 3: Define symbolic variables corresponding to each variables in the zonotope
    n_generators = size(generators)[2]
    nf_generators=size(fgenerators)[2]
    n_vars2=size(exponents)[1]
    n_vars1=size(fexponents)[1]
    @show(n_vars1,n_vars2)
    
    #create variables
    @variables var[1:n_vars1+n_vars2]
    y=[var...]

    # Step 4: Construct symbolic polynomial expression
    # Initialize with the center
    symbolic_expr = []
    
    # Add terms for each generator
    for i in 1:dim
        temp_symbol_expr=fcenter[i]-center[i]
        for j in 1:nf_generators
            poly_term = fgenerators[i,j] * prod(var[k]^fexponents[k,j] for k in 1:n_vars1)
            temp_symbol_expr+=poly_term
        end
        for j in 1:n_generators
            poly_term = generators[i,j] * prod(var[k+n_vars1]^exponents[k,j] for k in 1:n_vars2)
            temp_symbol_expr-=poly_term
        end
        push!(symbolic_expr,temp_symbol_expr)
    end

    #quantifiers=vcat(["forall" for i in 1:n_vars1],["exists" for i in 1:dim])

    quantifier=quant_and_varsorders(n_vars1,n_vars2)
    listfunctions=[]
    listgradient=[]
    for i in 1:dim
        f=build_function(symbolic_expr[i],[var[j] for j=1:n_vars1+n_vars2],expression=Val{false})
        push!(listfunctions,f)
        gradfsymbo=Symbolics.gradient(symbolic_expr[i],y)
        gradfun=build_function(gradfsymbo,[y[j] for j=1:n_vars1+n_vars2],expression=Val{false})[1]
        push!(listgradient,gradfun)
    end

    #@show(symbolic_expr)
    return listfunctions,listgradient,quantifier,n_vars1,n_vars2
    return symbolic_expr
end

#res= polynomial_zonotopes_to_function(fP,P)


# Example Usage:
# Define a Polynomial Zonotope
##SP = SimpleSparsePolynomialZonotope([2.0, 0], [1 0;0 2.], [1 4;1 2])

"""R=RealField()
S,(x,y,s,t)=polynomial_ring(R,["x","y","s","t"])

P=get_SSPZ_from_polynomials([x+t,y+s])
fP=get_SSPZ_from_polynomials([x^2+x^3+s+t,y^2+s])
fP.E
fP.c 

geometricalinclusion_polynomial_zonotopes_to_function(fP,P)
pozon=SimpleSparsePolynomialZonotope([0 , 0.0],[1 0.0; 0 1],[1 0; 0 0; 1 0; 1 0; 0 1; 1 1])"""

# Convert to symbolic representation"""

