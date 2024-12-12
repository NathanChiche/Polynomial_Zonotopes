using LazySets
using Symbolics
using Nemo

include("conversions.jl")

function functionalinclusion_polynomial_zonotopes_to_function(fpzonotope::PolynomialZonotope,pzonotope::PolynomialZonotope)
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
    y=[vars...]

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

    quantifiers=vcat(["forall" for i in 1:n_vars],["exists" for i in 1:dim])
    listfunctions=[]
    listgradient=[]
    for i in 1:dim
        f=build_function(symbolic_expr[i],[var[j] for j=1:n_vars+dim],expression=Val{false})
        push!(listfunctions,f)
        gradfsymbo=Symbolics.gradient(symbolic_expr[i],y)
        gradfun=build_function(gradfsymbo,[y[j] for j=1:n_vars+dim],expression=Val{false})[1]
        push!(listgradient,gradfun)

    end
    @show(symbolic_expr)

    return listfunctions,listgradient,quantifiers,n_vars+dim
    return symbolic_expr
end

function geometricalinclusion_polynomial_zonotopes_to_function(fpzonotope::PolynomialZonotope,pzonotope::PolynomialZonotope)
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

    # Step 4: Construct symbolic polynomial expression
    # Initialize with the center
    symbolic_expr = []
    
    # Add terms for each generator
    for i in 1:dim
        temp_symbol_expr=fcenter[i]
        for j in 1:nf_generators
            poly_term = fgenerators[i,j] * prod(var[k]^fexponents[k,j] for k in 1:n_vars1)
            temp_symbol_expr+=poly_term
        end
        push!(symbolic_expr,temp_symbol_expr)
    end

    #@show(n_vars-dim)
    for i in 1:dim
        symbolic_expr[i] =  symbolic_expr[i] - center[i]
        for j in 1:n_generators 
            symbolic_expr[i]= symbolic_expr[i] - (generators[i,j] *prod(var[k+n_vars1]^exponents[k,j] for k in 1:n_vars2))
        end
    end

    #quantifiers=vcat(["forall" for i in 1:n_vars1],["exists" for i in 1:dim])

    quantifier=quant_and_varsorders(n_vars1,n_vars2)
    listfunctions=[]
    for i in 1:dim
        f=build_function(symbolic_expr[i],[var[j] for j=1:n_vars1+n_vars2],expression=Val{false})
        push!(listfunctions,f)
    end
    @show(symbolic_expr)
    return listfunctions,quantifier,n_vars1,n_vars2
    return symbolic_expr
end

#res= polynomial_zonotopes_to_function(fP,P)


# Example Usage:
# Define a Polynomial Zonotope
##SP = SimpleSparsePolynomialZonotope([2.0, 0], [1 0;0 2.], [1 4;1 2])

"""R=RealField()
S,(x,y,s,t)=PolynomialRing(R,["x","y","s","t"])

P=get_SSPZ_from_polynomials([x+t,y+s])
fP=get_SSPZ_from_polynomials([x^2+x^3+s+t,y^2+s])
fP.E
fP.c 

geometricalinclusion_polynomial_zonotopes_to_function(fP,P)
pozon=SimpleSparsePolynomialZonotope([0 , 0.0],[1 0.0; 0 1],[1 0; 0 0; 1 0; 1 0; 0 1; 1 1])"""

# Convert to symbolic representation"""

