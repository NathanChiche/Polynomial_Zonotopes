using LazySets
using Symbolics
using Nemo

include("conversions.jl")

function polynomial_zonotopes_to_function(fpzonotope::PolynomialZonotope,pzonotope::PolynomialZonotope)
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


    for i in 1:dim
        symbolic_expr[i] =  symbolic_expr[i] - center[i]
        for j in 1:n_generators
            symbolic_expr[i]= symbolic_expr[i] - (generators[i,j] * prod(var[k]^exponents[k,j] for k in 1:n_vars-dim)*prod(var[k+dim]^exponents[k,j] for k in n_vars-dim+1:n_vars))
        end
    end

    quantifiers=vcat(["forall" for i in 1:n_vars],["exists" for i in 1:dim])
    listfunctions=[]
    for i in 1:dim
        f=build_function(symbolic_expr[i],[var[j] for j=1:n_vars+dim],expression=Val{false})
        push!(listfunctions,f)
    end
    @show(symbolic_expr)
    return listfunctions,quantifiers,n_vars+dim
    return symbolic_expr
end

#res= polynomial_zonotopes_to_function(fP,P)


# Example Usage:
# Define a Polynomial Zonotope
##SP = SimpleSparsePolynomialZonotope([2.0, 0], [1 0;0 2.], [1 4;1 2])

"""R=RealField()
S,(x,y,s,t)=PolynomialRing(R,["x","y","s","t"])

P=get_SSPZ_from_polynomials([1+x+t,y+s])
fP=get_SSPZ_from_polynomials([x^2+x^3+s+t,y^2+s])
fP.E
fP.c 
SP.G
@show(SP)

prod(2 for i in 1:3)*prod(2 for j in 1+1:5-1)
nv=8
@variables xva[1:nv]
pozon=SimpleSparsePolynomialZonotope([0 , 0.0],[1 0.0; 0 1],[1 0; 0 0; 1 0; 1 0; 0 1; 1 1])
prod(xva[k]^pozon.E[k,2] for k in 1:4)
prod(xva[k+2]^pozon.E[k,1] for k in 5:6)
# Convert to symbolic representation"""

