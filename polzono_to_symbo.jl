using LazySets
using Symbolics

function polynomial_zonotope_to_symbolic(pzonotope::PolynomialZonotope,var)
    # Step 2: Extract center and generators
    center = pzonotope.c  # Center of the zonotope
    generators = pzonotope.G  # Polynomial generators
    exponents = pzonotope.E # Exponents
    dim=size(generators)[1]
    # Step 3: Define symbolic variables corresponding to each variables in the zonotope
    n_generators = size(generators)[2]
    n_vars=size(exponents)[1]
    

    # Step 4: Construct symbolic polynomial expression
    # Initialize with the center
    symbolic_expr = []
    
    # Add terms for each generator
    for i in 1:dim
        @show(exponents)
        temp_symbol_expr=center[i]
        for j in 1:n_generators
            #@show(exponents[j])
            #@show(generators[i][j])
            #@show(var[k] for k in 1:n_vars)
            poly_term = generators[i,j] * prod(var[k]^exponents[k,j] for k in 1:n_vars)
            temp_symbol_expr+=poly_term
        end
        push!(symbolic_expr,temp_symbol_expr)
    end
    
    return symbolic_expr
end


# Example Usage:
# Define a Polynomial Zonotope
S = SimpleSparsePolynomialZonotope([2.0, 0], [1 0;0 2.], [1 4;1 2])
E=S.E
E[2,1]
@variables z[1:2]

# Convert to symbolic representation
symbolic_expression = polynomial_zonotope_to_symbolic(S,z)


# Print the resulting symbolic function
println(symbolic_expression)
