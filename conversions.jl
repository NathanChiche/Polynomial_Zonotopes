
function get_polynomials_from_SSPZ(PZ::SimpleSparsePolynomialZonotope,field::Nemo.Field)#checked

    # on récupère les polynomes P1,...,Pn issus de la forme PZ={(P1(x1,...xp),...,Pn(x1,...xp)) pour x dans la boule unité pour la distance max}
    #attention il ne faut pas de termes redondants
    c=LazySets.center(PZ)
    G=genmat(PZ)
    E=expmat(PZ)
    nb_vars=size(E)[1]
    anneau,(x)=polynomial_ring(field,nb_vars) 
    #monomials=get_monomials_from_expmat(E,anneau)
    #Polynomes=Array{Any}(undef,size(c)[1])
    Polynomes=[anneau(0) for k in 1:size(c)[1]]
    for i in 1:size(G)[1]
        Polynomes[i]+=c[i]
        for j in 1:size(G)[2]
            setcoeff!(Polynomes[i], E[:,j] , G[i,j])
        end
    end
    return Polynomes
end

function get_SSPZ_from_polynomials(Polynomes)#il faut faire qqchose pour éviter les polynomes de dimesion 1
    """récupérer la forme PZ=<c,G,E> à partir de  PZ={(P1(x1,...xp),...,Pn(x1,...xp)) """

    liste_exp=Vector{Int64}[]
    n=length(Polynomes)
    #deg_list=Array{Int}(undef,n)
    for i in 1:n
        for j in 1:length(Polynomes[i])
            exp=exponent_vector(Polynomes[i],j)
            push!(liste_exp,exp)
        end
    end
    unique!(liste_exp)#on enlève les doublons dans la liste qui viennent du fait que des monomes interviennent dans plusieurs polynomes
    
    nb_vars=length(liste_exp[1])#on veut connaitre le nombre de variables dans l'anneau de polynomes pour les lignes de E 

    zero=zeros(Int64,nb_vars)
    deleteat!(liste_exp,findall(x->x==zero,liste_exp))
    nb_monomes=length(liste_exp)#on compte le nombre de monômes pour la taille de G après avoir enlevé le monôme constant
    #println("nb monomes dans getSSPZ from poly: ",nb_monomes)
    G=zeros(n,nb_monomes)
    #E=zeros([T={Int64}],nb_vars,nb_monomes)
    E=Array{Int64}(undef,nb_vars,nb_monomes)
    c=zeros(n)
    tmp=0
    for i in 1:nb_monomes
        for k in 1:nb_vars
            E[k,i]=liste_exp[i][k] # remplissage de la matrice des exposants
        end
        
        for j in 1:n
            exp=liste_exp[i]
            a=Float64(coeff(Polynomes[j],exp))
            if a==0
                tmp=tmp+1
            end  
            G[j,i]+=a #remplissage de la matrice génératrice 
        end
    end
    #println("nombre de coeffs nuls sur exposants: ",tmp)
    for i in 1:n
        c[i]=Float64(coeff(Polynomes[i],zero)) #remplissage ici du vecteur de centre 
    end

    #println("taille de E apres modifs dans getSSPZ: ",size(E)[2])
    return SimpleSparsePolynomialZonotope(c,G,E)
end

function SimpletoSPZ(SSPZ)
    nb_vars=size(SSPZ.E)[1]
    id=collect(1:nb_vars)
    zer=zeros(Float64,length(SSPZ.c),1)
    return SparsePolynomialZonotope(SSPZ.c,SSPZ.G,zer,SSPZ.E,id)
end

function SPZ_to_SimpleSPZ(PZ::SparsePolynomialZonotope)
    dim,nb_indep=size(PZ.GI)
    @assert dim==length(PZ.c) "erreur delogique avec la dimension"
    z=zeros(Float64,dim,nb_indep)
    if nb_indep==0 || z==PZ.GI
        return SimpleSparsePolynomialZonotope(PZ.c,PZ.G,PZ.E)
    end
    println("Mauvaise partie CONVERSION GI NON VIDE")
    Gnew=hcat(PZ.G,PZ.GI)
    n,m=size(PZ.E)
    Z=zeros(Int64,nb_indep,m)
    Enew=vcat(PZ.E,Z)
    
    v=zeros(Int64,nb_indep+n)
    v_temp=zeros(Int64,nb_indep+n)

    for j in 1:nb_indep
        v_temp=copy(v)
        v_temp[n+j]=1
        Enew=hcat(Enew,v_temp)
    end
    return SimpleSparsePolynomialZonotope(PZ.c,Gnew,Enew)
end


function sparse_poly_zono_to_dynamic(z)
    d = size(z.E, 1)
    @polyvar x[1:d]  # Create polynomial variables x[1], x[2], ..., x[d]
    n = length(z.c)

    num_generators = size(z.G, 2)
    poly_list = Vector{Any}(undef, n)
    
    for i in 1:n
        poly = z.c[i]
        for j in 1:num_generators
            # Construct the monomial: x_1^(E[j,1]) * x_2^(E[j,2]) * ... * x_d^(E[j,d])
            monomial = prod(x[k]^z.E[k,j] for k in 1:d)
            poly += z.G[i, j] * monomial
        end
        poly_list[i] = poly
    end
    
    return poly_list, x
end


function polynomial_transformation(poly_list, T)
    # Get the union of independent variables from all transformation polynomials.
    y_vars_set = Set{Any}()
    for Ti in T
        union!(y_vars_set, variables(Ti))
    end
    # Sort them to fix an order.
    #@show(y_vars_set)
    y_vars_ordered = sort(collect(y_vars_set), by=string)
    #@show(y_vars_ordered)
    
    if length(y_vars_ordered) != length(poly_list)
        error("The number of independent variables in the transformation ($(length(y_vars_ordered))) must match the dimension of the zonotope ($(length(poly_list))).")
    end
    
    # Build substitution pairs: each variable is replaced by the corresponding zonotope polynomial.
    sub_pairs = [y_vars_ordered[i] => poly_list[i] for i in 1:length(y_vars_ordered)]
    
    # Apply the substitution using splatting with subs.
    transformed_zonotope = [subs(Ti, sub_pairs...) for Ti in T]
    return transformed_zonotope
end

"""using MultivariatePolynomials
@polyvar x1 x2
p1 = 1 + x1       # p₁(x) = 1 + x
p2 = 2 + 2*x2     # p₂(x) = 2 + 2x
poly_list = [p1, p2]
@polyvar y1 y2
T1 = 1    # T₁(y₁, y₂) = y₁² + y₂
T2 = y1 + y2^2   # T₂(y₁, y₂) = y₁ + y₂²
T = [T1, T2]
transformed = polynomial_transformation(poly_list, T)"""


function get_exponent(m, v)
    # Get the variables that appear in m.
    vs = variables(m)
    # Get the full exponent vector corresponding to those variables.
    exps = exponents(m)
    # Find the position of v in vs.
    idx = findfirst(x -> x == v, vs)
    return idx === nothing ? 0 : exps[idx]
end


function dynamic_to_sparse_poly_zono(P_list)
    n = length(P_list)
    # Collect the union of variables from all polynomials.
    vars_set = Set{Any}()
    for p in P_list
        union!(vars_set, variables(p))
    end
    # Sort the variables (using their string representation for consistency).
    vars = sort(collect(vars_set), by=string)
    d = length(vars)
    
    # Prepare the center vector and a temporary storage for nonconstant coefficients.
    c = Vector{Float64}(undef, n)
    poly_coeffs = Vector{Dict{NTuple{d,Int}, Number}}(undef, n)
    for i in 1:n
        poly_coeffs[i] = Dict{NTuple{d,Int}, Number}()
    end
    
    # Dictionary to collect the union of nonconstant monomials (represented as exponent tuples).
    monomial_to_index = Dict{NTuple{d,Int}, Int}()
    
    # Process each polynomial.
    for i in 1:n
        p = P_list[i]
        mons = MultivariatePolynomials.monomials(p)
        coeffs = MultivariatePolynomials.coefficients(p)
        constant_found = false
        for (m, a) in zip(mons, coeffs)
            # Compute the total degree of monomial m over all variables in our sorted union.
            deg = sum(get_exponent(m, v) for v in vars)
            if deg == 0
                c[i] = a
                constant_found = true
            else
                # Build the exponent tuple for m with respect to our sorted order.
                exp_tuple = ntuple(j -> get_exponent(m, vars[j]), d)
                poly_coeffs[i][exp_tuple] = a
                # Add this monomial (if not already present) to the union.
                if !haskey(monomial_to_index, exp_tuple)
                    monomial_to_index[exp_tuple] = 0  # placeholder; will assign index later.
                end
            end
        end
        # If no constant term was found, set it to 0.
        if !constant_found
            c[i] = 0
        end
    end
    
    # Create a sorted list of the unique nonconstant monomials (their exponent tuples).
    gen_keys = sort(collect(keys(monomial_to_index)))
    m = length(gen_keys)
    # Assign each unique monomial an index (column in the generator matrix).
    for (j, key) in enumerate(gen_keys)
        monomial_to_index[key] = j
    end
    
    # Build the generator matrix G (n×m).
    G = zeros(Float64, n, m)
    for i in 1:n
        for (exp_tuple, a) in poly_coeffs[i]
            j = monomial_to_index[exp_tuple]
            G[i, j] = a
        end
    end
    
    # Build the exponent matrix E (d×m): each column is the exponent tuple of a generator.
    E = zeros(Int, d, m)
    for (j, exp_tuple) in enumerate(gen_keys)
        for k in 1:d
            E[k, j] = exp_tuple[k]
        end
    end
    
    return SimpleSparsePolynomialZonotope(c,G,E)
end


"""PT=SimpleSparsePolynomialZonotope([0.0 , 0],[1 0 2; 0 1.0 1],[2 0 1; 0 1 1])
#pol=Array{DynamicPolynomials.Polynomial{DynamicPolynomials.Commutative{DynamicPolynomials.CreationOrder}, Graded{LexOrder}, Int64}}(undef,2)
@polyvar z1 z2
pol=z1+z1*z2
MultivariatePolynomials.monomials(pol)
PT=SimpleSparsePolynomialZonotope([1.0 , 0],[1 3 2; 0 1.0 1],[2 0 1; 0 1 1; 0 1 0])
sparse_poly_zono_to_dynamic(PT)[1]
sparse_poly_zono_to_dynamic((dynamic_to_sparse_poly_zono(sparse_poly_zono_to_dynamic(PT)[1])))[1]"""

