function reduce_order(P::SimpleSparsePolynomialZonotope, r::Real,
    method::AbstractReductionMethod=GIR05())
@assert r ≥ 1
n = dim(P)
h = ngens_dep(P)
q = ngens_indep(P)

c = center(P)
G = genmat_dep(P)
GI = genmat_indep(P)
E = expmat(P)
idx = indexvector(P)

a = max(0, min(h + q, ceil(Int, h + q - n * (r - 1))))
Gbar = hcat(G, GI)
norms = [norm(g) for g in eachcol(Gbar)]
th = sort(norms)[a]

# TODO: case a = 0
# TODO is constructing an array of booleans the most efficient way?
K = [norms[i] ≤ th for i in 1:h]
Kbar = .!K

H = [norms[h + i] ≤ th for i in 1:q]
Hbar = .!H

PZ = SparsePolynomialZonotope(c, G[:, K], GI[:, H], E[:, K], idx)
Z = reduce_order(overapproximate(PZ, Zonotope), 1, method)

Ebar = E[:, Kbar]
N = [!iszero(e) for e in eachrow(Ebar)]

cz = center(Z)
Gz = genmat(Z)
return SparsePolynomialZonotope(cz, G[:, Kbar], hcat(GI[:, Hbar], Gz),
                  Ebar[N, :], idx[N])
end