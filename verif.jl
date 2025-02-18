using Nemo
using LazySets


include("testinclusion_robrange.jl")
include("inclusiongeom_robustrange.jl")
include("conversions.jl")
include("joins.jl")
include("plotsample.jl")
include("polynomap.jl")
include("matlabmatrix.jl")
include("reduction.jl")

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
