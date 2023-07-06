using RangeEnclosures
using IntervalArithmetic
using ForwardDiff
using IntervalOptimisation
#using MosekTools

g(x) = (x[1] + 2x[2] - 7)^2 + (2x[1] + x[2] - 5)^2;
f(x)=x[1]+x[2]x[3]+7
f2(x)=x[1]^2
f3(x)=f(x)*f2(x)
f3([1,1,1])
f2(x)=1/16*x[1]
f2(16)
f
dom = IntervalBox(-1..1, 3)
enclose(f,dom,BranchAndBoundEnclosure())

enclose(x -> 1+ 2*x+ x^2, -1 .. 1)


a(x)=1/16*x[1]^2+1/16*x[2]^2
b(x)=-1/16*(3*x[1]^2+3*x[2]^2 -x[1]^2*x[3]-x[2]^2*x[3]) 
c(x)=(9+6*x[3]+x[3]^2)*1/64*(x[1]^2+x[2]^2)

I_1(x) = a(x) - b(x) + c(x)
I_2(x)= -2*a(x) + b(x)
I_3(x)=b(x)^2/(4*a(x)) - c(x)
I_4(x)= a(x) + b(x) + c(x)
I_5(x) = 2*a(x) + b(x)

prod_1(x)=I_1(x)*I_2(x)
prod_2(x)=I_2(x)*I_3(x)
prod_3(x)=I_4(x)*I_5(x)
prod_4(x)=I_3(x)*I_5(x)

enclose(prod_1,dom,BranchAndBoundEnclosure())
enclose(prod_2,dom,BranchAndBoundEnclosure())
enclose(prod_3,dom,BranchAndBoundEnclosure())
enclose(prod_4,dom,BranchAndBoundEnclosure())
enclose(prod_3,dom,NaturalEnclosure())
enclose(prod_3,dom,MooreSkelboeEnclosure())
