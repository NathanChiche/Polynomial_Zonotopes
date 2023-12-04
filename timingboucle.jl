using Dates

function boucleparligne(A)
    n,m=size(A)
    sum=0
    for i in 1:n
        for j in 1:m
            sum+=A[i,j]
        end
    end
    return sum
end

function boucleparcolonne(A)
    n,m=size(A)
    sum=0
    for j in 1:m
        for i in 1:n
            sum+=A[i,j]
        end
    end
    return sum
end

@time M=zeros(Float64,10,10000000)
@time boucleparligne(M)
@time boucleparcolonne(M)

start_time = now()
fin=iterate_polynomials_over_PZ([chatal1,chatal2],PChatal,5,0,R,"bary",max_order=1000000000)
end_time = now()
elapsed = end_time - start_time
println("temps des iterations:", elapsed)