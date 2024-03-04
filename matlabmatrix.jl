function MatlabMatrix(SSPZ::SimpleSparsePolynomialZonotope,file,name)
    c=SSPZ.c 
    E=SSPZ.E 
    G=SSPZ.G 
    open(file, "w") do f
        write(f, "centre"*name*" = [")
        for i in 1:length(c)
            if i!=length(c)
                write(f,"$(c[i]) ;")
            else
                write(f,"$(c[i])] \n\n")
            end
        end
        write(f," \n")
        write(f, "generateurs"*name*" = [")
        for i in 1:size(G)[1]
            for j in 1:size(G)[2]
                write(f,"$(G[i,j]) ")
            end
            if i!=size(G)[1]
                write(f,";")
            else
                write(f,"] \n\n")
            end
        end

        write(f, "exposants"*name*" = [")
        println(size(E)[1])
        for i in 1:size(E)[1]
            
            for j in 1:size(E)[2]
                write(f,"$(E[i,j]) ")
            end
            if i!=size(E)[1]
                write(f,";")
            else
                write(f,"] \n")
            end
        end

    end
end




