function Lt=MixMatrixt

    global N;
    
    w=Weights;
    Dtilde=DtildeMatrix;

    for i=1:N+1
        for j=1:N-1
            Lt(i,j)=(w(i).*Dtilde(i,j+1));
        end
    end

return