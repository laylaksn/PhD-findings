function At = StiffnessMatrixt

    global N;

    w=Weights;
    D=DMatrix;
    Dtilde=DtildeMatrix;
    
    n=1:N+1;
    for i=1:N+1
            for j=1:N-1
                At(i,j)=sum(w(n).*D(n,i).*Dtilde(n,j+1));
            end
    end

return