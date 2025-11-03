function Att = StiffnessMatrixtt

    global N;
    
    w=Weights;
    Dtilde=DtildeMatrix;
    
    n=1:N+1;
    for i=1:N-1
        for j=1:N-1
            Att(i,j)=sum(w(n).*Dtilde(n,i+1).*Dtilde(n,j+1));
        end
    end

return
