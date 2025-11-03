function A = StiffnessMatrix

    global N;

    w=Weights;
    D=DMatrix;
    n=1:N+1;
   
    for i=1:N+1
        for j=1:N+1
        A(i,j)=sum(w(n).*D(n,i).*D(n,j));
        end
    end
  
return