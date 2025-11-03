function wE = MassMatrixt

    global N;

    w=Weights;
    E=Ematrix;
 
    for i=1:N+1
        for j=1:N-1
            wE(i,j) = (w(i).*E(i,j));
        end
    end

return
