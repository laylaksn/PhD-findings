function wEE = MassMatrixtt

    global N;

    w=Weights;
    E=Ematrix;
    mm=1:N+1;  

    for i=1:N-1
            for j=1:N-1
                wEE(i,j) = sum(w(mm).*E(mm,i).*E(mm,j));

            end
    end

return