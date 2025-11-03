function wED = Gmatrix

    global N;
    syms x
    
    w=Weights;
    D=DMatrix;
    E=Ematrix;
    mm=1:N+1;
    for i=1:N-1
            for j=1:N+1
                wED(i,j) = sum(w(mm).*E(mm,i).*D(mm,j));
            end
    end

return
