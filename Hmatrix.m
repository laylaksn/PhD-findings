function H = Hmatrix

    global N;
   
    w=Weights;
    E=Ematrix;

    for k=1:N-1
        for m=1:N+1
            H(k) = sum(w(m)*E(m,k));
        end
    end
    
end