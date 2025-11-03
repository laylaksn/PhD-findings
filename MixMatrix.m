function L=MixMatrix

    global N;

    w=Weights;
    D=DMatrix;

    for i=1:N+1
            for j=1:N+1
                  L(i,j)=(w(i).*D(i,j));
            end
    end

return
