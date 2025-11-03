function M =MassMatrix
 
    global N;

    w=Weights;
   
    for i=1:N+1
        M(i,i)=w(i);
    end
    
    
end
 