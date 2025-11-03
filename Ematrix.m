function E = Ematrix

global N;

N=N+1;
syms x

NM1=N-1;
    % Chebyshev Gauss Lobatto nodes
    xc=cos(pi*(0:NM1)/NM1)';
    
    % Uniform nodes
    xu=linspace(-1,1,N)';

    % Make a close initial guess
    if NM1<3
        xx=xc;
    else
        xx=xc+sin(pi*xu)./(4*NM1);
    end

    P=zeros(N,N);% Use to compute the Legendre Vandermonde Matrix

    xold=2;
while max(abs(xx-xold))>eps

    xold=xx;
        
    P(:,1)=1;    P(:,2)=xx; %first two columns 1s and xs
    
    for k=2:NM1 %following columns follow recursion relation:
        
        P(:,k+1)=( (2*k-1)*xx.*P(:,k)-(k-1)*P(:,k-1) )/k;
    end
     
    %Update x using Newton-Raphson method
    xx=xold-( xx.*P(:,N)-P(:,NM1) )./( N*P(:,N) );
end

xx=flipud(xx);
X=repmat(xx,1,N); %Creates a square matrix with column entries x
   
    Xdiff=X-X'+eye(N); %eye is for Kronecker Delta
    
    L=repmat(P(:,N),1,N); %Replicates Legendre Vandemonde Matrix
    L(1:(N+1):N*N)=1;
    
    D=(L./(Xdiff.*L'));
    D(1:(N+1):N*N)=0;
    D(1)=-(N*NM1)/4;
    D(N*N)=(N*NM1)/4; %Gives differentiation matrix
 
   for l=1:N
        massmatrix(l,l)=1; 
   end
    
  E=massmatrix;
 
  for i=1:N
    E(1,i) = 0.5*(1-(xx(i))^2)*D(1,i);
  end
  for i=1:N
    E(N,i) = -0.5*(1-(xx(i))^2)*D(N,i);
  end
  
  
   E=E(:,2:N-1);
      
   N=N-1;
 
end
    