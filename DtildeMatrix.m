function Dtilde=DtildeMatrix
global N
 
    %   C. Canuto, M. Y. Hussaini, A. Quarteroni, T. A. Tang,
    %   "Spectral Methods in Fluid Dynamics," Section 2.3,
    %   Springer-Verlag, 1987.
    %
    % Implementation originally written by:
    %   Greg von Winckel, 05/26/2004
    %   Contact: gregvw@chtm.unm.edu
    
    NP1=N+1;
       % Chebyshev Gauss Lobatto nodes
    xc=cos(pi*(0:N)/N)';
    
    % Uniform nodes
    xuu=linspace(0,1,NP1)';

    % Make a close initial guess
    if N<3
        xx=xc;
    else
        xx=xc+sin(pi*xuu)./(4*N);
    end

    P=zeros(NP1,NP1);% Use to compute the Legendre Vandermonde Matrix

    xold=2;
while max(abs(xx-xold))>eps

    xold=xx;
        
    P(:,1)=1;    P(:,2)=xx; %first two columns 1s and xs
    
    for k=2:N %following columns follow recursion relation:
        P(:,k+1)=( (2*k-1)*xx.*P(:,k)-(k-1)*P(:,k-1) )/k;
    end
     
    %Update x using Newton-Raphson method
    xx=xold-( xx.*P(:,NP1)-P(:,N) )./( NP1*P(:,NP1) );
end
 xx=flipud(xx);

X=repmat(xx,1,NP1); %Creates a square matrix with column entries x
    Xdiff=X-X'+eye(NP1); %eye is for Kronecker Delta
    
    LL=repmat(P(:,NP1),1,NP1); %Replicates Legendre Vandemonde Matrix
    LL(1:(NP1+1):NP1*NP1)=1;
    D=(LL./(Xdiff.*LL'));
    D(1:(NP1+1):NP1*NP1)=0;
    D(1)=-(NP1*N)/4;
    D(NP1*NP1)=(NP1*N)/4; %Gives differentiation matrix

    %========= Pressure Differentiation Matrix ===============================
    % Written by Layla Sadeghi Namaghi

    %All entries for j not 0 or N
   for i=2:N
       for q=2:N
            Dtilde(q,i)=((1-xx(i)^2)/(1-xx(q)^2))*D(q,i);
       end
   end

   %Diagonal entries i=j
   for i=2:N
      Dtilde(i,i)=D(i,i)+2*xx(i)/(1-xx(i)^2) ;
   end
  Dsq=zeros(N+1);

  %second differentiation matrix
ii=1:N+1;
for q=1:N+1
    for k=1:N+1
        Dsq(q,k)=sum(D(q,ii)*D(ii,k));
    end
end

% Dsq
   %When j=0 and j=N
   for i=1:N+1
     Dtilde(1,i)=(1-(xx(i))^2)*(D(1,i)+(Dsq(1,i)))/4;
     Dtilde(N+1,i)=(1-(xx(i))^2)*(D(N+1,i)-(Dsq(N+1,i)))/4;
   end
  
end