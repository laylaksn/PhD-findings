function w=Weights
global N

%   C. Canuto, M. Y. Hussaini, A. Quarteroni, T. A. Tang,
%   "Spectral Methods in Fluid Dynamics," Section 2.3,
%   Springer-Verlag, 1987.
%
% Implementation originally written by:
%   Greg von Winckel, 05/26/2004
%   Contact: gregvw@chtm.unm.edu

%setting up nodes
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

    w=2./(N*NP1*(LL(:,NP1).^2)); 
   end