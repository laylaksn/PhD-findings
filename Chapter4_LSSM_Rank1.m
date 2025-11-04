% =====================================================================
% The first section of this code (computation of nodes, differentiation
% matrix, and weights) is based on:
%
%   C. Canuto, M. Y. Hussaini, A. Quarteroni, T. A. Tang,
%   "Spectral Methods in Fluid Dynamics," Section 2.3,
%   Springer-Verlag, 1987.
%
% Implementation originally written by:
%   Greg von Winckel, 05/26/2004
%   Contact: gregvw@chtm.unm.edu
%
% ================================================================
% All subsequent code beyond this is written by Layla Sadeghi Namaghi
% ================================================================

clc 
clear
close

format short

% Number of enrichments is 1
no_enrichment=1;

% N GLL points
Nvalues=[8,12,16,20] ;

% Max number of iterations
max_iterations = 1000;

% Stopping criterion
epsilon=10^(-10);

cumTime=0;
executionTimes = [];
 
L2phi=zeros(length(Nvalues),1);
L2U=zeros(length(Nvalues),1);
L2V=zeros(length(Nvalues),1); 
for number=1:length(Nvalues)
    each_count_value=[];
    N=Nvalues(number);
    N1=N+1;

    NN=N;
 
    % Chebyshev Gauss Lobatto nodes
    xc=cos(pi*(0:NN)/NN)';
    
    % Uniform nodes
    xu=linspace(-1,1,N1)';

    % Make a close initial guess
    if NN<3
        xx=xc;
    else
        xx=xc+sin(pi*xu)./(4*NN);
    end

    P=zeros(N1,N1);% Use to compute the Legendre Vandermonde Matrix

    xold=2;
    while max(abs(xx-xold))>eps
    
        xold=xx;
            
        P(:,1)=1;    P(:,2)=xx; %first two columns 1s and xs
        
        for k=2:NN %following columns follow recursion relation:
            P(:,k+1)=( (2*k-1)*xx.*P(:,k)-(k-1)*P(:,k-1) )/k;
        end
         
        %Update x using Newton-Raphson method
        xx=xold-( xx.*P(:,N1)-P(:,NN) )./( N1*P(:,N1) );
    end

    xx=flipud(xx);
    node=xx;

%============= Differentiation matrix and weights ==================

    XX=repmat(xx,1,N1); %Creates a square matrix with column entries x
    Xdiff=XX-XX'+eye(N1); %eye is for Kronecker Delta
        
    L=repmat(P(:,N1),1,N1); %Replicates Legendre Vandemonde Matrix
    L(1:(N1+1):N1*N1)=1;
    D=(L./(Xdiff.*L'));
    D(1:(N1+1):N1*N1)=0;
    D(1)=-(N1*NN)/4;
    D(N1*N1)=(N1*NN)/4; %Gives differentiation matrix
    
    w=2./(NN*N1*(L(:,N1).^2)); %Computation of weights

%=======================================================================
% Set up matrices
%=======================================================================   

    A = zeros(N1); M = zeros(N1); L = zeros(N1);

    for i=1:N1
        for j=1:N1
            sum=0;
            for n=1:N1
                sum=sum+(w(n)*D(n,i)*D(n,j));
            end
            A(i,j)=sum;
        end
    end
    for i=1:N1
        M(i,i)=w(i);
    end
    for i=1:N1
        for j=1:N1
            L(i,j)=(w(i).*D(i,j));
        end
    end

    % source term
    F = (2)*(pi^2)*sin(pi*node).*(sin(pi*node'));

    % true solutions
    syms sx sy
    phi=@(sx,sy) sin(pi*sx)*sin(pi*sy); 
    U=@(x,y) - pi*cos(pi*x).*sin(pi*y); 
    V=@(x,y) - pi*sin(pi*x).*cos(pi*y);

    B =zeros(N1);
    for i=1:N1
        for j=1:N1
            sum=0;
            for n=1:N1
                sum=sum+(w(n)*w(j)*D(n,i)*F(n,j));
            end
            B(i,j)=sum;
        end
    end

    C=zeros(N1);
    for i=1:N1
        for j=1:N1
            sum=0;
            for n=1:N1
                sum=sum+w(i)*w(n)*F(i,n)*D(n,j);
            end
            C(i,j)=sum;
        end
    end

    % Introducing (N−1)×(N−1) matrices
    A2 = A(2:N, 2:N);
    M2 = M(2:N, 2:N);
    L2 = L(2:N, 2:N);
    B2 = B(2:N, 2:N);
    C2 = C(2:N, 2:N);
    
    % Introducing (N−1)×(N+1) matrices
    A3 = A(2:N, :);
    M3 = M(2:N, :);
    L3 = L(2:N, :);
    B3 = B(2:N, :);
    C3 = C(2:N, :);
    
    % Introducing (N+1)×(N−1) matrices
    A4 = A(:, 2:N);
    M4 = M(:, 2:N);
    L4 = L(:, 2:N);
    B4 = B(:, 2:N);
    C4 = C(:, 2:N);

%=======================================================================
%       PGD begins
%=======================================================================
    %Start off Q
    Q=0;

    %Start approximations
    Xphi = zeros(N+1,Q ) ;    Xu = zeros(N+1,Q ) ;    Xv = zeros(N+1,Q ) ;
    Yphi = zeros(N+1,Q ) ;    Yu = zeros(N+1,Q ) ;    Yv = zeros(N+1,Q ) ;

    %Qth enrichments
    for term = 1:no_enrichment
        
        startTime(number) = cputime;
        
       
        %Unknown functions for enrichment
        R1 = ones( 3*N-1 ,1) ;
        S1 = ones( 3*N-1 ,1) ;
        oldR = zeros( 3*N-1 ,1) ;
        oldS = zeros( 3*N-1 ,1) ;
        Rphi = rand(N-1 ,1) ;        Ru = rand(N+1, 1 ) ;        Rv = rand(N-1, 1 ) ;
        Sphi = rand(N-1, 1 ) ;       Su = rand(N-1 ,1) ;         Sv = rand(N+1 ,1) ;
        counter = 0 ;
        
        % Iteration loop
        while max(max( abs (R1*S1' - oldR* oldS') ) ) >epsilon
            oldR = R1 ;
            oldS = S1 ;

            Mr11 = (Rphi'*A2*Rphi)*M2 + (Rphi'*M2*Rphi)*A2;
            Mr12 = (Rphi'*L4'*Ru)*M2 ;
            Mr13 = (Rphi'*M2*Rv)*L4';
            Mr21 = (Ru'*L4*Rphi)*M2;
            Mr22 = (Ru'*A*Ru) *M2 + (Ru'*M*Ru) *A2 + (Ru'*M*Ru) *M2 ;
            Mr23 = (Ru'*L3'*Rv)*L3 - (Ru'*L4*Rv)*L4';
            Mr31 = (Rv'*M2*Rphi)*L4;
            Mr32 = (Rv'*L3*Ru)*L3' - (Rv'*L4'*Ru)*L4 ;
            Mr33 = (Rv'*A2*Rv)*M + (Rv'*M2*Rv)*A + (Rv'*M2*Rv)*M;

            Matrix2 = [Mr11, Mr12, Mr13; Mr21, Mr22, Mr23; Mr31, Mr32, Mr33] ;
            Vr1 = zeros(N-1,1);
            Vr2 = B4'*Ru;
            Vr3 = C3'*Rv;
            Vector2 = [ Vr1 ; Vr2 ; Vr3 ] ;

            S1 = Matrix2\Vector2; 
            Sphi = S1(1:N-1);
            Su = S1(N:2*N-2);
            Sv = S1(2*N-1:3*N-1);

            Ms11 = (Sphi'*M2*Sphi)*A2+(Sphi'*A2*Sphi)*M2;
            Ms12 = (Sphi'*M2*Su)*L4';
            Ms13 = (Sphi'*L4'*Sv)*M2;
            Ms21 = (Su'*M2*Sphi)*L4;
            Ms22 = (Su'*M2*Su)*A+(Su'*A2*Su)*M+(Su'*M2*Su)*M;
            Ms23 = (Su'*L3*Sv)*L3'-(Su'*L4'*Sv)*L4;
            Ms31 = (Sv'*L4*Sphi)*M2;
            Ms32 = (Sv'*L3'*Su)*L3-(Sv'*L4*Su)*L4';
            Ms33 = (Sv'*A*Sv)*M2+(Sv'*M*Sv)*M2+(Sv'*M*Sv)*A2;
            Matrix1 = [Ms11,Ms12,Ms13;Ms21,Ms22,Ms23;Ms31,Ms32,Ms33];

            Vs1 = zeros(N-1,1) ;
            Vs2 = B4*Su; 
            Vs3 = C3*Sv;
            Vector1 = [ Vs1 ; Vs2 ; Vs3 ] ;

            R1 = Matrix1\Vector1 ;
            Rphi = R1(1:N-1);
            Ru = R1(N:2*N);
            Rv = R1(2*N+1:3*N-1);
          
            counter=counter+1;
            if counter == max_iterations
                break
            end
    
        end
        endTime(number) = cputime;     
        each_count_value = [each_count_value, counter];

        % Calculate execution time for this iteration
        enrichmentTime(number) = endTime(number) - startTime(number);

        % cumTime = cumTime + enrichmentTime;
       
        %setting XQ YQ from the calculated values
        xphi = Rphi ;        xu= Ru;        xv = Rv ;
        yphi = Sphi ;        yu = Su ;      yv = Sv ;

        %adding the known boundary values
        Xphi = [Xphi, [0; xphi; 0]];
        Xu   = [Xu, xu];
        Xv   = [Xv, [0; xv; 0]];
        
        Yphi = [Yphi, [0; yphi; 0]];
        Yu   = [Yu, [0; yu; 0]];
        Yv   = [Yv, yv];
        Q= Q+1;
        
        %calculating the approximations
        approxphi = zeros(N+1) ;
        approxu = zeros(N+1) ;
        approxv = zeros(N+1) ;
        for k = 1 :N+1
            for l = 1 :N+1
                for j = 1 : Q
                    approxphi(l,k) = approxphi(l,k) + Yphi(k,j)*Xphi(l,j);
                    approxu(l,k) = approxu(l,k) + Yu(k,j)*Xu(l,j) ;
                    approxv(l,k) = approxv(l,k) + Yv(k,j)*Xv(l,j) ;
                end
            end
        end
            
        %working out the error
        PHIerror= zeros(N+1) ;            PHImesh= zeros(N+1) ;
        Umesh = zeros(N+1) ;              Vmesh = zeros(N+1) ;
        Uerror = zeros(N+1) ;             Verror = zeros(N+1) ;

        for k = 1 :N+1
            for l = 1 :N+1
                Uerror(k,l) = approxu(k,l) - U(node(k),node(k));
                Verror(k,l) = approxv(k,l) - V(node(k),node(l));
                PHIerror(k,l) = approxphi(k,l) - phi(node(k),node(l));
                PHImesh(k,l) = phi(node(k),node(l));
                Umesh(k,l) = U(node(k),node(l));
                Vmesh(k,l) = V(node(k),node(l));
            end
        end
    
        QL2phierror=0;
        QL2uerror=0;
        QL2verror=0;
        
        for num=1:N1
            for m=1:N1
                QL2phierror = QL2phierror + ((w(m)*w(num)*(approxphi(m,num)-PHImesh(m,num))^2)) ;
                QL2uerror = QL2uerror + ((w(m)*w(num)*(approxu(m,num)-Umesh(m,num))^2)) ;
                QL2verror = QL2verror + ((w(m)*w(num)*(approxv(m,num)-Vmesh(m,num))^2)) ;
            end
        end
        QL2phierror = QL2phierror^(1/2);
        QL2uerror = QL2uerror^(1/2);
        QL2verror = QL2verror^(1/2);
        
        QL2phi(Q,number)=QL2phierror;
 
    end


    each_count_b(number,:) = each_count_value;
    L2phi(number)=QL2phierror;
    L2U(number)=QL2uerror;
    L2V(number)=QL2verror;
end

NP1 = N+1; 
outersum1=0; 
for m = 1:NP1
innersum1=0;    
   for n =1:NP1
       sum1=0;
       sum2=0; 
            for k=1:NP1   
               sum1 = sum1+(D(m,k)*approxu(k,n));
               sum2 = sum2+(D(n,k)*approxv(m,k));  
            end
            dvdy(m,n) = sum1;
            dudx(m,n) = sum2; 
    innersum1 = innersum1 + w(n)*(dudx(m,n)+dvdy(m,n)-F(m,n))^2; 
    end
    outersum1 = outersum1 + innersum1*w(m); 
end
J1 = outersum1 ;
 

outersum1=0;
outersum2=0;
for m = 1:NP1
innersum1=0; 
innersum2=0;
    for n =1:NP1
        sum1=0;
        sum2=0; 
            for k=1:NP1
               sum1 = sum1+(D(m,k)*approxphi(k,n));
               sum2 = sum2+(D(n,k)*approxphi(m,k)); 
            end
               dphidx(m,n) = sum2;
               dphidy(m,n) = sum1;  
    innersum1 = innersum1 + w(n)*(dphidy(m,n) + approxu(m,n) )^2; 
    innersum2 = innersum2 + w(n)*(dphidx(m,n) + approxv(m,n) )^2; 
    end
    outersum1 = outersum1 + innersum1*w(m); 
    outersum2 = outersum2 + innersum2*w(m); 
end
J2 = outersum1;
J3 = outersum2 ;

outersum1=0;
for m = 1:NP1
innersum1=0;    
    for n =1:NP1
        sum1=0;
        sum2=0; 
            for k=1:NP1
               sum1 = sum1+(D(m,k)*approxu(k,n));
               sum2 = sum2+(D(n,k)*approxv(m,k)); 
            end
               dudx(m,n) = sum1;
               dvdy(m,n) = sum2; 
            
    innersum1 = innersum1 + w(n)*(-dudx(m,n)+dvdy(m,n))^2; 
    end
    outersum1 = outersum1 + innersum1*w(m); 
end
J4 = outersum1 ;
 

J1=J1^2;
J2=J2^2;
J3=J3^2;
J4=J4^2;
funcval = [J1,J2,J3,J4];

N_values=Nvalues'; 
 
fprintf('L2 errors for Example 4, i=1, using the LS SM PGD, after one enrichment:')
fprintf('\n');
fprintf('\n'); 
fprintf('N =       ');
fprintf('%8d ', Nvalues);
fprintf('\n-----------------------------------------------------\n');
fprintf('L2 error p   ');
fprintf('%6.2e ', L2phi); 
fprintf('\n');
fprintf('\n');

fprintf('L2 errors for u when N=20 : ')
disp(L2U(end))
fprintf('L2 errors for v when N=20 : ')
disp(L2V(end))

each_count_bT = each_count_b';
fprintf('\n');
fprintf('Number of ADFPA iterations required at the first enrichment stage for Example 9 using the LS SM PGD, with epsilon =  %6.2e', epsilon)
fprintf('\n');
fprintf('\n'); 
fprintf('N =          ');
fprintf('%8d ', Nvalues);
fprintf('\n-----------------------------------------------------\n');
fprintf('ADFPA iter   ');
fprintf('%8d ', each_count_bT); 
fprintf('\n');
fprintf('\n');



fprintf('Time taken (in seconds):')
fprintf('\n');
fprintf('N =          ');
fprintf('%8d ', Nvalues);
fprintf('\n');
fprintf('Time       ');
disp( enrichmentTime); 

fprintf('J1, J2, J3, J4 evaluated respectively (these are presented in Chapter 7 of the thesis):')
fprintf('   ');
fprintf('%6.2e ', funcval); 

