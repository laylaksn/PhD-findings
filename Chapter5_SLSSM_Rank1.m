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

clear all
clc
close all

syms x y
format short

% true solution
phi=@(x,y) sin(pi*x)*sin(pi*y); 
U=@(x,y) - pi*cos(pi*x).*sin(pi*y)  ;
V=@(x,y) - pi*sin(pi*x).*cos(pi*y)  ;
     
% b values
b_values = [0, 1;
            1, 0;
            0.01, 1;
            1, 0.01;
            0.1, 1;
            1, 0.1]; 

% epsilon bar values
eps_values = [1; 0.001 ]; 

% Number of enrichments
no_enrichment=1;

% N GLL points
Nvalues=[8,12,16,20,24] ;

% to turn off weak BS (0), or on (1)
weakBC =  1; 
 
% Max number of iterations
max_iterations = 100;

% Tolerance
epsilon=10^(-10);
  
L2phi=zeros(length(Nvalues),1);
L2U=zeros(length(Nvalues),1);
L2V=zeros(length(Nvalues),1); 
 
for ei = 1:size(eps_values,1)
    ep = eps_values(ei);

    for bi = 1:size(b_values,1)

        b = [b_values(bi,1) b_values(bi,2)];
        for number=1:length(Nvalues) 
            N=Nvalues(number);
           
            NN=N;
            N1=N+1;
               % Chebyshev Gauss Lobatto nodes
            xc=cos(pi*(0:NN)/NN)';

            % Uniform nodes
            xu=linspace(0,1,N1)';
        
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
            
            xs=node;
            ys=node;
         
            phi=@(x,y) sin(pi*x)*sin(pi*y) ; 
            U=@(x,y)- pi*cos(pi*x).*sin(pi*y)  ;
            V=@(x,y)- pi*sin(pi*x).*cos(pi*y)  ;
        
            F = zeros(N+1, N+1);  
            for i = 1:N+1
                for j = 1:N+1
                    F(i,j) = ep*((2)*(pi^2)*sin(pi*xs(i))*(sin(pi*ys(j))) ); 
                    F(i,j) = F(i,j) + b(1)*(pi*cos(pi*xs(i))*sin(pi*ys(j)) ); 
                    F(i,j) = F(i,j) + b(2)*(pi*sin(pi*xs(i))*cos(pi*ys(j)) ); 
                end
            end
     
       
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
            B=zeros(N1);
            for i=1:N1
                for j=1:N1
                    sum=0;
                    for n=1:N1
                    sum=sum+(w(n)*w(j)*D(n,i)*F(n,j));
                    end
                    B(i,j)=sum;
                end
            end
            
            DD=zeros(N1);
            for i=1:N1
                for j=1:N1
                    sum=0;
                    for n=1:N1
                     sum=sum+w(i)*w(n)*F(i,n)*D(n,j);
                    end
                    DD(i,j)=sum;
                end
            end
     
            %Introducing (N−1)ˆ2 matrices
            A2 = A( 2 :N, 2 :N) ;
            M2 = M( 2 :N, 2 :N) ;
            L2 = L ( 2 :N, 2 :N) ;
            B2 = B ( 2 :N, 2 :N) ;
            D2 = DD( 2 :N, 2 :N) ;
        
            %Introducing (N−1) *(N+1) sized matrices
            M3 = M( 2 :N, : ) ;
            A3 = A( 2 :N, : ) ;
            L3 = L ( 2 :N, : ) ;
            B3 = B ( 2 :N, : ) ;
            D3 = DD( 2 :N, : ) ;
        
            %Introducing (N+1) *(N−1) sized matrix
            L4 = L ( : , 2 :N) ;
            M4 = M( : , 2 :N) ;
            A4 = A ( : , 2 :N) ;
            B4 = B ( : , 2 :N) ;
            D4 = DD ( : , 2 :N) ;
            %Start off Q
            Q=0;
        
            %Start approximations
            Xphi = zeros(N+1,Q ) ;
            Xu = zeros(N+1,Q ) ;
            Xv = zeros(N+1,Q ) ;
            Yphi = zeros(N+1,Q ) ;
            Yu = zeros(N+1,Q ) ;
            Yv = zeros(N+1,Q ) ;

            alpha_1 = zeros(N+1,1);
            alpha_2 = zeros(N+1,1);
            alpha_3 = zeros(N+1,1);
            alpha_4 = zeros(N+1,1);
        
            alpha_1(1) = ep + max(b(2)/sqrt(2) - b(1)/sqrt(2),0);
            for n = 2:N+1
                alpha_1(n) = ep;
            end
        
            alpha_2(1) = ep + max(b(1)/sqrt(2) - b(2)/sqrt(2),0);
            for n = 2:N+1
                alpha_2(n) = ep;
            end
        
            alpha_3(1) = ep + max(b(2)/sqrt(2) + b(1)/sqrt(2),0);
            for n = 2:N
                alpha_3(n) = ep + b(1);
            end
            alpha_3(N+1) = ep + max(b(1)/sqrt(2) - b(2)/sqrt(2),0);
        
            alpha_4(1) = ep + max(b(2)/sqrt(2) + b(1)/sqrt(2),0);
            for n = 2:N
                alpha_4(n) = ep + b(2);
            end
            alpha_4(N+1) = ep + max(b(2)/sqrt(2) - b(1)/sqrt(2),0);

            %enrichments begin
            for term = 1:no_enrichment

                startTime = cputime;
        
                Q=Q;
        
                %Initialize the Matrices containing the known Q−1 enrichments
                Pphi = zeros(N+1) ;
                Pu = zeros(N+1,N+1) ;
                Pv = zeros(N+1,N+1) ;
        
                XphiA = zeros(N+1,Q ) ;
                XphiM= zeros(N+1,Q ) ;
                XphiL= zeros(N+1,Q ) ;
                XphiLT =zeros(N+1,Q);
                XuA= zeros(N+1,Q ) ;
                XuM= zeros(N+1,Q ) ;
                XuL= zeros(N+1,Q ) ;
                XuLT = zeros(N+1,Q ) ;
                XvA= zeros(N+1,Q ) ;
                XvM= zeros(N+1,Q ) ;
                XvL= zeros(N+1,Q ) ;
                XvLT = zeros(N+1,Q ) ;
        
                YphiA= zeros(N+1,Q ) ;
                YphiM= zeros(N+1,Q ) ;
                YphiL= zeros(N+1,Q ) ;
                YphiLT=zeros(N+1, Q);
                YuA= zeros(N+1,Q ) ;
                YuM= zeros(N+1,Q ) ;
                YuL= zeros(N+1,Q ) ;
                YuLT = zeros(N+1,Q ) ;
                YvA= zeros(N+1,Q ) ;
                YvL= zeros(N+1,Q ) ;
                YvLT = zeros(N+1,Q ) ;
                YvM= zeros(N+1,Q ) ;

                for i = 1 :N+1
                    for j = 1 : Q
                         for k = 1 :N+1
                            XphiM( i , j ) = XphiM ( i , j ) + Xphi ( k , j ) *M( k , i ) ;
                            XuLT( i , j ) = XuLT( i , j ) + Xu( k , j ) *L( k , i);
                            YphiA ( i , j ) = YphiA ( i , j ) + Yphi ( k , j ) *A( k , i ) ;
                            YphiM( i , j ) = YphiM ( i , j ) + Yphi ( k , j ) *M( k , i ) ;
                            YuM( i , j ) = YuM( i , j ) + Yu( k , j ) *M( k , i ) ;
                            YvLT( i , j ) = YvLT( i , j ) + Yv( k , j ) *L( k , i ) ;
                             XuL( i , j ) = XuL( i , j ) + Xu( k , j ) *L ( i , k ) ;
                             XvA( i , j ) = XvA( i , j ) + Xv( k , j ) *A( k , i ) ;
                             YuA( i , j ) = YuA( i , j ) + Yu( k , j ) *A( k , i ) ;
                            XuA( i , j ) = XuA( i , j ) + Xu( k , j ) *A( k , i ) ;
                            XuM( i , j ) = XuM( i , j ) + Xu( k , j ) *M( k , i ) ;
                            XvL( i , j ) = XvL( i , j ) + Xv( k , j ) *L( i , k ) ;
                            XvLT( i , j ) = XvLT( i , j ) +Xv( k , j ) *L( k , i ) ;
                            YphiL ( i , j ) = YphiL ( i , j ) + Yphi ( k , j ) *L( i , k ) ;
                            YuL( i , j ) = YuL( i , j ) + Yu( k , j ) *L( i , k ) ;
                            YuLT( i , j ) = YuLT( i , j ) + Yu( k , j ) *L( k , i ) ;
                            YvA( i , j ) = YvA( i , j ) + Yv( k , j ) *A( k , i ) ;
                            YvM( i , j ) = YvM( i , j ) + Yv( k , j ) *M( k , i ) ;
                            XphiA ( i , j ) = XphiA ( i , j ) + Xphi ( k , j ) *A( i , k ) ;
                            XphiLT ( i , j ) = XphiLT ( i , j ) + Xphi ( k , j ) *L( k , i ) ;
                            YphiLT ( i , j ) = YphiLT ( i , j ) + Yphi ( k , j ) *L( k , i ) ;
                            XvM( i , j ) = XvM( i , j ) + Xv( k , j ) *M( k , i ) ;
                            XphiL ( i , j ) = XphiL ( i , j ) + Xphi ( k , j ) *L ( i , k ) ;
                            YvL ( i , j ) = YvL ( i , j ) + Yv( k , j ) *L( i , k ) ;
                         end 
                    end
                end
        
                for i = 1 :N+1
                    for j = 1 : Q
                        for k = 1 :N+1                                                                                                             %here                        %here                                                                                                   %here                        %here                                                                                                                                                
                            Pphi(k,i) = Pphi(k,i) + (XuLT(k,j)*YuM(i,j) + XphiA(k,j)*YphiM(i,j) + XvM(k,j)*YvLT(i,j) + XphiM(k,j)*YphiA(i,j) + ep*b(1)*XuA(k,j)*YuM(i,j) + ep*b(1)*XvLT(k,j)*YvL(i,j) + b(1)*b(1)*XphiA(k,j)*YphiM(i,j) + b(1)*b(2)*XphiLT(k,j)*YphiL(i,j) + ep*b(2)*XuL(k,j)*YuLT(i,j) + ep*b(2)*XvM(k,j)*YvA(i,j) + b(2)*b(1)*XphiL(k,j)*YphiLT(i,j) + b(2)*b(2)*XphiM(k,j)*YphiA(i,j));
                            Pu(k,i) = Pu(k,i) + (ep*b(1)*XphiA(k,j)*YphiM(i,j)  + ep*b(2)*XphiLT(k,j)*YphiL(i,j) + ep*ep*XuA(k,j)*YuM(i,j) + ep*ep*XvLT(k,j)*YvL(i,j) + XuM(k,j)*YuM(i,j) + XphiL(k,j)*YphiM(i,j) - XvL(k,j)*YvLT(i,j) + XuM(k,j)*YuA(i,j));
                            Pv(k,i) = Pv(k,i) + (ep*b(1)*XphiL(k,j)*YphiLT(i,j) + ep*b(2)*XphiM(k,j)*YphiA(i,j)  + ep*ep*XvM(k,j)*YvA(i,j) + ep*ep*XuL(k,j)*YuLT(i,j) + XvM(k,j)*YvM(i,j) + XphiM(k,j)*YphiL(i,j) - XuLT(k,j)*YuL(i,j) + XvA(k,j)*YvM(i,j));
                        end
                    end
                end


                Pu = Pu(:,2:end-1);
                Pv = Pv(2:end-1,:);
        
                %Initial guesses for first iteration
                R1 = ones( 3*N+1 ,1) ;
                S1 = ones( 3*N+1 ,1) ;
                oldR = zeros( 3*N+1 ,1) ;
                oldS = zeros( 3*N+1 ,1) ;
                Rphi = rand(N+1 ,1) ;
                Ru = rand(N+1, 1 ) ;
                Rv = rand(N-1, 1 ) ;

                if term >1 
                    Rphi = Xphi(:,term-1);
                    Ru = Xu(:,term-1);
                    Rv = Xv(2:end-1,term-1);
                else
                end
        
                counter = 0 ; 

                %iterations begin
                while max(max( abs (R1*S1' - oldR* oldS') ) ) >epsilon

                    oldR = R1 ;
                    oldS = S1 ;
    
                    %LHS computation for first system
                    % contributions from dI/ds_0^phi
                    djdS_0a = alpha_1(1)*(Rphi(N+1))^2*w(1);
                    djdS_0b = alpha_3(1)*(Rphi(1))^2*w(1);
    
                    djdS_0c = 0;
                    for n=1:N+1
                        djdS_0c = djdS_0c + alpha_4(1)*w(n)*(Rphi(n))^2; 
                    end
    
                    djdSphi_0 =  djdS_0a + djdS_0b + djdS_0c;
                    djdSphi_0 = weakBC*djdSphi_0;

                    % dI/dSphi_k, k=1,...,N-1
                    for k =1:N-1
                        djdSphi_k(k) = alpha_1(k+1)*(Rphi(N+1))^2*w(k+1) + alpha_3(k+1)*(Rphi(1))^2*w(k+1); %+-
                    end
                    djdSphi_k = weakBC*djdSphi_k;
    
                    %from dI/ds_N^phi
                    djdS_Na = alpha_1(N+1)*(Rphi(N+1))^2*w(N+1);
                    djdS_Nb = alpha_3(N+1)*(Rphi(1))^2*w(N+1);
    
                    djdS_Nc = 0;
                    for n=1:N+1
                        djdS_Nc = djdS_Nc + alpha_2(N+1)*w(n)*(Rphi(n))^2;
                    end
    
                    djdSphi_N =  djdS_Na + djdS_Nb + djdS_Nc;
                    djdSphi_N = weakBC*djdSphi_N;
    
                    Mr11 = (Rphi'*A*Rphi)*M + (Rphi'*M*Rphi)*A + b(1)*b(1)*(Rphi'*A*Rphi)*M + b(1)*b(2)*(Rphi'*L'*Rphi)*L+b(2)*b(1)*(Rphi'*L*Rphi)*L'+ b(2)*b(2)*(Rphi'*M*Rphi)*A ;
     
                    %add weak BCs down the diagonal
                    Mr11(1,1) = Mr11(1,1) + djdSphi_0;
                    Mr11(N+1,N+1) = Mr11(N+1,N+1) + djdSphi_N;

                    for d=1:N-1
                        Mr11(d+1,d+1) = Mr11(d+1,d+1) + djdSphi_k(d);
                    end
    
                    Mr12 = (Rphi'*L'*Ru)*M4       + b(1)*(Rphi'*A*Ru)*M4*ep  + b(2)*(Rphi'*L*Ru)*L3'*ep;
                    Mr13 = (Rphi'*M4*Rv)*L'       + b(1)*(Rphi'*L3'*Rv)*L*ep + b(2)*( Rphi'*M4*Rv)*A*ep;
                    Mr21 = (Ru'*L*Rphi)*M3        + b(1)*(Ru'*A*Rphi)*M3*ep  + b(2)*(Ru'*L'*Rphi)*L3*ep ; 
                    Mr22 = (Ru'*A*Ru)*M2*ep*ep    + (Ru'*M*Ru)*A2            + (Ru'*M*Ru)*M2 ; 
                    Mr23 = (Ru'*L3'*Rv)*L3*ep*ep  - (Ru'*L4*Rv)*L4';
                    Mr31 = (Rv'*M3*Rphi)*L        + b(1)*(Rv'*L3*Rphi)*L'*ep + b(2)*(Rv'*M3*Rphi)*A*ep;
                    Mr32 = (Rv'*L3*Ru)*L3'*ep*ep  - (Rv'*L4'*Ru)*L4 ;
                    Mr33 = (Rv'*A2*Rv)*M          + (Rv'*M2*Rv)*A*ep*ep         + (Rv'*M2*Rv)*M;  
                    Matrix2 = [Mr11, Mr12, Mr13; Mr21, Mr22, Mr23; Mr31, Mr32, Mr33] ;
       
                    %weak bcs for RHS enrichments >1
                    EVr_0a = 0;
                        for q=1:Q
                            EVr_0a = EVr_0a + Xphi(N+1,q)*Yphi(1,q);
                        end
                    Vr_0a = alpha_1(1)*EVr_0a*w(1)*Rphi(N+1);
    
                    EVr_0b = 0;
                        for q=1:Q
                            EVr_0b = EVr_0b + Xphi(1,q)*Yphi(1,q);
                        end
                    Vr_0b = alpha_3(1)*EVr_0b*w(1)*Rphi(1); 

                    Vr_0c = 0;
                    for n = 1:N+1
                        EVr_0c = 0;
                        for q=1:Q
                            EVr_0c = EVr_0c + Xphi(n,q)*Yphi(1,q);
                        end
                        Vr_0c = Vr_0c + alpha_4(1)*w(n)*Rphi(n)*EVr_0c; 
                    end
    
                    Vr_0 = Vr_0a + Vr_0b + Vr_0c;
                    Vr_0 = weakBC*Vr_0;
    
                    for k =1:N-1      
                        EVr_ka = 0;
                        for q=1:Q
                            EVr_ka = EVr_ka + Xphi(N+1,q)*Yphi(k+1,q);
                        end
                        Vr_ka(k) = EVr_ka*w(k+1)*alpha_1(k+1); 
                    end
                    Vr_ka = Vr_ka*Rphi(N+1);
    
                    for k =1:N-1    
                        EVr_kb = 0;
                        for q=1:Q
                            EVr_kb = EVr_kb + Xphi(1,q)*Yphi(k+1,q);
                        end
                        Vr_kb(k) = EVr_kb*w(k+1)*alpha_3(k+1);
                    end
                    Vr_kb = Vr_kb*Rphi(1);

                    Vr_k = Vr_ka + Vr_kb;
                    Vr_k = weakBC*Vr_k;
    
                    EVr_Na = 0;
                        for q=1:Q
                            EVr_Na = EVr_Na + Xphi(N+1,q)*Yphi(N+1,q);
                        end
                    Vr_Na = alpha_1(N+1)*w(N+1)*Rphi(N+1)*EVr_Na; 
    
                    EVr_Nb = 0;
                        for q=1:Q
                            EVr_Nb = EVr_Nb + Xphi(1,q)*Yphi(N+1,q);
                        end
                    Vr_Nb = alpha_3(N+1)*w(N+1)*Rphi(1)*EVr_Nb; 
    
                    Vr_Nc = 0;
                    for n = 1:N+1
                        EVr_Nc = 0;
                        for q=1:Q
                            EVr_Nc = EVr_Nc + Xphi(n,q)*Yphi(N+1,q);
                        end
                        Vr_Nc = Vr_Nc + alpha_2(N+1)*w(n)*Rphi(n)*EVr_Nc;
                    end 
    
                    Vr_N = Vr_Na + Vr_Nb + Vr_Nc;
                    Vr_N = weakBC*Vr_N;
    
                    %RHS computation for first system
                    Vr1= b(1)*B'*Rphi + b(2)*DD'*Rphi - Pphi'*Rphi;

                    Vr1(1) = Vr1(1) - Vr_0;
                    for d = 1:N-1
                        Vr1(d+1) = Vr1(d+1) - Vr_k(d);
                    end
                    Vr1(N+1) = Vr1(N+1) - Vr_N;
    
                    Vr2 = B4'*Ru*ep - Pu'*Ru;
                    Vr3 = D3'*Rv*ep - Pv'*Rv ;
    
                    Vector2 = [ Vr1 ; Vr2 ; Vr3 ] ;
    
                    %Use first system to compute S values = LHS \ RHS
                    S1 = Matrix2\Vector2; 
                    Sphi = S1(1:N+1);
                    Su = S1(N+2:2*N);
                    Sv = S1(2*N+1:3*N+1);
    
                    %LHS computation for second system
                    %weak BC terms
                    djdRphi_0a = alpha_2(1)*(Sphi(N+1))^2 *w(1);
                    djdRphi_0b = alpha_4(1)*(Sphi(1))^2 *w(1);
    
                    djdRphi_0c = 0;
                    for n=1:N+1
                        djdRphi_0c = djdRphi_0c + alpha_3(1)*w(n)*(Sphi(n))^2;
                    end 
    
                    djdRphi_0 = djdRphi_0a + djdRphi_0b + djdRphi_0c;
                    djdRphi_0 = weakBC*djdRphi_0;

                    for i =1:N-1
                          djdRphi_j(i) = alpha_2(i+1)*(Sphi(N+1))^2 *w(i+1) + alpha_4(i+1)*(Sphi(1))^2 *w(i+1); %-+
                    end
                    djdRphi_j=weakBC*djdRphi_j; 
    
                    djdRphi_Na = alpha_2(N+1)*(Sphi(N+1))^2 *w(N+1);
                    djdRphi_Nb = alpha_4(N+1)*(Sphi(1))^2 *w(N+1);
                    djdRphi_Nc = 0;
                    for n=1:N+1
                        djdRphi_Nc = djdRphi_Nc + alpha_1(N+1)*w(n)*(Sphi(n))^2;
                    end 
    
                    djdRphi_N = djdRphi_Na + djdRphi_Nb + djdRphi_Nc;
                    djdRphi_N = weakBC*djdRphi_N;
    
                    Ms11 = (Sphi'*M*Sphi)*A + (Sphi'*A*Sphi)*M + b(1)*b(1)*(Sphi'*M*Sphi)*A + b(1)*b(2)*(Sphi'*L*Sphi)*L'+b(2)*b(1)*(Sphi'*L'*Sphi)*L + b(2)*b(2)*(Sphi'*A*Sphi)*M;
    
                    Ms11(1,1) = Ms11(1,1) + djdRphi_0;
                    for d=1:N-1
                        Ms11(d+1,d+1) = Ms11(d+1,d+1) + djdRphi_j(d); %only adding to diag
                    end
                    Ms11(N+1,N+1) = Ms11(N+1,N+1) + djdRphi_N;
    
                    Ms12 = (Sphi'*M4*Su)*L'           + b(1)*(Sphi'*M4*Su)*A*ep    + b(2)*(Sphi'*L3'*Su)*L*ep;  
                    Ms13 = (Sphi'*L'*Sv)*M4           + b(1)*(Sphi'*L*Sv)*L3'*ep   + b(2)*(Sphi'*A*Sv)*M4*ep;       
                    Ms21 = (Su'*M3*Sphi)*L            + b(1)*(Su'*M3*Sphi)*A*ep    + b(2)*(Su'*L3*Sphi)*L'*ep;
                    Ms22 = (Su'*M2*Su)  *A   *ep *ep  +      (Su'*A2*Su)*M              + (Su'*M2*Su)*M;         
                    Ms23 = (Su'*L3*Sv)  *L3' *ep *ep  -      (Su'*L4'*Sv)*L4;                                         
                    Ms31 = (Sv'*L*Sphi) *M3           + b(1)*(Sv'*L'*Sphi)*L3*ep   + b(2)*(Sv'*A*Sphi)*M3*ep;
                    Ms32 = (Sv'*L3'*Su) *L3  *ep *ep  - (Sv'*L4*Su)*L4';                                     
                    Ms33 = (Sv'*A*Sv)   *M2  *ep *ep  + (Sv'*M*Sv)*M2              + (Sv'*M*Sv)*A2;          
                    Matrix1 = [Ms11,Ms12,Ms13;Ms21,Ms22,Ms23;Ms31,Ms32,Ms33];
     
                    EVs_0a = 0;
                        for q=1:Q
                            EVs_0a = EVs_0a + Xphi(1,q)*Yphi(N+1,q);
                        end
                    Vs_0a = alpha_2(1)*EVs_0a*w(1)*Sphi(N+1); 
    
                    EVs_0b = 0;
                        for q=1:Q
                            EVs_0b = EVs_0b + Xphi(1,q)*Yphi(1,q);
                        end
                    Vs_0b = alpha_4(1)*EVs_0b*w(1)*Sphi(1); 

                    Vs_0c = 0;
                    for n = 1:N+1
                        EVs_0c = 0;
                        for q=1:Q
                            EVs_0c = EVs_0c + Xphi(1,q)*Yphi(n,q);
                        end
                        Vs_0c = Vs_0c + alpha_3(1)*w(n)*Sphi(n)*EVs_0c;
                    end
    
                    Vs_0 = Vs_0a + Vs_0b + Vs_0c;
                    Vs_0 = weakBC*Vs_0;
    
                    for j =1:N-1      
                        EVs_ja = 0;
                        for q=1:Q
                            EVs_ja = EVs_ja + Xphi(j+1,q)*Yphi(N+1,q);
                        end
                        Vs_ja(j) = EVs_ja*w(j+1)*alpha_2(j+1); 
                    end
                    Vs_ja = Vs_ja*Sphi(N+1);

                    for j =1:N-1    
                        EVs_jb = 0;
                        for q=1:Q
                            EVs_jb = EVs_jb + Xphi(j+1,q)*Yphi(1,q);
                        end
                        Vs_jb(j) = EVs_jb*w(j+1)*alpha_4(j+1);
                    end
                    Vs_jb = Vs_jb*Sphi(1);
    
                    Vs_j = Vs_ja + Vs_jb;
                    Vs_j = weakBC*Vs_j;
    
                    EVs_Na = 0;
                        for q=1:Q
                            EVs_Na = EVs_Na + Xphi(N+1,q)*Yphi(N+1,q);
                        end  
                    Vs_Na = alpha_2(N+1)*w(N+1)*Sphi(N+1)*EVs_Na; %was a2 should be a1
    
                    EVs_Nb = 0;
                        for q=1:Q
                            EVs_Nb = EVs_Nb + Xphi(N+1,q)*Yphi(1,q);
                        end  
                    Vs_Nb = alpha_4(N+1)*w(N+1)*Sphi(1)*EVs_Nb; %was a4 should be a3

                    Vs_Nc = 0;
                    for n = 1:N+1
                        EVs_Nc = 0;
                        for q=1:Q
                            EVs_Nc = EVs_Nc + Xphi(N+1,q)*Yphi(n,q);
                        end
                        Vs_Nc = Vs_Nc + alpha_1(N+1)*w(n)*Sphi(n)*EVs_Nc; 
                    end
                    Vs_Nc = Vs_Nc;
    
                    Vs_N = Vs_Na + Vs_Nb + Vs_Nc;
                    Vs_N = weakBC*Vs_N;
    
                    %RHS computation for first system
                    Vs1 = b(1)*B*Sphi + b(2)*DD*Sphi - Pphi*Sphi; 
    
                    Vs1(1) = Vs1(1) - Vs_0;
                    for d = 1:N-1
                        Vs1(d+1) = Vs1(d+1) - Vs_j(d);
                    end
                    Vs1(N+1) = Vs1(N+1) - Vs_N;
    
                    Vs2 = B4*Su*ep-Pu*Su;
                    Vs3 = D3*Sv*ep-Pv*Sv;

                    Vector1 = [ Vs1 ; Vs2 ; Vs3 ] ;
    
    
                    %USe second system for find R values = LHS \ RHS
                    R1 = Matrix1\Vector1 ;
                    Rphi = R1 (1 :N+1) ;
                    Ru = R1(N+2: 2*N+2) ;
                    Rv = R1(2*N+3:3*N+1);
    
                    %count and cap number of iterations
                    counter=counter+1;

                    if counter == max_iterations
                        break
                    end
     
                    NP1=N+1;
                
                    Alphau = 0;
                    Alphav = 0;
                    Alphaphi = 0;
                    
                    for i=1:N-1
                        Alphav = Alphav + (Rv(i)*Rv(i)*w(i+1));
                        Alphaphi = Alphaphi + (Rphi(i)*Rphi(i)*w(i+1));
                    end
                    for i =1:NP1
                        Alphau = Alphau + (Ru(i)*Ru(i)*w(i));
                    end
                    
                    Alphau = Alphau^0.5;
                    Alphav = Alphav^0.5;
                    Alphaphi = Alphaphi^0.5;
                    
                    Betau = 0;
                    Betav = 0;
                    Betaphi = 0;
                    
                    for i=1:N-1
                        Betau = Betau + (Su(i)*Su(i)*w(i+1));
                        Betaphi = Betaphi + (Sphi(i)*Sphi(i)*w(i+1));
                    end
                    for i =1:NP1
                        Betav = Betav + (Sv(i)*Sv(i)*w(i));
                    end
                    
                    Betau = Betau^0.5;
                    Betav = Betav^0.5;
                    Betaphi = Betaphi^0.5;
                    
                    Ru = Ru*(Betau/Alphau);
                    Rv = Rv*(Betav/Alphav);
                    Rphi = Rphi*(Betaphi/Alphaphi);
                    Su = Su*(Alphau/Betau);
                    Sv = Sv*(Alphav/Betav);
                    Sphi = Sphi*(Alphaphi/Betaphi);
                    
                    R1 = [Rphi; Ru; Rv];
                    S1 = [Sphi; Su; Sv];


                end

                %setting XQ YQ from the calculated values
                xphi = Rphi ;
                xu= Ru;
                xv = Rv ;
                yphi = Sphi ;
                yu = Su ;
                yv = Sv ; 

                %adding the known boundary v al u e s
                Xphi = [ Xphi , [  xphi ] ] ;
                Xu = [ Xu, xu ] ;
                Xv = [ Xv, [ 0 ; xv ; 0 ] ] ;
                Yphi = [ Yphi , [ yphi ] ] ;
                Yu = [ Yu , [ 0 ; yu ; 0 ] ] ;
                Yv = [ Yv, yv ] ;
                Q= Q+1; 

                %calculating the approximations
                approxphi = zeros(N+1) ;
                approxu = zeros(N+1) ;
                approxv = zeros(N+1) ;
                for k = 1 :N+1
                    for l = 1 :N+1
                        for j = 1 : Q
                            approxphi(l,k) = approxphi(l,k)+Yphi(k,j)*Xphi(l,j);
                            approxu (l,k ) = approxu ( l,k ) + Yu( k , j ) *Xu( l , j ) ;
                            approxv (l,k ) = approxv ( l,k ) + Yv( k , j ) *Xv( l , j ) ;
                        end
                    end
                end

                %working out the error
                E = 0 ;
                Emax = 0 ;
                PHIerror= zeros(N+1) ;
                PHImesh= zeros(N+1) ;
                Umesh = zeros(N+1) ;
                Vmesh = zeros(N+1) ;
                Uerror = zeros(N+1) ;
                Verror = zeros(N+1) ;
    
                for k = 1 :N+1
                    for l = 1 :N+1
                        err = approxphi(l,k)- phi(node(l),node(k));
                        Uerror(l,k) = approxu(l,k) - U(node(l),node(k));
                        Verror(l,k) = approxv(l,k) - V(node(l),node(k));
                        PHIerror(l,k) = approxphi(l,k) - phi(node(l),node(k));
                        PHImesh(l,k) = phi(node(k),node(l));
                        Umesh(l,k) = U(node(l),node(k));
                        Vmesh(l,k) = V(node(l),node(k)); 
                    end
                end
         
                QL2phierror=0;
                QL2uerror=0;
                QL2verror=0;
            
                for num=1:N1
                    for m=1:N1
                        QL2phierror=QL2phierror+((w(m)*w(num)*(PHIerror(m,num))^2)) ;
                        QL2uerror=QL2uerror+((w(m)*w(num)*(Uerror(m,num))^2)) ;
                        QL2verror=QL2verror+((w(m)*w(num)*(Verror(m,num))^2)) ;
                    end
                end
                QL2phierror= QL2phierror^(1/2);
                QL2uerror= QL2uerror^(1/2);
                QL2verror= QL2verror^(1/2);
            
                QL2phi(Q,number)=QL2phierror;
        
            
            end
            Maximum_Error(number)=Emax; 
            RMS_Error(number)=E;
            
            L2phierror=0;
            L2Uerror=0;
            L2Verror=0;
            for num=1:N1
                for m=1:N1
                   L2Uerror=L2Uerror+((w(m)*w(num)*(Uerror(m,num))^2)) ;
                   L2Verror=L2Verror+((w(m)*w(num)*(Verror(m,num))^2)) ;
                   L2phierror=L2phierror+((w(m)*w(num)*(PHIerror(m,num))^2)) ;
                end
            end
            L2Uerror= L2Uerror^(1/2);
            L2Verror= L2Verror^(1/2);
            L2phierror= L2phierror^(1/2);
            
            L2phi(number)=L2phierror;
            L2U(number)=L2Uerror;
            L2V(number)=L2Verror;
    
            NBL2phi(number,bi,ei)=L2phierror;
        
        end
        
    end
    
    N_values=Nvalues';


    figure('Units', 'normalized', 'Position', [0.2 0.2 0.5 0.5]); 
    set(gcf, 'Color', 'white');
    handles = semilogy(N_values,NBL2phi(:,:,ei), 'Linewidth', 2);
    legend('$\mathbf{b}=(0,1)^T$','$\mathbf{b}=(1,0)^T$','$\mathbf{b}=(0.01,1)^T$', '$\mathbf{b}=(1,0.01)^T$', '$\mathbf{b}=(0.1,1)^T$', '$\mathbf{b}=(1,0.1)^T$', 'Interpreter', 'latex','FontSize', 14 )
    grid on
    xlabel('$N$','interpreter','latex','fontsize',20)
    xlim([N_values(1) N_values(end)])
    xticks(N_values)
    ylabel('$L^{2}$ Error in $\phi$','interpreter','latex','fontsize',20)
    % title(sprintf('$\\bar{\\epsilon} = %.3f$', ep), ...
    %      'Interpreter', 'latex', 'fontsize', 14);
    set(gca,'fontsize',14);

end
 

A = squeeze(NBL2phi(end, :, :)); 
rowLabels = { ...
    'b=(0,1)^T', ...
    'b=(1,0)^T', ...
    'b=(0.01,1)^T', ...
    'b=(1,0.01)^T', ...
    'b=(0.1,1)^T', ...
    'b=(1,0.1)^T'};

T = array2table(A, ...
    'VariableNames', {'epsilon bar = 1', 'epsilon bar = 0.001'}, ...  
    'RowNames', rowLabels);
fprintf('\n');
fprintf('L2 errors for φ for Example 6 using the SLS SM PGD, after one enrichment with N=24:')
fprintf('\n');
disp(T)

