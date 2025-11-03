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

%% USER INPUT
%  Depending which figure you want to replicate, set plot equal to one of
%  the following

% plot = "Figure 5.6 (a)"; 
% plot = "Figure 5.6 (b)"; 
% plot = "Figure 5.6 (c)";  
% plot = "Figure 5.6 (d)";
% plot = "Figure 5.9 (a)"; 
% plot = "Figure 5.9 (b)"; 
% plot = "Figure 5.9 (c)";
plot = "Figure 5.9 (d)";


if plot == "Figure 5.6 (a)"

    b_values = [0, 1]; 
    eps_values = 1;  
    ax = [-1, 1]; % spectral method
    cy = [-1, 1]; 
    N=48;
    
elseif plot == "Figure 5.6 (b)"

    b_values = [0, 1]; 
    eps_values = 0.1;  
    ax = [-1,  1];
    cy = [-1 , 1]; 
    N=48;

elseif plot == "Figure 5.6 (c)"
     
    b_values = [0, 1]; 
    eps_values = 0.01;
    ax = [-1 ,-0.85,0.85,  1];
    cy = [-1 ,-0.999,0.9999, 1];   
    N=48; 

elseif plot == "Figure 5.6 (d)" 
    
    b_values = [0,1]; 
    eps_values = 0.001;
    cy = [-1 ,-0.99,0.99, 1]; 
    ax = [-1 ,-0.93 ,0.93, 1];  
    N=48;

    elseif plot == "Figure 5.9 (a)"

    b_values = [1, 1]; 
    eps_values = 1;  
    ax = [-1, 1];
    cy = [-1, 1]; 
    N=48;

    elseif plot == "Figure 5.9 (b)"

    b_values = [1, 1]; 
    eps_values = 0.1;  
    ax = [-1, 1];
    cy = [-1, 1]; 
    N=48;

elseif plot == "Figure 5.9 (c)"

    b_values = [1, 1]; 
    eps_values = 0.01;  
    ax = [-1, 0.95, 1];
    cy = [-1, 0.95, 1]; 
    N=48;
    
elseif plot == "Figure 5.9 (d)"

    b_values = [1, 1]; 
    eps_values = 0.001;  
    ax = [-1, 0.99, 1];
    cy = [-1, 0.99, 1]; 
    N=48;

else
end
 
% tolerance
epsilon = 1e-8;

% Maximum number of enrichments
no_enrichment = 3; 

%Maximum number of iterations in the fixed point loop
stopcrit=20;
Qerrors=[0,0,0];

aa=-1;
bb=1;
cc=-1;
dd=1;
syms x y
executionTimes = [];
 
for bi = 1:size(b_values,1)
    b = [b_values(bi,1) b_values(bi,2)];

    for ei = 1:size(eps_values,1)
        ep = eps_values(ei);
    
        Kx=length(ax)-1;
        Ky=length(cy)-1; 
        
        % Truncation + 1
        NN=N-1;
        N1=N;
     
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
        A = zeros(N1); MM = zeros(N1); L = zeros(N1);
        
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
            MM(i,i)=w(i);
        end
        for i=1:N1
            for j=1:N1
                L(i,j)=(w(i).*D(i,j));
            end
        end
        Ay=A;
        Ax=A;
        Mx=MM;
        My=MM;
        Ly=L;
        Lx=L;
         
        dx=zeros(Kx,1)';
        dy=zeros(Ky,1)'; 
        localFx=zeros(1,N);
        localFy=zeros(1,N);
        % Setting up the global matrices
        GlobalAx=zeros((N)*Kx-(Kx-1),(N)*Kx-(Kx-1));
        GlobalMx=zeros((N)*Kx-(Kx-1),(N)*Kx-(Kx-1));
        GlobalLx=zeros((N)*Kx-(Kx-1),(N)*Kx-(Kx-1));
        GlobalAy=zeros((N)*Ky-(Ky-1),(N)*Ky-(Ky-1));
        GlobalMy=zeros((N)*Ky-(Ky-1),(N)*Ky-(Ky-1));
        GlobalLy=zeros((N)*Ky-(Ky-1),(N)*Ky-(Ky-1));
        
        localDiffMatrix=zeros(N,N);
        GlobalDiff=zeros((N)*Ky-(Ky-1),(N)*Ky-(Ky-1));
%========== Local and Global matrix computations ========================
%========================================================================
    
        physicalx=[];
        for i=1:Kx %for each element in x
        
            %compute delta x
            dx(i)= ax(i+1)-ax(i);
            
            %map to physical points in each element
            physx=ax(i)+(xx+1)*(ax(i+1)-ax(i))/2; % this is called xi invverse in write up
            
            %all points in physical domain x
            physicalx=[physicalx; physx];
    
            xdomain =physicalx;
            for ik = 1:(Kx - 1)
                index_to_delete = ((N * ik) - (ik - 1));
    
                % Check if the index is within the vector bounds
                if index_to_delete <= length(xdomain)
                    xdomain(index_to_delete) = [];
                end
            end
             
            %Local Stiffness matrix for each element in x
            localAx=2/(dx(i))*Ax;
            %Global Stiffness matrix for x
            GlobalAx([(i-1)*N-(i-2):(i)*N-(i-1)],[(i-1)*N-(i-2):(i)*N-(i-1)]) = GlobalAx([(i-1)*N-(i-2):(i)*N-(i-1)],[(i-1)*N-(i-2):(i)*N-(i-1)]) + localAx;
            
            %Local Mass matrix for each element in x
            localMx=(dx(i))/2*Mx;
            %Global Mass matrix for x
            GlobalMx([(i-1)*N-(i-2):(i)*N-(i-1)],[(i-1)*N-(i-2):(i)*N-(i-1)])=GlobalMx([(i-1)*N-(i-2):(i)*N-(i-1)],[(i-1)*N-(i-2):(i)*N-(i-1)]) + localMx;    
            
            %Local L matrix for each element in x
            localLx = Lx;
            %Global Mass matrix for x
            GlobalLx([(i-1)*N-(i-2):(i)*N-(i-1)],[(i-1)*N-(i-2):(i)*N-(i-1)])=GlobalLx([(i-1)*N-(i-2):(i)*N-(i-1)],[(i-1)*N-(i-2):(i)*N-(i-1)]) + localLx;    
              
        end
    
        physicaly=[];
        for i=1:Ky %for each element in y
        
            %compute delta y
            dy(i)= cy(i+1)-cy(i);
            
            %map to physical points in each element
            physy=cy(i)+(xx+1)*(cy(i+1)-cy(i))/2;
            
            physicaly=[physicaly; physy];
            ydomain =physicaly;
              
            for ik = 1:(Ky - 1)
                index_to_delete = ((N * ik) - (ik - 1));
    
                % Check if the index is within the vector bounds
                if index_to_delete <= length(ydomain)
                    ydomain(index_to_delete) = [];
                end
            end
            %Local Stiffness matrix for each element in y
            localAy=2/(dy(i))*Ay;
            %Global Stiffness matrix for y
            GlobalAy([(i-1)*N-(i-2):(i)*N-(i-1)],[(i-1)*N-(i-2):(i)*N-(i-1)]) = GlobalAy([(i-1)*N-(i-2):(i)*N-(i-1)],[(i-1)*N-(i-2):(i)*N-(i-1)]) + localAy;
            
            %Local Mass matrix for each element in y
            localMy=(dy(i))/2*My;
            %Global Mass matrix for y
            GlobalMy([(i-1)*N-(i-2):(i)*N-(i-1)],[(i-1)*N-(i-2):(i)*N-(i-1)]) = GlobalMy([(i-1)*N-(i-2):(i)*N-(i-1)],[(i-1)*N-(i-2):(i)*N-(i-1)]) + localMy;    
            
            %Local L matrix for each element in y
            localLy=Ly;
            %Global L matrix for y
            GlobalLy([(i-1)*N-(i-2):(i)*N-(i-1)],[(i-1)*N-(i-2):(i)*N-(i-1)]) = GlobalLy([(i-1)*N-(i-2):(i)*N-(i-1)],[(i-1)*N-(i-2):(i)*N-(i-1)]) + localLy;    
    
        end
    
        GlobalBx= zeros((N)*Kx-(Kx-1));
        GlobalCx= zeros((N)*Kx-(Kx-1));
        Bf=zeros(N1);
        Cf=zeros(N1);
    
        for i=1:Kx %for each element in x
            %compute delta x
            dx(i)= ax(i+1)-ax(i);
                  
            %map to physical points in each element
            physx=ax(i)+(xx+1)*(ax(i+1)-ax(i))/2;
         
            for k=1:Ky   
                %compute delta y
                dy(k)= cy(k+1)-cy(k);
                
                %map to physical points in each element
                physy=cy(k)+(xx+1)*(cy(k+1)-cy(k))/2;
                
                Fxy = zeros(N, N);  
                for iw = 1:N
                    for jw = 1:N
                        Fxy(iw,jw) = 1;
                    end
                end
                  
                for ii=1:N1
                    for jj=1:N1
                    summ=0;
                        for nn=1:N1
                            summ=summ+(w(nn)*D(nn,ii)*Fxy(nn,jj));
                        end
                        Bf(ii,jj)=(dy(k))/2*w(jj)*summ;
                    end
                end   
        
                for ii=1:N1
                    for jj=1:N1
                        suum=0;
                        for nn=1:N1
                            suum=suum+w(nn)*Fxy(ii,nn)*D(nn,jj);
                        end
                        Cf(ii,jj)=(dx(i))/2*w(ii)*suum;
                    end
                end 
                Bx=Bf;
                Cx=Cf;
                
                %Local B matrix for each element in x
                localBx = Bx;
                %Global B matrix for x
                GlobalBx([(i-1)*N-(i-2):(i)*N-(i-1)],[(k-1)*N-(k-2):(k)*N-(k-1)])=GlobalBx([(i-1)*N-(i-2):(i)*N-(i-1)],[(k-1)*N-(k-2):(k)*N-(k-1)]) + localBx;    
                
                %Local C matrix for each element in x
                localCx=Cx;
                %Global C matrix for x
                GlobalCx([(i-1)*N-(i-2):(i)*N-(i-1)],[(k-1)*N-(k-2):(k)*N-(k-1)])=GlobalCx([(i-1)*N-(i-2):(i)*N-(i-1)],[(k-1)*N-(k-2):(k)*N-(k-1)]) + localCx;    
            end 
        end
    
        B=GlobalBx;
        C=GlobalCx;
 
%========================= PGD begins ====================================
%=========================================================================

        My=GlobalMy;
        Mx=GlobalMx;
        Ay=GlobalAy;
        Ax=GlobalAx;
        Ly=GlobalLy;
        Lx=GlobalLx;
     
        A=Ax;
        M=Mx;
        L=Lx;
        
        x=xx;
        y=xx; 
        N0 = size(N*Kx-(Kx-1),2);
        X=zeros(N*Kx-(Kx-1),1);
        Y=zeros(N*Kx-(Kx-1),1);
        
        %Start off Q
        Q=0;
        
        %Start approximations
        Xphi = zeros(N*Kx-(Kx-1),Q ) ;
        Xu = zeros(N*Kx-(Kx-1),Q ) ;
        Xv = zeros(N*Kx-(Kx-1),Q ) ;
        Yphi = zeros(N*Kx-(Kx-1),Q ) ;
        Yu = zeros(N*Kx-(Kx-1),Q ) ;
        Yv = zeros(N*Kx-(Kx-1),Q ) ;

        %Qth enrichments
        for term = 1:no_enrichment
            startTime = cputime;
        
            %Initialize the Matrices containing the known Q−1 enrichments
            Pphi = zeros(N*Kx-(Kx-1)) ;
            Pu = zeros(N*Kx-(Kx-1),N*Kx-(Kx-1)) ;
            Pv = zeros(N*Kx-(Kx-1),N*Kx-(Kx-1)) ;
    
            XphiA = zeros(N*Kx-(Kx-1),Q ) ;
            XphiM= zeros(N*Kx-(Kx-1),Q ) ;
            XphiL= zeros(N*Kx-(Kx-1),Q ) ;
            XphiLT =zeros(N*Kx-(Kx-1),Q);
            XuA= zeros(N*Kx-(Kx-1),Q ) ;
            XuM= zeros(N*Kx-(Kx-1),Q ) ;
            XuL= zeros(N*Kx-(Kx-1),Q ) ;
            XuLT = zeros(N*Kx-(Kx-1),Q ) ;
            XvA= zeros(N*Kx-(Kx-1),Q ) ;
            XvM= zeros(N*Kx-(Kx-1),Q ) ;
            XvL= zeros(N*Kx-(Kx-1),Q ) ;
            XvLT = zeros(N*Kx-(Kx-1),Q ) ;
    
            YphiA= zeros(N*Kx-(Kx-1),Q ) ;
            YphiM= zeros(N*Kx-(Kx-1),Q ) ;
            YphiL= zeros(N*Kx-(Kx-1),Q ) ;
            YphiLT=zeros(N*Kx-(Kx-1), Q);
            YuA= zeros(N*Kx-(Kx-1),Q ) ;
            YuM= zeros(N*Kx-(Kx-1),Q ) ;
            YuL= zeros(N*Kx-(Kx-1),Q ) ;
            YuLT = zeros(N*Kx-(Kx-1),Q ) ;
            YvA= zeros(N*Kx-(Kx-1),Q ) ;
            YvL= zeros(N*Kx-(Kx-1),Q ) ;
            YvLT = zeros(N*Kx-(Kx-1),Q ) ;
            YvM= zeros(N*Kx-(Kx-1),Q ) ;
    
            for i = 1 :N*Kx-(Kx-1)
                for j = 1 : Q
                    for k = 1 :N*Kx-(Kx-1)
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
        
            for i = 1 :N*Kx-(Kx-1)
                for j = 1 : Q
                    for k = 1 :N*Kx-(Kx-1)                                                                                                           %here                        %here                                                                                                   %here                        %here                                                                                                                                                
                        Pphi(k,i) = Pphi(k,i) + (XuLT(k,j)*YuM(i,j) + XphiA(k,j)*YphiM(i,j) + XvM(k,j)*YvLT(i,j) + XphiM(k,j)*YphiA(i,j) + ep*b(1)*XuA(k,j)*YuM(i,j) + ep*b(1)*XvLT(k,j)*YvL(i,j) + b(1)*b(1)*XphiA(k,j)*YphiM(i,j) + b(1)*b(2)*XphiLT(k,j)*YphiL(i,j) + ep*b(2)*XuL(k,j)*YuLT(i,j) + ep*b(2)*XvM(k,j)*YvA(i,j) + b(2)*b(1)*XphiL(k,j)*YphiLT(i,j) + b(2)*b(2)*XphiM(k,j)*YphiA(i,j));
                        Pu(k,i) = Pu(k,i) + (ep*b(1)*XphiA(k,j)*YphiM(i,j)  + ep*b(2)*XphiLT(k,j)*YphiL(i,j) + ep*ep*XuA(k,j)*YuM(i,j) + ep*ep*XvLT(k,j)*YvL(i,j) + XuM(k,j)*YuM(i,j) + XphiL(k,j)*YphiM(i,j) - XvL(k,j)*YvLT(i,j) + XuM(k,j)*YuA(i,j));
                        Pv(k,i) = Pv(k,i) + (ep*b(1)*XphiL(k,j)*YphiLT(i,j) + ep*b(2)*XphiM(k,j)*YphiA(i,j)  + ep*ep*XvM(k,j)*YvA(i,j) + ep*ep*XuL(k,j)*YuLT(i,j) + XvM(k,j)*YvM(i,j) + XphiM(k,j)*YphiL(i,j) - XuLT(k,j)*YuL(i,j) + XvA(k,j)*YvM(i,j));
                    end
                end
            end
        
            %Unknown functions for enrichment
            oldR = zeros( 3*((N)*Kx-(Kx-1)),1) ;
            oldS = zeros( 3*((N)*Kx-(Kx-1)),1) ;
            
            Rphi = rand(N*Kx-(Kx-1) ,1) ;
            Ru = rand(N*Kx-(Kx-1), 1 ) ;
            Rv = rand(N*Kx-(Kx-1), 1 ) ;
            if term >1
            Rphi = Xphi(:,term-1);
            Ru = Xu(:,term-1);
            Rv = Xv(:,term-1);
            
            else
            end
    
            counter =0;
            %make an initial guess
            R1 = randn(3*((N)*Kx-(Kx-1)),1);
            S1 = randn(3*((N)*Kx-(Kx-1)),1);
            Rphi(1) = 0;
            Rphi(end)=0;
            Rv(1) = 0;
            Rv(end)=0;
            Sphi(1) = 0;
            Sphi(end)=0;
            Su(1) = 0;
            Su(end)=0;
     
            %fixed point iterations
            while max(max( abs (R1*S1' - oldR* oldS') ) ) >epsilon     
                oldR = R1; 
                oldS = S1; 
                
                % solve for the first matrix system
                %LHS
                M2_11 = (Rphi'*A*Rphi)*M + (Rphi'*M*Rphi)*A + b(1)*b(1)*(Rphi'*A*Rphi)*M + b(1)*b(2)*(Rphi'*L'*Rphi)*L + b(2)*b(1)*(Rphi'*L*Rphi)*L'+ b(2)*b(2)*(Rphi'*M*Rphi)*A ;
                M2_12 = (Rphi'*L'*Ru)*M + b(1)*ep*(Rphi'*A*Ru)*M+ b(2)*ep*(Rphi'*L*Ru)*L';
                M2_13 = (Rphi'*M*Rv)*L' + b(1)*ep*(Rphi'*L'*Rv)*L + b(2)*ep*(Rphi'*M*Rv)*A;
                M2_21 = (Ru'*L*Rphi)*M + b(1)*ep*(Ru'*A*Rphi)*M + b(2)*ep*(Ru'*L'*Rphi)*L;
                M2_22 = (Ru'*A*Ru)*M*ep*ep + (Ru'*M*Ru)*A + (Ru'*M*Ru)*M;
                M2_23 = (Ru'*L'*Rv)*L*ep*ep - (Ru'*L*Rv)*L';
                M2_31 = (Rv'*M*Rphi)*L + b(1)*ep*(Rv'*L*Rphi)*L' + b(2)*ep*(Rv'*M*Rphi)*A;
                M2_32 = (Rv'*L*Ru)*L'*ep*ep - (Rv'*L'*Ru)*L;
                M2_33 = (Rv'*A*Rv)*M + (Rv'*M*Rv)*A*ep*ep + (Rv'*M*Rv)*M; 
    
                Matrix2 = [M2_11, M2_12, M2_13; M2_21, M2_22, M2_23; M2_31, M2_32, M2_33] ;
                % RHS
                V2_1 =  - Pphi'*Rphi +b(1)*B'*Rphi + b(2)*C'*Rphi;
                V2_2 = B'*Ru*ep-Pu'*Ru;
                V2_3 = C'*Rv*ep-Pv'*Rv ;
                Vector2 = [ V2_1 ; V2_2 ; V2_3 ] ;
    
                Matrix2(2*(N*Kx-(Kx-1)),:) = [];
                Matrix2(N*Kx-(Kx-2),:) = [];
                Matrix2(N*Kx-(Kx-1),:) = [];
                Matrix2(1,:) = [];
                
                Matrix2(:,2*(N*Kx-(Kx-1))) = [];
                Matrix2(:,N*Kx-(Kx-2)) = [];
                Matrix2(:,N*Kx-(Kx-1)) = [];
                Matrix2(:,1) = [];
                
                Vector2(2*(N*Kx-(Kx-1))) = [];
                Vector2(N*Kx-(Kx-2)) = [];
                Vector2(N*Kx-(Kx-1)) = [];
                Vector2(1) = [];
                
                S1 = Matrix2\Vector2;
                Sphi = S1(1:N*Kx-(Kx-1)-2);
                Su = S1(N*Kx-(Kx-2)-2:2*(N*Kx-(Kx-1))-4);
                Sv = S1(2*(N*Kx-(Kx-1))+1 - 4: end);
                
                Sphi = [0;Sphi;0];
                Su= [0;Su;0];
    
                S1 = [Sphi;Su;Sv]; 
                
                % Solve for the second system
                % LHS
                M1_11 = (Sphi'*M*Sphi)*A+(Sphi'*A*Sphi)*M + b(1)*b(1)*(Sphi'*M*Sphi)*A + b(1)*b(2)*(Sphi'*L*Sphi)*L' + b(2)*b(1)*(Sphi'*L'*Sphi)*L + b(2)*b(2)*(Sphi'*A*Sphi)*M;
                M1_12 = (Sphi'*M*Su)*L' + b(1)*ep*(Sphi'*M*Su)*A + b(2)*ep*(Sphi'*L'*Su)*L;
                M1_13 = (Sphi'*L'*Sv)*M + b(1)*ep*(Sphi'*L*Sv)*L'+b(2)*ep*(Sphi'*A*Sv)*M;
                M1_21 = (Su'*M*Sphi)*L + b(1)*ep*(Su'*M*Sphi)*A + b(2)*ep*(Su'*L*Sphi)*L';
                M1_22 = (Su'*M*Su)*A*ep*ep + (Su'*A*Su)*M + (Su'*M*Su)*M;
                M1_23 = (Su'*L*Sv)*L'*ep*ep - (Su'*L'*Sv)*L;
                M1_31 = (Sv'*L*Sphi)*M + b(1)*ep*(Sv'*L'*Sphi)*L + b(2)*ep*(Sv'*A*Sphi)*M;
                M1_32 = (Sv'*L'*Su)*L*ep*ep-(Sv'*L*Su)*L';
                M1_33 = (Sv'*A*Sv)*M*ep*ep+(Sv'*M*Sv)*M+(Sv'*M*Sv)*A;
                Matrix1 = [M1_11,M1_12,M1_13;M1_21,M1_22,M1_23;M1_31,M1_32,M1_33];
                
                % RHS
                V1_1 = -Pphi*Sphi + b(1)*B*Sphi + b(2)*C*Sphi;
                V1_2 = -Pu*Su + B*Su*ep;
                V1_3 = -Pv*Sv + C*Sv*ep;
                
                Vector1 = [V1_1 ; V1_2 ; V1_3 ] ;
                 
                Matrix1(end,:) = [];
                Matrix1(2*(N*Kx-(Kx-1))+1,:) = [];
                Matrix1(N*Kx-(Kx-1),:) = [];
                Matrix1(1,:) = [];
                Matrix1(:,end) = [];
                Matrix1(:,2*(N*Kx-(Kx-1))+1) = [];
                Matrix1(:,N*Kx-(Kx-1)) = [];
                Matrix1(:,1) = [];
                
                Vector1(end) = [];
                Vector1(2*(N*Kx-(Kx-1))+1) = [];
                Vector1(N*Kx-(Kx-1)) = [];
                Vector1(1) = [];
                
                R1 = Matrix1\Vector1;
                Rphi = R1(1:N*Kx-(Kx-1)-2);
                Ru = R1(N*Kx-(Kx-2)-2:2*(N*Kx-(Kx-1))-2);
                Rv = R1(2*(N*Kx-(Kx-1))+1 - 2: end);
                
                Rphi = [0;Rphi;0];
                Rv= [0;Rv;0];
                R1 = [Rphi;Ru;Rv];
              
                Alphau = 0;
                Alphav = 0;
                Alphaphi = 0;
                
                for i=1:(N*Kx-(Kx-1))
                    Alphav = Alphav + (Rv(i)*Rv(i));
                    Alphaphi = Alphaphi + (Rphi(i)*Rphi(i));
                    Alphau = Alphau + (Ru(i)*Ru(i));
                end
                    
                Alphau = Alphau^0.5;
                Alphav = Alphav^0.5;
                Alphaphi = Alphaphi^0.5;
                
                Ru = Ru*(1/Alphau);
                Rv = Rv*(1/Alphav);
                Rphi = Rphi*(1/Alphaphi);
                R1=[Rphi;Ru;Rv];
                
                counter=counter+1 ;
                if counter == stopcrit
                    break 
                end
                
            end 
       
            %setting XQ YQ from the calculated values
            xphi = Rphi ;
            xu= Ru;
            xv = Rv ;
            yphi = Sphi ;
            yu = Su ;
            yv = Sv ;
    
            %adding the known boundary values
            Xphi = [ Xphi , xphi ] ;
            Xu = [ Xu, xu ] ;
            Xv = [ Xv, xv ] ;
            Yphi = [ Yphi , yphi ] ;
            Yu = [ Yu , [ yu ] ] ;
            Yv = [ Yv, yv ] ;
            Q= Q+1; 
    
            %calculating the approximations
            approxphi = zeros(N*Kx-(Kx-1), N*Ky-(Ky-1)) ;
            approxu = zeros(N*Kx-(Kx-1) ) ;
            approxv = zeros(N*Kx-(Kx-1) ) ;
            for k = 1 :N*Kx-(Kx-1) 
                for l = 1 :N*Ky-(Ky-1)
                    for j = 1 : Q
                        approxphi(k,l) = approxphi(k,l) + Yphi(l,j)*Xphi(k,j);
                        approxu(k,l) = approxu(k,l) + Yu(l,j)*Xu(k,j) ;
                        approxv(k,l) = approxv(k,l) + Yv(l,j)*Xv(k,j) ;
                    end
                end
            end
    
            approxphi_1 = zeros(N*Kx-(Kx-1) ) ;
            
            for k = 1 :N*Kx-(Kx-1) 
                for l = 1 :N*Kx-(Kx-1) 
                    for j =Q
                        approxphi_1(k,l) = Yphi(l,j)*Xphi(k,j);
                    end
                end
            end
            approxphiQ(:,:, Q) = approxphi_1;     
     
            PhiNorm=0;
            PhiNorm_all(bi, ei) = PhiNorm;
            
            for i=1:Kx
                for j=1:Ky
                    for num=1:N1
                        for m=1:N1
                           PhiNorm=PhiNorm+(ax(i+1)-ax(i))*(cy(j+1)-cy(j))/(2^(Kx))*w(m)*w(num)*(approxphi_1(m,num))^2 ;
                        end
                    end
                end
            end
    
            PhiNorm= PhiNorm^(1/2);
            QeLphi(Q,ei)=PhiNorm; 
        end 
   
    end

end
 
figure('Units', 'normalized', 'Position', [0.2 0.2 0.5 0.5]);  
surf(ydomain, xdomain, approxphi);
title(sprintf(' $\\bar{\\epsilon}=$%.3g', eps_values(1)), 'Interpreter', 'Latex');
view(30,50);
set(gca,'fontsize',12);
xlabel('$x$','interpreter','latex', 'fontsize', 20);
ylabel('$y$','interpreter','latex', 'fontsize', 20);
