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

% Number of elements in each direction
Kx=3; Ky=Kx;

% Max number of iterations
max_iterations = 100;

% Tolerance
epsilon=10^(-10);

executionTimes = [];
cumTime = zeros(6,1);

% loops over epsilon bar values
for ei = 1: size(eps_values,1)
    ep = eps_values(ei);

    if ep == 1
        epsilon=10^(-10);
    else 
        epsilon=10^(-9);
    end

    % loop over b values
    for bi = 1:size(b_values,1) 
        b = [b_values(bi,1) b_values(bi,2)];

        % loop over N values
        for Num=1:length(Nvalues) 
            cumcounter =0;
            cumcounter_values =[]; 
            each_count_value=[];
           
       
%======== Define 2 1D arrays that are the physical domain ===========
            aa=-1;
            bb=1;
            cc=-1;
            dd=1;
            % For randomly distributed elements
            ax=aa + (bb-aa)*rand(Kx+1,1);
            ax=ax';
            
            % For linearly spaced elements
            %ax=linspace(aa,bb,Kx+1);
            ax(1)=aa;
            ax(Kx+1)=bb;
            ax=sort(ax);
            
            %For randomly distributed elements
            cy=cc + (dd-cc)*rand(Ky+1,1);
            cy=cy';
            
            %For linearly spaced elements
            %cy=linspace(cc,dd,Ky+1);
            cy(1)=cc;
            cy(Ky+1)=dd;
            cy=sort(cy);
             
            ax=linspace(aa,bb,Kx+1);
            cy=linspace(cc,dd,Ky+1); 
             
            N=Nvalues(Num);
        
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
            GlobalAx=zeros((N)*Kx-(Kx-1),(N)*Kx-(Kx-1));
            GlobalMx=zeros((N)*Kx-(Kx-1),(N)*Kx-(Kx-1));
            GlobalLx=zeros((N)*Kx-(Kx-1),(N)*Kx-(Kx-1));
            GlobalBx=zeros((N)*Kx-(Kx-1),(N)*Kx-(Kx-1));
            GlobalCx=zeros((N)*Kx-(Kx-1),(N)*Kx-(Kx-1));
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
                physx=ax(i)+(xx+1)*(ax(i+1)-ax(i))/2; 
                
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
    
                for ik = 1:(Kx - 1)
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
            
                    Fxy = ep*((2)*(pi^2)*sin(pi*physx).*(sin(pi*physy'))  );
                    Fxy = Fxy + b(1)*(pi*cos(pi*physx)*sin(pi*physy')   );
                    Fxy = Fxy + b(2)*(pi*sin(pi*physx)*cos(pi*physy')  );
     
                   % Compute B and D locally     
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
                    end %here weird in last col of first el and forst col of last el... j N+1 then j 1
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
        
                localDiffMatrix = 2/(dx(i)) * D; 
                GlobalDiff([(i-1)*N-(i-2):(i)*N-(i-1)],[(i-1)*N-(i-2):(i)*N-(i-1)]) = GlobalDiff([(i-1)*N-(i-2):(i)*N-(i-1)],[(i-1)*N-(i-2):(i)*N-(i-1)]) + localDiffMatrix;
          
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
                Pu = zeros(N*Kx-(Kx-1)) ;
                Pv = zeros(N*Kx-(Kx-1)) ;
                 
                XphiA = zeros(N*Kx-(Kx-1),Q ) ;
                XphiM= zeros(N*Kx-(Kx-1),Q ) ;
                XphiL= zeros(N*Kx-(Kx-1),Q ) ;
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
                YuA= zeros(N*Kx-(Kx-1),Q ) ;
                YuM= zeros(N*Kx-(Kx-1),Q ) ;
                YuL= zeros(N*Kx-(Kx-1),Q ) ;
                YuLT = zeros(N*Kx-(Kx-1),Q ) ;
                YvA= zeros(N*Kx-(Kx-1),Q ) ;
                YvL= zeros(N*Kx-(Kx-1),Q ) ;
                YvLT = zeros(N*Kx-(Kx-1),Q ) ;
                YvM= zeros(N*Kx-(Kx-1),Q ) ;
                cXphiA = zeros (N*Kx-(Kx-1),Q ) ;
                cXphiLT = zeros (N*Kx-(Kx-1),Q ) ;
                cYphiLT = zeros (N*Kx-(Kx-1),Q ) ;
                cYphiA = zeros (N*Kx-(Kx-1),Q ) ;
                cXuA = zeros(N*Kx-(Kx-1),Q ) ;
                cXvLT = zeros(N*Kx-(Kx-1),Q ) ;
                c2XphiLT = zeros(N*Kx-(Kx-1),Q ) ;
                cXuL = zeros(N*Kx-(Kx-1),Q ) ;
                cXvM = zeros(N*Kx-(Kx-1),Q ) ;
                cXphiL = zeros(N*Kx-(Kx-1),Q ) ;
                cYuM = zeros(N*Kx-(Kx-1),Q ) ;
                cYvL = zeros(N*Kx-(Kx-1),Q ) ;
                cYphiL = zeros(N*Kx-(Kx-1),Q ) ;
                cYuLT = zeros(N*Kx-(Kx-1),Q ) ;
                cYvA = zeros(N*Kx-(Kx-1),Q ) ;
                c2YphiLT = zeros(N*Kx-(Kx-1),Q ) ;
                    
                for i = 1 :N*Kx-(Kx-1)
                    for j = 1 : Q
                        for k = 1 :N*Kx-(Kx-1)
                            XphiA ( i , j ) = XphiA ( i , j ) + Xphi ( k , j ) *A( k , i ) ;
                            XphiM( i , j ) = XphiM ( i , j ) + Xphi ( k , j ) *M( k , i ) ;
                            XuL( i , j ) = XuL( i , j ) + Xu( k , j ) *L ( i , k ) ;
                            XuLT( i , j ) = XuLT( i , j ) + Xu( k , j ) *L ( k , i ) ;
                            XvA( i , j ) = XvA( i , j ) + Xv( k , j ) *A( k , i ) ;
                            XvM( i , j ) = XvM( i , j ) + Xv( k , j ) *M( k , i ) ;
                            YphiA ( i , j ) = YphiA ( i , j ) + Yphi ( k , j ) *A( k , i ) ;
                            YphiM( i , j ) = YphiM ( i , j ) + Yphi ( k , j ) *M( k , i ) ;
                            YuA( i , j ) = YuA( i , j ) + Yu( k , j ) *A( k , i ) ;
                            YuM( i , j ) = YuM( i , j ) + Yu( k , j ) *M( k , i ) ;
                            YvL( i , j ) = YvL( i , j ) + Yv( k , j ) *L ( i , k ) ;
                            YvLT( i , j ) = YvLT( i , j ) + Yv( k , j ) *L ( k , i ) ;
                            XphiL ( i , j ) = XphiL ( i , j ) + Xphi ( k , j ) *L( i , k ) ;
                            XuA( i , j ) = XuA( i , j ) + Xu( k , j ) *A( k , i ) ;
                            XuM( i , j ) = XuM( i , j ) + Xu( k , j ) *M( k , i ) ;
                            XvL( i , j ) = XvL( i , j ) + Xv( k , j ) *L( i , k ) ;
                            XvLT( i , j ) = XvLT( i , j ) +Xv( k , j ) *L( k , i ) ;
                            YphiL ( i , j ) = YphiL ( i , j ) + Yphi ( k , j ) *L( i , k ) ;
                            YuL( i , j ) = YuL( i , j ) + Yu( k , j ) *L( i , k ) ;
                            YuLT( i , j ) = YuLT( i , j ) + Yu( k , j ) *L( k , i ) ;
                            YvA( i , j ) = YvA( i , j ) + Yv( k , j ) *A( k , i ) ;
                            YvM( i , j ) = YvM( i , j ) + Yv( k , j ) *M( k , i ) ;
                            cXphiA ( i , j ) = cXphiA ( i , j ) + Xphi ( k , j ) *A( i , k ) ;
                            cXphiLT ( i , j ) = cXphiLT ( i , j ) + Xphi ( k , j ) *L( k , i ) ;
                            cYphiLT ( i , j ) = cYphiLT ( i , j ) + Yphi ( k , j ) *L( k , i ) ;
                            cYphiA ( i , j ) = cYphiA ( i , j ) + Yphi ( k , j ) *A( i , k ) ;
                            cXuA( i , j ) = cXuA( i , j ) + Xu( k , j ) *A( k , i ) ;
                            cXvLT( i , j ) = cXvLT( i , j ) + Xv( k , j ) *L ( k , i ) ;
                            c2XphiLT ( i , j ) = c2XphiLT ( i , j ) + Xphi ( k , j ) *L ( k , i ) ;
                            cXuL ( i , j ) = cXuL ( i , j ) + Xu( k , j ) *L ( i , k ) ;
                            cXvM( i , j ) = cXvM( i , j ) + Xv( k , j ) *M( k , i ) ;
                            cXphiL ( i , j ) = cXphiL ( i , j ) + Xphi ( k , j ) *L ( i , k ) ;
                            cYuM( i , j ) = cYuM( i , j ) + Yu( k , j ) *M( k , i ) ;
                            cYvL ( i , j ) = cYvL ( i , j ) + Yv( k , j ) *L ( i , k ) ;
                            cYphiL ( i , j ) = cYphiL ( i , j ) + Yphi ( k , j ) *L ( i , k ) ;
                            cYuLT( i , j ) = cYuLT( i , j ) + Yu( k , j ) *L( k , i ) ;
                            cYvA( i , j ) = cYvA( i , j ) + Yv( k , j ) *A( k , i ) ;
                            c2YphiLT ( i , j ) = c2YphiLT ( i , j ) + Yphi ( k , j ) *L( k , i ) ;
                        end
                    end
                end
              
                for i = 1 :N*Kx-(Kx-1)
                    for j = 1 : Q
                        for k = 1 :N*Kx-(Kx-1)
                            Pphi(k,i) = Pphi(k,i) + (XuLT(k,j)*YuM(i,j) + XphiA(k,j)*YphiM(i,j) + XvM(k,j)*YvLT(i,j) + XphiM(k,j)*YphiA(i,j) + ep*b(1)*cXuA(k,j)*cYuM(i,j) + ep*b(1)*cXvLT(k,j)*cYvL(i,j) + b(1)*b(1)*XphiA(k,j)*YphiM(i,j) + b(1)*b(2)*c2XphiLT(k,j)*cYphiL(i,j) + ep*b(2)*cXuL(k,j)*cYuLT(i,j) + ep*b(2)*cXvM(k,j)*cYvA(i,j) + b(2)*b(1)*cXphiL(k,j)*c2YphiLT(i,j) + b(2)*b(2)*XphiM(k,j)*YphiA(i,j));
                            Pu(k,i) = Pu(k,i) +(ep*b(1)*cXphiA(k,j)*YphiM(i,j) + ep*b(2)*cXphiLT(k,j)*cYphiL(i,j) + ep*ep*XuA(k,j)*YuM(i,j) + ep*ep*XvLT(k,j)*YvL(i,j) + XuM(k,j)*YuM(i,j) + XphiL(k,j)*YphiM(i,j) - XvL(k,j)*YvLT(i,j) + XuM(k,j)*YuA(i,j)) ;
                            Pv(k,i) = Pv(k,i) + (ep*b(1)*cXphiL(k,j)*cYphiLT(i,j) + ep*b(2)*XphiM(k,j)*cYphiA(i,j) + ep*ep*XvM(k,j)*YvA(i,j) + ep*ep*XuL(k,j)*YuLT(i,j) + XvM(k,j)*YvM(i,j) + XphiM(k,j)*YphiL(i,j) - XuLT(k,j)*YuL(i,j) + XvA(k,j)*YvM(i,j));
                        end
                    end
                end
            
                %Unknown functions for enrichment
                oldR = zeros( 3*((N)*Kx-(Kx-1)),1) ;
                oldS = zeros( 3*((N)*Kx-(Kx-1)),1) ;
                
                Rphi = rand(N*Kx-(Kx-1) ,1) ;
                Ru = rand(N*Kx-(Kx-1), 1 ) ;
                Rv = rand(N*Kx-(Kx-1), 1 ) ;
                
                counter =0; 
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
                    
                    M2_11 = (Rphi'*A*Rphi)*M + (Rphi'*M*Rphi)*A + b(1)*b(1)*(Rphi'*A*Rphi)*M + b(1)*b(2)*(Rphi'*L'*Rphi)*L + b(2)*b(1)*(Rphi'*L*Rphi)*L'+ b(2)*b(2)*(Rphi'*M*Rphi)*A ;
                    M2_12 = (Rphi'*L'*Ru)*M +ep*b(1)*(Rphi'*A*Ru)*M+ ep*b(2)*(Rphi'*L*Ru)*L';
                    M2_13 = (Rphi'*M*Rv)*L' + ep*b(1)*(Rphi'*L'*Rv)*L + ep*b(2)*(Rphi'*M*Rv)*A;
                    M2_21 = (Ru'*L*Rphi)*M + ep*b(1)*(Ru'*A*Rphi)*M + ep*b(2)*(Ru'*L'*Rphi)*L;
                    M2_22 = ep*ep*(Ru'*A*Ru)*M + (Ru'*M*Ru)*A + (Ru'*M*Ru)*M;
                    M2_23 = ep*ep*(Ru'*L'*Rv)*L - (Ru'*L*Rv)*L';
                    M2_31 = (Rv'*M*Rphi)*L + ep*b(1)*(Rv'*L*Rphi)*L' + ep*b(2)*(Rv'*M*Rphi)*A;
                    M2_32 = ep*ep*(Rv'*L*Ru)*L' - (Rv'*L'*Ru)*L;
                    M2_33 = (Rv'*A*Rv)*M + ep*ep*(Rv'*M*Rv)*A + (Rv'*M*Rv)*M; 
        
                    Matrix2 = [M2_11, M2_12, M2_13; M2_21, M2_22, M2_23; M2_31, M2_32, M2_33] ;
                    V2_1 =  - Pphi'*Rphi +b(1)*B'*Rphi + b(2)*C'*Rphi;
                    V2_2 = ep*B'*Ru-Pu'*Ru;
                    V2_3 = ep*C'*Rv-Pv'*Rv ;
                    Vector2 = [ V2_1 ; V2_2 ; V2_3 ] ;
        
                    % trial taking out BCs
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
                    %Solve for x
                    M1_11 = (Sphi'*M*Sphi)*A+(Sphi'*A*Sphi)*M + b(1)*b(1)*(Sphi'*M*Sphi)*A + b(1)*b(2)*(Sphi'*L*Sphi)*L' + b(2)*b(1)*(Sphi'*L'*Sphi)*L + b(2)*b(2)*(Sphi'*A*Sphi)*M;
                    M1_12 = (Sphi'*M*Su)*L' + ep*b(1)*(Sphi'*M*Su)*A + ep*b(2)*(Sphi'*L'*Su)*L;
                    M1_13 = (Sphi'*L'*Sv)*M + ep*b(1)*(Sphi'*L*Sv)*L'+ep*b(2)*(Sphi'*A*Sv)*M;
                    M1_21 = (Su'*M*Sphi)*L + ep*b(1)*(Su'*M*Sphi)*A + ep*b(2)*(Su'*L*Sphi)*L';
                    M1_22 = ep*ep*(Su'*M*Su)*A + (Su'*A*Su)*M + (Su'*M*Su)*M;
                    M1_23 = ep*ep*(Su'*L*Sv)*L' - (Su'*L'*Sv)*L;
                    M1_31 = (Sv'*L*Sphi)*M + ep*b(1)*(Sv'*L'*Sphi)*L + ep*b(2)*(Sv'*A*Sphi)*M;
                    M1_32 = ep*ep*(Sv'*L'*Su)*L-(Sv'*L*Su)*L';
                    M1_33 = ep*ep*(Sv'*A*Sv)*M+(Sv'*M*Sv)*M+(Sv'*M*Sv)*A;
                    Matrix1 = [M1_11,M1_12,M1_13;M1_21,M1_22,M1_23;M1_31,M1_32,M1_33];
                    
                    V1_1 = -Pphi*Sphi + b(1)*B*Sphi + b(2)*C*Sphi;
                    V1_2 = -Pu*Su + ep*B*Su;
                    V1_3 = -Pv*Sv + ep*C*Sv;
                     
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
                    
                    counter=counter+1;
                    if counter == max_iterations
                        break 
                    end
                
                end

                cumcounter = [cumcounter + counter]; 
                cumcounter_values = [cumcounter_values, cumcounter];
                each_count_value = [each_count_value, counter]; 
                endTime = cputime;       
        
                % Calculate execution time for this enrichment (the only enrichment!)
                enrichmentTime = endTime - startTime;
                cumTime(bi) = cumTime(bi) + enrichmentTime;
           
                %setting XQ YQ from the calculated values
                xphi = Rphi ;
                xu= Ru;
                xv = Rv ;
                yphi = Sphi ;
                yu = Su ;
                yv = Sv ;
        
                %adding the known boundary v al u e s
                Xphi = [ Xphi , xphi ] ;
                Xu = [ Xu, xu ] ;
                Xv = [ Xv, xv ] ;
                Yphi = [ Yphi , yphi ] ;
                Yu = [ Yu , [ yu ] ] ;
                Yv = [ Yv, yv ] ;
                Q= Q+1;
                executionTimes(Q,bi) = cumTime(bi);
    
                %calculating the approximations
                approxphi = zeros(N*Kx-(Kx-1) ) ;
                approxu = zeros(N*Kx-(Kx-1) ) ;
                approxv = zeros(N*Kx-(Kx-1) ) ;
                for k = 1 :N*Kx-(Kx-1) 
                    for l = 1 :N*Kx-(Kx-1) 
                        for j = 1 : Q
                            approxphi(k,l) = approxphi(k,l) + Yphi(l,j)*Xphi(k,j);
                            approxu(k,l) = approxu(k,l) + Yu(l,j)*Xu(k,j) ;
                            approxv(k,l) = approxv(k,l) + Yv(l,j)*Xv(k,j) ;
                        end
                    end
                end

                %working out the error
                E = 0 ;
                Emax = 0 ;
                PHIerror= zeros(N*Kx-(Kx-1) ) ;
                PHImesh= zeros(N*Kx-(Kx-1) ) ;
                Umesh = zeros(N*Kx-(Kx-1) ) ;
                Vmesh = zeros(N*Kx-(Kx-1) ) ;
                Uerror = zeros(N*Kx-(Kx-1) ) ;
                Verror = zeros(N*Kx-(Kx-1) ) ;
         
                for k = 1 :N*Kx-(Kx-1) 
                    for l = 1 :N*Kx-(Kx-1)  
                        Uerror(k,l) = approxu(k,l) - U(xdomain(k),ydomain(l));
                        Verror(k,l) = approxv(k,l) - V(xdomain(k),ydomain(l));
                        PHIerror(k,l) = approxphi(k,l) - phi(xdomain(k),ydomain(l));
                        PHImesh(k,l) = phi(xdomain(k),ydomain(l));
                        Umesh(k,l) = U(xdomain(k),ydomain(l));
                        Vmesh(k,l) = V(xdomain(k),ydomain(l)); 
                    end
                end
            end
         
            L2errorphi=0;
            L2erroru=0;
            L2errorv=0;
            for i=1:Kx
                for j=1:Ky
                    for num=1:N1
                        for m=1:N1
                            L2errorphi=L2errorphi+(ax(i+1)-ax(i))*(cy(j+1)-cy(j))/(2^(Kx))*w(m)*w(num)*(approxphi(m,num)-PHImesh(m,num))^2 ;
                            L2erroru=L2erroru+(ax(i+1)-ax(i))*(cy(j+1)-cy(j))/(2^(Kx))*w(m)*w(num)*(approxu(m,num)-Umesh(m,num))^2 ;
                            L2errorv=L2errorv+(ax(i+1)-ax(i))*(cy(j+1)-cy(j))/(2^(Kx))*w(m)*w(num)*(approxv(m,num)-Vmesh(m,num))^2 ;
                        end
                    end
                end
            end
            L2errorphi= L2errorphi^(1/2);
            L2erroru= L2erroru^(1/2);
            L2errorv= L2errorv^(1/2);
            
            NBL2phi(Num,bi,ei)=L2errorphi;
            NBL2u(Num,bi,ei)=L2erroru;
            NBL2v(Num,bi,ei)=L2errorv;

        end
        N_values=Nvalues';
        allb_cumcounter(bi,ei) = cumcounter_values;
    
    end

    figure('Units', 'normalized', 'Position', [0.2 0.2 0.5 0.5]); 
    set(gcf, 'Color', 'white');
    handles = semilogy(N_values,NBL2phi(:,:,ei), 'Linewidth', 2);
    legend('$\mathbf{b}=(0,1)^T$','$\mathbf{b}=(1,0)^T$','$\mathbf{b}=(0.01,1)^T$', '$\mathbf{b}=(1,0.01)^T$', '$\mathbf{b}=(0.1,1)^T$', '$\mathbf{b}=(1,0.1)^T$', 'Interpreter', 'latex','FontSize', 14 )
    grid on
    xlabel('$N$','interpreter','latex','fontsize',20)
    xlim([N_values(1) N_values(end)])
    xticks(N_values)
    ylabel('$L^{2}$ Error in $\phi$','interpreter','latex','fontsize',20)
    title(sprintf('$\\bar{\\epsilon} = %.3f$', ep), ...
         'Interpreter', 'latex', 'fontsize', 14);
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
fprintf('L2 errors for φ for Example 6 using the LS SEM PGD, after one enrichment with N=24:')
fprintf('\n');
disp(T)


rowLabels = { ...
    'b=(0,1)^T', ...
    'b=(1,0)^T', ...
    'b=(0.01,1)^T', ...
    'b=(1,0.01)^T', ...
    'b=(0.1,1)^T', ...
    'b=(1,0.1)^T'};

T = array2table(allb_cumcounter, ...
    'VariableNames', {'epsilon bar = 1', 'epsilon bar = 0.001'}, ...  
    'RowNames', rowLabels);
fprintf('\n');
fprintf('Number of ADFPA iterations required at the first enrichment stage for Example 6 using the LS SEM PGD, for N=24 and epsilon = %6.2e\n', epsilon )
fprintf('\n');
disp(T)

fprintf('Time taken for the last example (in seconds):')
disp(enrichmentTime)
