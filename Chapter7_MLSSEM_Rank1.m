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
clear all
close all
 
syms x y

% Number of enrichments
no_enrichment = 1; 
    
% True solution
phi=@(x,y)  sin(pi*x)*sin(pi*y);
U=@(x,y) - pi*cos(pi*x).*sin(pi*y); 
V=@(x,y) - pi*sin(pi*x).*cos(pi*y);
 
% N GLL points 
Nvalues = [8,12,16,20];
  
% Max number of iterations  
max_iterations = 100; 
   
% Tolerance
epsilon = 1e-8;
  
for number =1:length(Nvalues);
    
    clearvars -except number Nvalues epsilon max_iterations no_enrichment U V phi L2phi L2U L2V J1_all J2_all J3_all J4_all eq7_16_all eq7_17_all eq7_18_all eq7_19_all
    Kx=3;  
    Ky=Kx; 
    cumcounter =0;
    cumcounter_values =[]; 
    each_count_value=[];
         
    aa=-1;
    bb=1;
    cc=-1;
    dd=1;
    ax=linspace(aa,bb,Kx+1);
    cy=linspace(cc,dd,Ky+1); 
    
    N = Nvalues(number);
    N1=N;
    NN=N-1;
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
 
%============= Differentiation matrix and weights ==================

    X=repmat(xx,1,N1); %Creates a square matrix with column entries x
    Xdiff=X-X'+eye(N1); %eye is for Kronecker Delta
    
    LL=repmat(P(:,N1),1,N1); %Replicates Legendre Vandemonde Matrix
    LL(1:(N1+1):N1*N1)=1;
    D=(LL./(Xdiff.*LL'));
    D(1:(N1+1):N1*N1)=0;
    D(1)=-(N1*NN)/4;
    D(N1*N1)=(N1*NN)/4; %Gives differentiation matrix
    w=2./(NN*N1*(LL(:,N1).^2)); %Computation of weights

%=======================================================================
% Set up matrices
%=======================================================================   
    
    t_tilde=node; 
    M = N; 
    for i=1:N-1 
        for k=1:M
            x_tilde(i,k) = ((node(i+1) - node(i))*t_tilde(k) + (node(i+1) + node(i)))/2;
        end
    end
    t= x_tilde;
  
    for i =1:N-1 
        for j=1:N
            prodk = 1; 
            sumk = 0; 
            for k = 1:N
                prodk = 1; 
                for l =1:N
                if l ~= j
                    numerator = t(i, k) - node(l);
                    denominator = node(j) - node(l); 
                    prodk  =  prodk .* (numerator  ./ denominator)  ;
                end
                end  
                h(j,k)  = prodk ;
                sumk = sumk + w(k)*h(j,k) ; 
            end 
            T(i,j) = 0.5*(node(i+1) - node(i))*sumk  ;
           
        end 
    end 
    F = zeros(N,N-1);
    for m=1:N
        for j= 1:N-1 
            sumk =0;
            for k = 1:j 
                sumk = sumk + D(m,k);
            end
            F(m,j) = sumk;
        end
    end
    F=-F;

    G = zeros(N);
    G(:, 1) = -F(:,1);
    G(:, end) = F(:,end);
    G(:, 2:end-1) = D(:, 2:end-1);

    zeroN = zeros(N-1,1);
    zeroN1 = zeros(N,1);
    zerosN1=zeroN1;
    for i=1:N-1
        for j=1:N
        sumM =0;
        for m = 1:N
            sumM = sumM + w(m)*F(m,i)*G(m,j);
        end
        wFG(i,j)= sumM;
        end
    end
    H=wFG;
    HT=H'; 
    
    for i=1:N-1
        for j=1:N-1
        sumM =0;
        for m = 1:N
            sumM = sumM + w(m)*F(m,i)*F(m,j);
        end
        Ft(i,j)= sumM;
        end
    end
    
    Gt = zeros(N,N);
    for i=1:N
        for j=1:N
        sumM =0;
        for m = 1:N
            sumM = sumM + w(m)*G(m,i)*G(m,j);
        end
        Gt(i,j)= sumM;
        end
    end
    
    M=zeros(N,N);
    for i =1:N
        M(i,i) = w(i);
    end
    M2 = M(2:end-1,2:end-1);
    
    Z = zeros(N,N);
    M2 = M(2:end-1,2:end-1);
    M4 = M(:,2:end-1);
    Gt2 = Gt(2:end-1,2:end-1);
    Gt4 = Gt(:,2:end-1);
    HT2 = HT(2:end-1,:);
    H2 = H(:,2:end-1);
    T2 = T(:,2:end-1);
     
    TFtT = T'*Ft*T;
    TH = T'*H;
    HTT = HT*T;
    
    HTT2 = HTT(2:end-1,:);
    TH2 = TH(:,2:end-1);
    HTT3 = HTT(:,2:end-1);
    TH3 = TH(2:end-1,:);
    TFtT2=TFtT(2:end-1, 2:end-1);
     
    Mx=M;
    My=M; 
    THx=TH;
    THy=TH;
    TFtTx=TFtT;
    TFtTy=TFtT;
    Gtx=Gt;
    Gty=Gt;
    %  
    dx=zeros(Kx,1)';
    dy=zeros(Ky,1)'; 
    localFx=zeros(1,N);
    localFy=zeros(1,N);
    % % Setting up the global matrices
    GlobalTHx=zeros((N)*Kx-(Kx-1),(N)*Kx-(Kx-1));
    GlobalMx=zeros((N)*Kx-(Kx-1),(N)*Kx-(Kx-1));
    GlobalTFtTx=zeros((N)*Kx-(Kx-1),(N)*Kx-(Kx-1));
    GlobalGtx=zeros((N)*Kx-(Kx-1),(N)*Kx-(Kx-1));
    GlobalV1bx=zeros((N)*Kx-(Kx-1),(N)*Kx-(Kx-1));
    GlobalV1cx=zeros((N)*Kx-(Kx-1),(N)*Kx-(Kx-1));
    GlobalTHy=zeros((N)*Ky-(Ky-1),(N)*Ky-(Ky-1));
    GlobalMy=zeros((N)*Ky-(Ky-1),(N)*Ky-(Ky-1));
    GlobalTFtTy=zeros((N)*Ky-(Ky-1),(N)*Ky-(Ky-1));
    GlobalGty=zeros((N)*Ky-(Ky-1),(N)*Ky-(Ky-1));
    GlobalV2by=zeros((N)*Ky-(Ky-1),(N)*Ky-(Ky-1));
    GlobalV2cy=zeros((N)*Ky-(Ky-1),(N)*Ky-(Ky-1));
    
    GlobalT=zeros((N)*Kx-(Kx-1)-1,(N)*Kx-(Kx-1));
    GlobalTa=zeros((N)*Kx-(Kx-1)-1,(N)*Kx-(Kx-1));
    GlobalG=zeros((N)*Kx-(Kx-1),(N)*Kx-(Kx-1));
    GlobalF=zeros((N)*Kx-(Kx-1),(N)*Kx-(Kx-1)-1);

% %========== Local and Global matrix computations ========================
% %========================================================================

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
        localTHx=2/(dx(i))*THx;
        %Global Stiffness matrix for x
        GlobalTHx([(i-1)*N-(i-2):(i)*N-(i-1)],[(i-1)*N-(i-2):(i)*N-(i-1)]) = GlobalTHx([(i-1)*N-(i-2):(i)*N-(i-1)],[(i-1)*N-(i-2):(i)*N-(i-1)]) + localTHx;
        
        localTFtTx=2/(dx(i))*TFtTx;
        GlobalTFtTx([(i-1)*N-(i-2):(i)*N-(i-1)],[(i-1)*N-(i-2):(i)*N-(i-1)]) = GlobalTFtTx([(i-1)*N-(i-2):(i)*N-(i-1)],[(i-1)*N-(i-2):(i)*N-(i-1)]) + localTFtTx;
        
        localGtx=2/(dx(i))*Gtx;
        GlobalGtx([(i-1)*N-(i-2):(i)*N-(i-1)],[(i-1)*N-(i-2):(i)*N-(i-1)]) = GlobalGtx([(i-1)*N-(i-2):(i)*N-(i-1)],[(i-1)*N-(i-2):(i)*N-(i-1)]) + localGtx;
        
        localMx=(dx(i))/2*Mx;
        %Global Mass matrix for x
        GlobalMx([(i-1)*N-(i-2):(i)*N-(i-1)],[(i-1)*N-(i-2):(i)*N-(i-1)])=GlobalMx([(i-1)*N-(i-2):(i)*N-(i-1)],[(i-1)*N-(i-2):(i)*N-(i-1)]) + localMx;    
        
        localF=2/(dx(i))*F;
        GlobalF([(i-1)*N-(i-2):(i)*N-(i-1)],[(i-1)*N-(i-2):(i)*N-(i-1)-1]) = GlobalF([(i-1)*N-(i-2):(i)*N-(i-1)],[(i-1)*N-(i-2):(i)*N-(i-1)-1]) + localF;
        
        localG=2/(dx(i))*G;
        GlobalG([(i-1)*N-(i-2):(i)*N-(i-1)],[(i-1)*N-(i-2):(i)*N-(i-1)]) = GlobalG([(i-1)*N-(i-2):(i)*N-(i-1)],[(i-1)*N-(i-2):(i)*N-(i-1)]) + localGtx;
 
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
        localTHy=2/(dy(i))*THy;
        %Global Stiffness matriy for y
        GlobalTHy([(i-1)*N-(i-2):(i)*N-(i-1)],[(i-1)*N-(i-2):(i)*N-(i-1)]) = GlobalTHy([(i-1)*N-(i-2):(i)*N-(i-1)],[(i-1)*N-(i-2):(i)*N-(i-1)]) + localTHy;
        
        %Local Stiffness matrix for each element in y
        localTFtTy=2/(dy(i))*TFtTy;
        %Global Stiffness matriy for y
        GlobalTFtTy([(i-1)*N-(i-2):(i)*N-(i-1)],[(i-1)*N-(i-2):(i)*N-(i-1)]) = GlobalTFtTy([(i-1)*N-(i-2):(i)*N-(i-1)],[(i-1)*N-(i-2):(i)*N-(i-1)]) + localTFtTy;
        
        %Local Stiffness matrix for each element in y
        localGty=2/(dy(i))*Gty;
        %Global Stiffness matriy for y
        GlobalGty([(i-1)*N-(i-2):(i)*N-(i-1)],[(i-1)*N-(i-2):(i)*N-(i-1)]) = GlobalGty([(i-1)*N-(i-2):(i)*N-(i-1)],[(i-1)*N-(i-2):(i)*N-(i-1)]) + localGty;
        
        
        %Local Mass matrix for each element in y
        localMy=(dy(i))/2*My;
        %Global Mass matrix for y
        GlobalMy([(i-1)*N-(i-2):(i)*N-(i-1)],[(i-1)*N-(i-2):(i)*N-(i-1)]) = GlobalMy([(i-1)*N-(i-2):(i)*N-(i-1)],[(i-1)*N-(i-2):(i)*N-(i-1)]) + localMy;    
 
    end
 
    GlobalV1bx= zeros((N)*Kx-(Kx-1));
    GlobalV1cx= zeros((N)*Kx-(Kx-1));
    GlobalV2by= zeros((N)*Ky-(Ky-1));
    GlobalV2cy= zeros((N)*Ky-(Ky-1));
    V1bx=zeros(N1);
    V1cx=zeros(N1);
    V2by=zeros(N1);
    V2cy=zeros(N1);
    
    localDiffMatrix=zeros(N,N);
    GlobalDiff=zeros((N)*Ky-(Ky-1),(N)*Ky-(Ky-1));
  
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
            for ii = 1:N
                for jj = 1:N
                    Fxy(ii,jj) = 2 * pi^2 * sin(pi*physx(ii)) * sin(pi*physy(jj))  ;  
                end
            end 
         
            TT=T; 
            TT= 2/dx(i) *TT; 

            for ii=1:N-1
                 for jj=1:N-1
                     suma =0;
                     for aa =1:N
                         sumb = 0;
                         for bb=1:N
                             sumb = sumb + TT(jj,bb)*Fxy(aa,bb);
                         end
                         suma = suma + sumb*TT(ii,aa);
                     end
                     fxy(ii,jj) = suma;
                 end 
            end
    
            fxy = dx(i)*dy(k)/4 * fxy;
                  
            V2by =  T'*Ft*fxy*HT' ; 
            V2cy = HT*fxy*Ft'*T ;
    
            V1bx = T'*Ft*fxy*HT' ;  
            V1cx = HT*fxy*Ft'*T ;  
         
            localV1bx = V1bx  ;  
            localV1cx = V1cx  ;  
                                     
            localV2by = V2by  ;
            localV2cy = V2cy  ;
        
            GlobalV2by([(i-1)*N-(i-2):(i)*N-(i-1)],[(k-1)*N-(k-2):(k)*N-(k-1)])=GlobalV2by([(i-1)*N-(i-2):(i)*N-(i-1)],[(k-1)*N-(k-2):(k)*N-(k-1)]) + localV2by;    
            
            GlobalV1bx([(i-1)*N-(i-2):(i)*N-(i-1)],[(k-1)*N-(k-2):(k)*N-(k-1)])=GlobalV1bx([(i-1)*N-(i-2):(i)*N-(i-1)],[(k-1)*N-(k-2):(k)*N-(k-1)]) + localV1bx;    
            
            GlobalV2cy([(i-1)*N-(i-2):(i)*N-(i-1)],[(k-1)*N-(k-2):(k)*N-(k-1)])=GlobalV2cy([(i-1)*N-(i-2):(i)*N-(i-1)],[(k-1)*N-(k-2):(k)*N-(k-1)]) + localV2cy;    
            
            GlobalV1cx([(i-1)*N-(i-2):(i)*N-(i-1)],[(k-1)*N-(k-2):(k)*N-(k-1)])=GlobalV1cx([(i-1)*N-(i-2):(i)*N-(i-1)],[(k-1)*N-(k-2):(k)*N-(k-1)]) + localV1cx;    
     
            localT = T; 
            GlobalT([(i-1)*N-(i-2):(i)*N-(i-1)-1],[(k-1)*N-(k-2):(k)*N-(k-1)]) = GlobalT([(i-1)*N-(i-2):(i)*N-(i-1)-1],[(k-1)*N-(k-2):(k)*N-(k-1)]) + localT;

        end
        localDiffMatrix = 2/(dx(i)) * D; 
        GlobalDiff([(i-1)*N-(i-2):(i)*N-(i-1)],[(i-1)*N-(i-2):(i)*N-(i-1)]) = GlobalDiff([(i-1)*N-(i-2):(i)*N-(i-1)],[(i-1)*N-(i-2):(i)*N-(i-1)]) + localDiffMatrix;
           
        localT = T; 
        GlobalTa([(i-1)*N-(i-2):(i)*N-(i-1)-1],[(i-1)*N-(i-2):(i)*N-(i-1)]) = GlobalTa([(i-1)*N-(i-2):(i)*N-(i-1)-1],[(i-1)*N-(i-2):(i)*N-(i-1)]) + localT;

    end

    V1bx = GlobalV1bx;
    V1cx = GlobalV1cx; 
    
    V2by = GlobalV2by; 
    V2cy = GlobalV2cy; 

% %========================= PGD begins ====================================
% %=========================================================================
 
    M=GlobalMx;
    TH=GlobalTHx;
    HTT=TH';
    Gt=GlobalGtx;
    TFtT=GlobalTFtTx;
     
    %%
    Q=0; 
    %Start approximations
    Xphi = zeros(N*Kx-(Kx-1),Q ) ;
    Xu = zeros(N*Kx-(Kx-1),Q ) ;
    Xv = zeros(N*Kx-(Kx-1),Q ) ;
    Yphi = zeros(N*Kx-(Kx-1),Q ) ;
    Yu = zeros(N*Kx-(Kx-1),Q ) ;
    Yv = zeros(N*Kx-(Kx-1),Q ) ;
 
    startTime = cputime;
    for term = 1:no_enrichment

        %Initialize the Matrices containing the known Q−1 enrichments
        Pphi_1 = zeros(N*Kx-(Kx-1)) ;
        Pu_1 = zeros(N*Kx-(Kx-1)) ;
        Pv_1 = zeros(N*Kx-(Kx-1)) ;
        Pphi_2 = zeros(N*Kx-(Kx-1)) ;
        Pu_2 = zeros(N*Kx-(Kx-1)) ;
        Pv_2 = zeros(N*Kx-(Kx-1)) ;
        
        XphiM = zeros(N*Kx-(Kx-1),Q ) ;
        XphiGt = zeros(N*Kx-(Kx-1),Q ) ;
        XvHTT = zeros(N*Kx-(Kx-1),Q ) ;
        XuM = zeros(N*Kx-(Kx-1),Q ) ;
        XphiTH = zeros(N*Kx-(Kx-1),Q ) ;
        XvTFtT = zeros(N*Kx-(Kx-1),Q ) ;
        XvGt = zeros(N*Kx-(Kx-1),Q ) ;
        XuTH = zeros(N*Kx-(Kx-1),Q ) ;
        XvTH = zeros(N*Kx-(Kx-1),Q ) ;
        XuGt = zeros(N*Kx-(Kx-1),Q ) ;
        XuTFtT = zeros(N*Kx-(Kx-1),Q ) ;
        YphiGt = zeros(N*Kx-(Kx-1),Q ) ;
        YphiM = zeros(N*Kx-(Kx-1),Q ) ; 
        YvM = zeros(N*Kx-(Kx-1),Q ) ;
        YvGt = zeros(N*Kx-(Kx-1),Q ) ;
        YvTFtT = zeros(N*Kx-(Kx-1),Q ) ;
        YuHTT = zeros(N*Kx-(Kx-1),Q ) ;
        YuTH = zeros(N*Kx-(Kx-1),Q ) ;
        YphiTH = zeros(N*Kx-(Kx-1),Q ) ;
        YvTH = zeros(N*Kx-(Kx-1),Q ) ;
        YvHTT = zeros(N*Kx-(Kx-1),Q ) ;
        YuTFtT = zeros(N*Kx-(Kx-1),Q ) ;
        YuGt = zeros(N*Kx-(Kx-1),Q ) ;
        
        XuHTT = zeros(N*Kx-(Kx-1),Q ) ; 
        YuM = zeros(N*Kx-(Kx-1),Q ) ; 
        XvM = zeros(N*Kx-(Kx-1),Q ) ; 
        YvHT = zeros(N*Kx-(Kx-1),Q ) ;  
        YuHT = zeros(N*Kx-(Kx-1),Q ) ; 
        
        for i = 1 : N*Kx-(Kx-1)
            for j = 1 : Q
                for k = 1 : N*Kx-(Kx-1)
                    XphiM ( i , j ) = XphiM ( i , j ) + Xphi ( k , j ) *M( i , k ) ;
                    XphiGt( i , j ) = XphiGt( i , j ) + Xphi( k , j ) *Gt( i , k ) ;
                    XvHTT( i , j ) = XvHTT( i , j ) + Xv( k , j ) *HTT( i , k ) ; 
                    XuM( i , j ) = XuM( i , j ) + Xu ( k , j ) *M( i , k);
                    XphiTH ( i , j ) = XphiTH ( i , j ) + Xphi ( k , j ) *TH( i , k ) ;  
                    XvTFtT ( i , j ) = XvTFtT ( i , j ) + Xv ( k , j ) *TFtT( i , k ) ;  
                    XvGt ( i , j ) = XvGt ( i , j ) + Xv ( k , j ) *Gt( i , k ) ;
                    XuTH ( i , j ) = XuTH ( i , j ) + Xu ( k , j ) *TH( i , k ) ;
                    XvTH ( i , j ) = XvTH ( i , j ) + Xv ( k , j ) *TH( i , k ) ;
                    XuGt ( i , j ) = XuGt ( i , j ) + Xu ( k , j ) *Gt( i , k ) ;
                    XuTFtT ( i , j ) = XuTFtT ( i , j ) + Xu ( k , j ) *TFtT( i , k ) ;
                    
                    YphiGt ( i , j ) = YphiGt ( i , j )  + Yphi ( k , j ) *Gt( i , k ) ;
                    YphiM ( i , j ) = YphiM ( i , j )  + Yphi ( k , j ) *M( i , k ) ;
                    YuHTT ( i , j ) = YuHTT ( i , j )  + Yu ( k , j ) *HTT( i , k ) ;
                    YvM ( i , j ) = YvM ( i , j )  + Yv ( k , j ) *M( i , k ) ;
                    YvGt ( i , j ) = YvGt ( i , j )  + Yv ( k , j ) *Gt( i , k ) ;
                    YvTFtT ( i , j ) = YvTFtT ( i , j )   + Yv ( k , j ) *TFtT( i , k ) ;
                    YuTH ( i , j ) = YuTH ( i , j )  + Yu ( k , j ) *TH( i , k ) ;
                    YphiTH ( i , j ) = YphiTH ( i , j )  + Yphi ( k , j ) *TH( i , k ) ;
                    YvTH ( i , j ) = YvTH ( i , j )  + Yv ( k , j ) *TH( i , k ) ;
                    YvHTT ( i , j ) = YvHTT ( i , j )  + Yv ( k , j ) *HTT( i , k ) ;
                    YuTFtT ( i , j ) = YuTFtT ( i , j ) + Yu ( k , j ) *TFtT( i , k ) ;
                    YuGt ( i , j ) = YuGt ( i , j ) + Yu ( k , j ) *Gt( i , k ) ;
                    
                    XuHTT ( i , j ) =  XuHTT ( i , j ) + Xu ( k , j ) *HTT( i , k ) ;
                    YuM ( i , j ) = YuM ( i , j ) + Yu ( k , j ) *M( i , k ) ; 
                    XvM ( i , j ) = XvM ( i , j ) + Xv ( k , j ) *M( i , k ) ; 
                 end
            end
        end
         
        for i =1: N*Kx-(Kx-1)
            for k =1: N*Kx-(Kx-1)
                for j=1:Q
                    %system 1
                    Pphi_1(k,i) = Pphi_1(k,i) +  XphiM(k,j)*YphiGt(i,j) + XphiGt(k,j)*YphiM(i,j) + XvM(k,j)*YvHTT(i,j) + XuHTT(k,j)*YuM(i,j); 
                    Pv_1(k,i) = Pv_1(k,i) + XphiM(k,j) * YphiTH(i,j) + XvM(k,j) * YvTFtT(i,j) + XvTFtT(k,j) * YvGt(i,j) + XvGt(k,j) * YvTFtT(i,j) - XuHTT(k,j) * YuTH(i,j) + XuTH(k,j) * YuHTT(i,j) ;
                    Pu_1(k,i) = Pu_1(k,i) + XphiTH(k,j) * YphiM(i,j)- XvTH(k,j) * YvHTT(i,j)  + XvHTT(k,j) * YvTH(i,j) + XuGt(k,j) * YuTFtT(i,j) + XuTFtT(k,j) * YuM(i,j) + XuTFtT(k,j) * YuGt(i,j);
                    
                    %system 2
                    Pphi_2(k,i) = Pphi_2(k,i) +  XphiM(k,j) * YphiGt(i,j) + XphiGt(k,j) * YphiM(i,j) + XuHTT(k,j)*YuM(i,j) + XvM(k,j)*YvHTT(i,j); 
                    Pu_2(k,i) = Pu_2(k,i) + XphiTH(k,j) * YphiM(i,j) + XuTFtT(k,j) * YuM(i,j) + XuTFtT(k,j) * YuGt(i,j) + XuGt(k,j) * YuTFtT(i,j) - XvTH(k,j) * YvHTT(i,j) + XvHTT(k,j) * YvTH(i,j) ;
                    Pv_2(k,i) = Pv_2(k,i) + XphiM(k,j) * YphiTH(i,j)- XuHTT(k,j) * YuTH(i,j)  + XuTH(k,j) * YuHTT(i,j) + XvGt(k,j) * YvTFtT(i,j) + XvM(k,j) * YvTFtT(i,j) + XvTFtT(k,j) * YvGt(i,j);
                
                end
            end
        end
 
        %Unknown functions for enrichment
        oldR = zeros( 3*((N)*Kx-(Kx-1)),1) ;
        oldS = zeros( 3*((N)*Kx-(Kx-1)),1) ;
        
        Sphi = rand(N*Kx-(Kx-1) ,1) ;
        Su = rand(N*Kx-(Kx-1), 1 ) ;
        Sv = rand(N*Kx-(Kx-1), 1 ) ;
        if term >1
            Rphi = Xphi(:,term-1);
            Ru = Xu(:,term-1);
            Rv = Xv(:,term-1);
           
        else
        end
        counter = 0;
        if term >1
            Sphi = Yphi(:,term-1);
            Su = Yu(:,term-1);
            Sv = Yv(:,term-1);
        else
        end
        
        Rphi(1) = 0;
        Rphi(end)=0;
        Rv(1) = 0;
        Rv(end)=0;
        Sphi(1) = 0;
        Sphi(end)=0;
        Su(1) = 0;
        Su(end)=0;
        R1 = randn(3*((N)*Kx-(Kx-1)),1);  
        S1 = randn(3*((N)*Kx-(Kx-1)),1);
          
        while max(max( abs (R1*S1' - oldR* oldS') ) ) >epsilon
            oldR = R1 ;
            oldS = S1 ;
             
            M1_11 = (Sphi'*M*Sphi)*Gt + (Sphi'*Gt*Sphi)*M ;  
            M1_12 = (Sphi'*HTT*Sv)*M ; 
            M1_13 = (Sphi'*M*Su)*HTT; 
            M1_21 = (Sv'*TH*Sphi)*M ; 
            M1_22 = (Sv'*TFtT*Sv)*M + (Sv'*TFtT*Sv)*Gt + (Sv'*Gt*Sv)*TFtT; 
            M1_23 = -(Sv'*TH*Su)*HTT + (Sv'*HTT*Su)*TH; 
            M1_31 = (Su'*M*Sphi)*TH;  
            M1_32 = -(Su'*HTT*Sv)*TH + (Su'*TH*Sv)*HTT; 
            M1_33 = (Su'*Gt*Su)*TFtT + (Su'*M*Su)*TFtT + (Su'*TFtT*Su)*Gt; 
            
            Matrix1 = [  M1_11, M1_12, M1_13;       
                         M1_21, M1_22, M1_23  ;
                         M1_31, M1_32, M1_33  ] ;   
                                                     
            v1_1 = - Pphi_1*Sphi;    
            v1_2 = V1bx*Sv - Pv_1*Sv;  
            v1_3 = V1cx*Su - Pu_1*Su; 
          
            Vector1 = [ v1_1 ; v1_2 ; v1_3 ];
         
            Matrix1(end,:) = [];
            Matrix1(2*N*Kx-(Kx),:) = []; 
            Matrix1(N*Kx-(Kx-1),:) = [];              
            Matrix1(1,:) = [];   
            
            Matrix1(:,2*(N*Kx-(Kx-1))) = [];
            Matrix1(:,N*Kx-(Kx-2)) = [];
            Matrix1(:,N*Kx-(Kx-1)) = [];
            Matrix1(:,1) = [];
                
            Vector1(end) = [];
            Vector1(2*N*Kx-(Kx)) = []; 
            Vector1(N*Kx-(Kx-1)) = [];
            Vector1(1) = [];
                
            R1 = Matrix1\Vector1;
            Rphi = R1(1:N*Kx-(Kx-1)-2);
            Rv = R1(N*Kx-(Kx-2)-2:2*(N*Kx-(Kx-1))-4);
            Ru = R1(2*(N*Kx-(Kx-1))+1 - 4: end);
                  
            Rphi = [0;Rphi;0];
            Rv= [0;Rv;0];
    
            R1 = [Rphi;Rv;Ru];
                 
            M2_11 = (Rphi'*Gt*Rphi)*M + (Rphi'*M*Rphi)*Gt;  
            M2_12 = (Rphi'*HTT*Ru)*M;
            M2_13 = (Rphi'*M*Rv)*HTT;
            M2_21 = (Ru'*TH*Rphi)*M; 
            M2_22 = (Ru'*TFtT*Ru)*M + (Ru'*TFtT*Ru)*Gt + (Ru'*Gt*Ru)*TFtT ;  
            M2_23 =  -(Ru'*TH*Rv)*HTT + (Ru'*HTT*Rv)*TH;  
            M2_31 = (Rv'*M*Rphi)*TH; 
            M2_32 = -(Rv'*HTT*Ru)*TH + (Rv'*TH*Ru)*HTT; 
            M2_33 = (Rv'*M*Rv)*TFtT + (Rv'*Gt*Rv)*TFtT + (Rv'*TFtT*Rv)*Gt;  
     
            Matrix2 = [  M2_11, M2_12, M2_13; 
                         M2_21, M2_22, M2_23  ;
                         M2_31, M2_32, M2_33  ] ;
         
            v2_1 = - Pphi_2'*Rphi;
            v2_2 = V2by*Ru - Pu_2'*Ru; 
            v2_3 = V2cy*Rv - Pv_2'*Rv; 
            
            Vector2 = [v2_1; v2_2 ; v2_3];
    
                
            Matrix2(end,:) = [];
            Matrix2(2*N*Kx-(Kx),:) = [];
            Matrix2(N*Kx-(Kx-1),:) = [];
            Matrix2(1,:) = [];
                
            Matrix2(:,2*(N*Kx-(Kx-1))) = [];
            Matrix2(:,N*Kx-(Kx-2)) = [];
            Matrix2(:,N*Kx-(Kx-1)) = [];
            Matrix2(:,1) = [];
                
            Vector2(end) = [];
            Vector2(2*N*Kx-(Kx)) = [];
            Vector2(N*Kx-(Kx-1)) = [];
            Vector2(1) = [];
                
            S1 = Matrix2\Vector2;
            Sphi = S1(1:N*Kx-(Kx-1)-2);
            Su = S1(N*Kx-(Kx-2)-2:2*(N*Kx-(Kx-1))-4);
            Sv = S1(2*(N*Kx-(Kx-1))+1 - 4: end);
                
            Sphi = [0;Sphi;0];
            Su= [0;Su;0];
    
            S1 = [Sphi;Su;Sv];
                
            counter=counter+1;
            disp(counter);
            if counter == max_iterations
                break
            end  
        end 
        cumcounter = [cumcounter + counter];
        % Append the current cumcounter value to the array
        cumcounter_values = [cumcounter_values, cumcounter];
        each_count_value = [each_count_value, counter];
        
        
        %setting XQ YQ from the calculated values
        xphi = Rphi ;
        xu= Ru;
        xv = Rv ;
        yphi = Sphi ;
        yu = Su ;
        yv = Sv ;
        
        
        %adding the known boundary values
        Xphi = [ Xphi , [  xphi ] ] ; 
        Xu = [ Xu, xu ] ;
        Xv = [ Xv, [  xv  ] ] ;
        Yphi = [ Yphi , [ yphi  ] ] ;
        Yu = [ Yu , [  yu  ] ] ;
        Yv = [ Yv, yv ] ;
        
        Q=Q+1;
             
        %calculating the approximations
        approxphi = zeros(N*Kx-(Kx-1)) ;
        approxu = zeros(N*Kx-(Kx-1)) ;
        approxv = zeros(N*Kx-(Kx-1)) ;
        for k = 1 :N*Kx-(Kx-1) 
            for l = 1 :N*Kx-(Kx-1) 
                for j = 1 : Q
                    approxphi(k,l) = approxphi(k,l) + Yphi(l,j)*Xphi(k,j); 
                    approxu(k,l) = approxu(k,l) + Yu(l,j)*Xu(k,j)*Kx ;
                    approxv(k,l) = approxv(k,l) + Yv(l,j)*Xv(k,j)*Kx ;
                end
            end
        end
            
        %working out the error
        PHIerror= zeros(N*Kx-(Kx-1)) ;
        PHImesh= zeros(N*Kx-(Kx-1)) ;
        Umesh = zeros(N*Kx-(Kx-1)) ;
        Vmesh = zeros(N*Kx-(Kx-1)) ;
        Uerror = zeros(N*Kx-(Kx-1)) ;
        Verror = zeros(N*Kx-(Kx-1)) ;
    
        for k = 1 : N*Kx-(Kx-1) 
            for l = 1 : N*Kx-(Kx-1) 
                
                Umesh(k,l) = U(xdomain(k),ydomain(l));
                Vmesh(k,l) = V(xdomain(k),ydomain(l));
                Uerror(k,l) = approxu(k,l) - Umesh(k,l);
                Verror(k,l) = approxv(k,l) - Vmesh(k,l);
                PHIerror(k,l) = approxphi(k,l) - phi(xdomain(k),ydomain(l));
                PHImesh(k,l) = phi(xdomain(k),ydomain(l));
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
        L2phi(Q,number) = L2errorphi^(1/2);
        L2U(Q,number) = L2erroru^(1/2);
        L2V(Q,number) = L2errorv^(1/2);
        
    end
    
    endTime = cputime;     
    enrichmentTime = endTime-startTime; 
    allb_cumcounter(number,:) = cumcounter_values;
    each_count_b(number,:) = each_count_value;
 

    T=GlobalTa; 
    
    Fxy = (2)*(pi^2)*sin(pi*xdomain).*(sin(pi*ydomain'));
    fxy = T*Fxy*T';
      
     D=GlobalDiff;
    
    for i=1:N*Kx-(Kx-1)
        w(i)= dx(1)/2 * M(i,i);
    end
    
    w(N)=0;
    w(2*N-1)=0;
     
    
    F=GlobalF;
    G=GlobalG;
    Ta=GlobalTa;
    Ybv = Ta*Yv;
    Xbu = Ta*Xu;
    Ybu = Ta*Yu;
    Xbv = Ta*Xv;
     
    for i=1:N*Kx-(Kx-1) -1
        for j=1:N*Kx-(Kx-1) 
            eq7_11(i,j) = Xbu(i)*Yu(j) + Yphi(j)*(Xphi(i+1)-Xphi(i));
        end
    end
    
    for i=1:N*Kx-(Kx-1)
        for j=1:N*Kx-(Kx-1)-1
            eq7_12(i,j) = Xv(i)*Ybv(j) + Xphi(i)*(Yphi(j+1)-Yphi(j));
        end
    end
    
    for i =1:N*Kx-(Kx-1)-1
        for j=1:N*Kx-(Kx-1)-1
            eq7_13(i,j) = (Xv(i+1)-Xv(i))*Ybv(j) - Xbu(i)*(Yu(j+1)- Yu(j));
            eq7_14(i,j) = (Xu(i+1)-Xu(i))*Ybu(j) + Xbv(i)*(Yv(j+1)- Yv(j)) - fxy(i,j)/Kx^2;
        end
    end
    
    eq7_16_all(number) = max(max((eq7_11)));
    eq7_17_all(number) = max(max((eq7_12)));
    eq7_18_all(number) = max(max((eq7_13)));
    eq7_19_all(number) = max(max((eq7_14)));
    

    summ = 0;
    for m = 1:N*Kx-(Kx-1)
        sumn = 0;
        for n=1:N*Kx-(Kx-1)
            sumj=0;
            for j =1:N*Kx-(Kx-1)
                sumj = sumj + G(m,j)*Xu(j);
            end
            suml = 0;
            for l =1:N*Kx-(Kx-1)-1
                suml = suml + F(n,l)*Ybu(l);
            end
            sumj2 = 0;
            for j =1:N*Kx-(Kx-1)-1
                sumj2 = sumj2 + F(m,j)*Xbv(j);
            end
            suml2 = 0;
            for l =1:N*Kx-(Kx-1)
                suml2 = suml2 + G(n,l)*Yv(l);
            end
            sumj3 = 0;
            for j =1:N*Kx-(Kx-1)-1
                suml3 = 0;
                for l=1:N
                    suml3 = suml3 + fxy(j,l);
                end
                sumj3 = sumj3 + suml3;
            end
            sumn = sumn + w(n)*(suml*sumj + sumj2*suml2 - sumj3);
        end
        summ = summ + sumn*w(m);
    end
    J1 = summ;
    J1= J1^2;

    summ = 0;
    for m = 1:N*Kx-(Kx-1)
        sumn = 0;
        for n=1:N*Kx-(Kx-1)
            sumj=0;
            for j =1:N*Kx-(Kx-1)-1
                sumj = sumj + F(m,j)*Xbu(j);
            end
            suml = 0;
            for l =1:N*Kx-(Kx-1)
                suml = suml + Yu(l);
            end
            sumj2 = 0;
            for j =1:N*Kx-(Kx-1)-1
                sumj2 = sumj2 + G(m,j)*Xphi(j);
            end
            suml2 = 0;
            for l =1:N*Kx-(Kx-1)
                suml2 = suml2 + Yphi(l);
            end
        sumn = sumn + w(n)*(suml*sumj + sumj2*suml2);   
        end
        summ = summ + w(m)*sumn;
    end
    J2 = summ;
    J2 = J2^2;
    
    summ = 0;
    for m = 1:N*Kx-(Kx-1)
        sumn = 0;
        for n=1:N*Kx-(Kx-1)
            sumj=0;
            for j =1:N*Kx-(Kx-1)
                sumj = sumj + Xphi(j);
            end
            suml = 0;
            for l =1:N*Kx-(Kx-1)
                suml = suml + G(n,l)*Yphi(l);
            end
            sumj2 = 0;
            for j =1:N*Kx-(Kx-1)
                sumj2 = sumj2 + Xv(j);
            end
            suml2 = 0;
            for l =1:N*Kx-(Kx-1)-1
                suml2 = suml2 + Ybv(l)*F(n,l);
            end
        sumn = sumn + w(n)*(suml*sumj + sumj2*suml2);   
        end
        summ = summ + w(m)*sumn;
    end
    J3 = summ;
    J3 = J3^2;
    
    
    summ = 0;
    for m = 1:N*Kx-(Kx-1)
        sumn = 0;
        for n=1:N*Kx-(Kx-1)
            sumj=0;
            for j =1:N*Kx-(Kx-1)
                sumj = sumj + Xv(j)*G(m,j);
            end
            suml = 0;
            for l =1:N*Kx-(Kx-1)-1
                suml = suml + F(n,l)*Ybv(l);
            end
            sumj2 = 0;
            for j =1:N*Kx-(Kx-1)-1
                sumj2 = sumj2 + Xbu(j)*F(m,j);
            end
            suml2 = 0;
            for l =1:N*Kx-(Kx-1)
                suml2 = suml2 + Yu(l)*G(n,l);
            end
            sumn = sumn + w(n)*(suml*sumj - sumj2*suml2);   
        end
        summ = summ + w(m)*sumn;
    end
    J4 = summ;
    J4 = J4^2;

    J1_all(number) = J1;
    J2_all(number) = J2;
    J3_all(number) = J3;
    J4_all(number) = J4;

end

L2Uvec = sqrt(L2U.^2 + L2V.^2);
 
fprintf('L2 errors for Example 11, i = 1, using the MLS SEM PGD, after one enrichment:')
fprintf('\n');
fprintf('\n'); 
fprintf('N =       ');
fprintf('%8d ', Nvalues);
fprintf('\n-----------------------------------------------------\n');
fprintf('L2 error φ   ');
fprintf('%6.2e ', L2phi(end,:));
fprintf('\nL2 error u   ');
fprintf('%6.2e ', L2Uvec(end,:));
fprintf('\n');
fprintf('\nJ1           ');
fprintf('%6.2e ', J1_all);
fprintf('\nJ2           ');
fprintf('%6.2e ', J2_all);
fprintf('\nJ3           ');
fprintf('%6.2e ', J3_all);
fprintf('\nJ4           ');
fprintf('%6.2e ', J4_all);
fprintf('\n');

fprintf('\n');
 fprintf('L2 errors for φ, u, and v for N = %d are given by %.3e, %.3e, and %.3e, respectively.\n', ...
    Nvalues(4), L2phi(4), L2U(4), L2V(4));

fprintf('\n');
fprintf('Values obtained for expressions (7.16)-(7.19) by the MLS SEM PGD for Example 11, i=1:')
fprintf('\n');
fprintf('\n'); 
fprintf('N = ');
fprintf('%8d ', Nvalues);
fprintf('\n-----------------------------------------------------\n');
fprintf('7.16   ');
fprintf('%6.2e ', eq7_16_all);
fprintf('\n7.17   ');
fprintf('%6.2e ', eq7_17_all);
fprintf('\n7.18   ');
fprintf('%6.2e ', eq7_18_all);
fprintf('\n7.19   ');
fprintf('%6.2e ', eq7_19_all);
fprintf('\n');

fprintf('\n');
 fprintf('For N=20, and epsilon = %.e, only %d ADFPA iterations were used for the first and only enrichment', ...
    epsilon, each_count_value);

fprintf('\n');
 fprintf('Time taken for the example with N=20 (in seconds):')
disp(enrichmentTime)
