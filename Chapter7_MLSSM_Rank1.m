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
phi=@(x,y) sin(pi*x)*sin(pi*y); 
U=@(x,y) - pi*cos(pi*x).*sin(pi*y) ; 
V=@(x,y) - pi*sin(pi*x).*cos(pi*y) ;
 
% N GLL points 
Nvalues = [8,12,16,20];

% Max number of iterations
max_iterations = 20; 
 
% Tolerance
epsilon = 1e-8;
 
warning('off', 'MATLAB:nearlySingularMatrix');

for number =1:length(Nvalues);
    cumcounter =0;
    cumcounter_values =[]; 
    each_count_value=[];

    N =Nvalues(number);
    N1=N+1;
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

    t_tilde=node; 
    for i=1:N 
        for k=1:N+1
            x_tilde(i,k) = ((node(i+1) - node(i))*t_tilde(k) + (node(i+1) + node(i)))/2;
        end
    end
    t= x_tilde;
      
    % T
    for i =1:N 
        for j=1:N+1
            prodk = 1; 
            sumk = 0; 
            for k = 1:N+1
                prodk = 1; 
                for l =1:N+1
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

    %sm, tn are GLL nodes -1 to 1
    sm = node;
    tn = node;
    
    % f tilde
    for i = 1:N
        for j = 1:N
            summ = 0;
            for m=1:N+1
                sumn = 0;
                for n=1:N+1
                    xco = ((node(i+1)-node(i))*sm(m) + (node(i+1)+node(i))) /2 ;
                    yco = ((node(j+1)-node(j))*tn(n) + (node(j+1)+node(j))) /2 ;
                    sumn = sumn + w(n)*( 2*pi^2*sin(pi*xco)*sin(pi*yco) );  % source term
                end
                summ = summ + sumn*w(m);
            end
            fxy(i,j) = ( ((node(i+1)-node(i))*(node(j+1)-node(j))) /4) *summ;
        end
    end
    
    % F
    F = zeros(N+1,N);
    for m=1:N+1
        for j= 1:N 
            sumk =0;
            for k = 1:j
                sumk = sumk + D(m,k);
            end
            F(m,j) = sumk;
        end
    end
    F=-F;

    % G
    G = zeros(N+1);
    G(:, 1) = -F(:,1);
    G(:, end) = F(:,end);
    G(:, 2:end-1) = D(:, 2:end-1);
    
    zeroN = zeros(N,1);
    zeroN1 = zeros(N+1,1);
    zerosN1=zeroN1;
    % H
    for i=1:N
        for j=1:N+1
            sumM =0;
            for m = 1:N+1
                sumM = sumM + w(m)*F(m,i)*G(m,j);
            end
            wFG(i,j)= sumM;
        end
    end
    H=wFG;
    HT=H'; 
    
    %F tilde
    for i=1:N
        for j=1:N
            sumM =0;
            for m = 1:N+1
                sumM = sumM + w(m)*F(m,i)*F(m,j);
            end
            Ft(i,j)= sumM;
        end
    end
    
    % G tilde
    for i=1:N+1
        for j=1:N+1
        sumM =0;
        for m = 1:N+1
            sumM = sumM + w(m)*G(m,i)*G(m,j);
        end
        Gt(i,j)= sumM;
        end
    end
    
    % M
    M=zeros(N+1,N+1);
    for i =1:N+1
        M(i,i) = w(i);
    end 
    Z = zeros(N+1,N+1);
     
    M2 = M(2:end-1,2:end-1);
    M4 = M(:,2:end-1);
    Gt2 = Gt(2:end-1,2:end-1);
    Gt4 = Gt(:,2:end-1);
    HT2 = HT(2:end-1,:);
    H2 = H(:,2:end-1);
    T2 = T(:,2:end-1);
     
    Q=0; 
    %Start approximations
    Xphi = zeros(N+1,Q ) ;
    Xu = zeros(N+1,Q ) ;
    Xv = zeros(N+1,Q ) ;
    Yphi = zeros(N+1,Q ) ;
    Yu = zeros(N+1,Q ) ;
    Yv = zeros(N+1,Q ) ;
    
    % Start enrichment loop
    for term = 1:no_enrichment
    
        startTime = cputime;
        
        %Initialize the Matrices containing the known max_iterations−1 enrichments
        Pphi_1 = zeros(N+1) ;
        Pu_1 = zeros(N+1) ;
        Pv_1 = zeros(N+1) ;
        Pphi_2 = zeros(N+1) ;
        Pu_2 = zeros(N+1) ;
        Pv_2 = zeros(N+1) ;
        XphiM = zeros(N+1,Q ) ;
        XphiGt = zeros(N+1,Q ) ;
        XvHTT = zeros(N+1,Q ) ;
        XuM = zeros(N+1,Q ) ;
        XphiTH = zeros(N+1,Q ) ;
        XvTFtT = zeros(N+1,Q ) ;
        XvGt = zeros(N+1,Q ) ;
        XuTH = zeros(N+1,Q ) ;
        XvTH = zeros(N+1,Q ) ;
        XuGt = zeros(N+1,Q ) ;
        XuTFtT = zeros(N+1,Q ) ;
        YphiGt = zeros(N+1,Q ) ;
        YphiM = zeros(N+1,Q ) ; 
        YvM = zeros(N+1,Q ) ;
        YvGt = zeros(N+1,Q ) ;
        YvTFtT = zeros(N+1,Q ) ;
        YuHTT = zeros(N+1,Q ) ;
        YuTH = zeros(N+1,Q ) ;
        YphiTH = zeros(N+1,Q ) ;
        YvTH = zeros(N+1,Q ) ;
        YvHTT = zeros(N+1,Q ) ;
        YuTFtT = zeros(N+1,Q ) ;
        YuGt = zeros(N+1,Q ) ;
        XuHTT = zeros(N+1,Q) ; 
        YuM = zeros(N+1,Q) ; 
        XvM = zeros(N+1,Q) ; 
        YvHT = zeros(N+1,Q) ; 
        YuHT = zeros(N+1,Q) ; 
        
        HTT = HT*T;
        TH = T'*H;
        TFtT = T'*Ft*T;
         
        for i = 1 :N+1
            for j = 1 : Q
                for k = 1 :N+1
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
         
         for i =1:N+1
             for k =1:N+1
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

        Pphi2_1 = Pphi_1(2:end-1,2:end-1);    
        Pv2_1 = Pv_1(2:end-1,2:end-1);
        Pu2_1 = Pu_1(2:end-1,2:end-1);
        Pphi2_2 = Pphi_2(2:end-1,2:end-1);    
        Pv2_2 = Pv_2(2:end-1,2:end-1);
        Pu2_2 = Pu_2(2:end-1,2:end-1);
             
            
        
        Sphi = rand(N-1,1); 
        Sv = rand(N+1,1); 
        Su=rand(N-1,1); 
        Rphi = rand(N-1,1); 
        Ru = rand(N+1,1); 
        Rv=rand(N-1,1); 
        
        N1=N+1;
        R1 = ones(3*N-1,1);
        S1 = ones(3*N-1,1);
        oldR = zeros(3*N-1,1);
        oldS = zeros(3*N-1,1);
        counter = 0;
        if term >1
            Sphi = Yphi(2:end-1,term-1);
            Su = Yu(2:end-1,term-1);
            Sv = Yv(:,term-1);
        else 
        end
        
        while max(max( abs (R1*S1' - oldR* oldS') ) ) >epsilon
            oldR = R1 ;
            oldS = S1 ;
            
            % Linear system 1
            Sbv = T*Sv;
            Sbu = T2*Su;
            
            M1_11 = (Sphi'*M2*Sphi)*Gt + (Sphi'*Gt2*Sphi)*M ; 
            M1_12 = (Sphi'*HT2*Sbv)*M ;
            M1_13 = (Sphi'*M2*Su)*HT *T;
            M1_21 = (Sbv'*H2*Sphi)*M ;
            M1_22 = (Sbv'*Ft*Sbv)*M + (Sbv'*Ft*Sbv)*Gt + (Sv'*Gt*Sv)*T'*Ft*T;
            M1_23 = -(Sbv'*H2*Su)*HT *T + (Sbu'*H*Sv)*T'*H;  
            M1_31 = (Su'*M2*Sphi)*T'*H; 
            M1_32 = -(Su'*HT2*Sbv)*T'*H + (Sbu'*H*Sv)*HT*T; 
            M1_33 = (Su'*Gt2*Su)*T'*Ft*T + (Su'*M2*Su)*T'*Ft*T + (Sbu'*Ft*Sbu)*Gt; 
            
            M1_11 = M1_11(2:end-1,2:end-1);
            M1_12 = M1_12(2:end-1,2:end-1);
            M1_13 = M1_13(2:end-1,:);
            M1_21 = M1_21(:,2:end-1);
            M1_22 = M1_22(:,2:end-1);
            M1_31 = M1_31(2:end-1,2:end-1); 
            M1_32 = M1_32(2:end-1,2:end-1);
            M1_33 = M1_33(2:end-1,:);
            
            Matrix1 = [  M1_11, M1_12, M1_13; 
                         M1_21, M1_22, M1_23  ;
                         M1_31, M1_32, M1_33  ] ;
    
            for b=1:N+1
                sumi=0;
                for i =1:N
                    sumj=0;
                    for j=1:N
                        sumk=0;
                        for k =1:N+1
                            suml =0;
                            for l=1:N
                                suml = suml + HT(k,l)*fxy(j,l);
                            end
                            sumk = sumk + suml*Sv(k);
                        end
                        sumj = sumj + sumk*Ft(i,j);
                    end
                    sumi = sumi + sumj*T(i,b);
                end
                V1_2(b) = sumi;
            end
            V1_2 = V1_2(:); 
           
            for i =1:N-1
                sumj=0;
                for j=1:N
                    sumb=0;
                    for b = 1:N-1
                        sumk=0;
                        for k =1:N
                            suml =0;
                            for l=1:N
                                suml = suml + Ft(k,l)*fxy(j,l);
                            end
                            sumk = sumk + suml*T(k,b+1);
                        end
                        sumb = sumb + sumk*Su(b);
                    end
                    sumj = sumj + sumb*HT(i+1,j);
                end
                V1_3(i) = sumj;
            end
           
            V1_3 = V1_3(:);
          
            V1_1 = zeros((N-1),1);
            v1_1 = V1_1 - Pphi2_1*Sphi;
            v1_2 = V1_2 - Pv_1*Sv;
            v1_3 = V1_3 - Pu2_1*Su;
           
            Vector1 = [ v1_1 ; v1_2 ; v1_3 ];
          
            R1 = Matrix1\Vector1;
            Rphi = R1(1:N-1,:);
            Rv = R1(N:2*N-2,:);
            Ru = R1(2*N-1:end,:);
        
            % Linear system 2
            Rbv = T2*Rv;
            Rbu = T*Ru;
           
            M2_11 = (Rphi'*Gt2*Rphi)*M + (Rphi'*M2*Rphi)*Gt; 
            M2_12 = (Rphi'*HT2*Rbu)*M;  
            M2_13 = (Rphi'*M2*Rv)*HT*T;  
            M2_21 = (Rbu'*H2*Rphi)*M;    
            M2_22 = (Rbu'*Ft*Rbu)*M + (Rbu'*Ft*Rbu)*Gt + (Ru'*Gt*Ru)*T'*Ft*T ; 
            M2_23 =  -(Rbu'*H2*Rv)*HT*T + (Ru'*HT*Rbv)*T'*H; 
            M2_31 = (Rv'*M2*Rphi)*T'*H; 
            M2_32 = -(Rv'*HT2*Rbu)*T'*H + (Rbv'*H*Ru)*HT*T; 
            M2_33 = (Rv'*M2*Rv)*T'*Ft*T + (Rv'*Gt2*Rv)*T'*Ft*T + (Rbv'*Ft*Rbv)*Gt; 
 
            M2_11 = M2_11(2:end-1,2:end-1);
            M2_12 = M2_12(2:end-1,2:end-1);
            M2_13 = M2_13(2:end-1,:);
            M2_21 = M2_21(:,2:end-1);
            M2_22 = M2_22(:,2:end-1);
            M2_31 = M2_31(2:end-1,2:end-1); 
            M2_32 = M2_32(2:end-1,2:end-1);
            M2_33 = M2_33(2:end-1,:);
            
            Matrix2 = [  M2_11, M2_12, M2_13; 
                         M2_21, M2_22, M2_23  ;
                         M2_31, M2_32, M2_33  ] ;
            
            for b =1:N+1 
                sumk=0;
                for k =1:N
                    suml=0;
                    for l=1:N
                        sumi=0;
                        for i =1:N1
                            sumj =0;
                            for j=1:N
                                sumj = sumj + HT(i,j)*fxy(j,l);
                            end
                            sumi = sumi + sumj*Ru(i);
                        end
                        suml = suml + sumi*Ft(k,l);
                    end
                    sumk = sumk + suml*T(k,b);
                end
                V2_2(b) = sumk;
            end
            V2_2 = V2_2(:); 
            
            for k =1:N-1
                suml=0;
                for l=1:N
                    sumb=0;
                    for b =1:N-1
                        sumi =0;
                        for i=1:N
                            sumj =0;
                            for j=1:N
                                sumj = sumj + Ft(i,j)*fxy(j,l);
                            end
                            sumi = sumi + sumj*T(i,b+1);
                        end
                        sumb = sumb + sumi*Rv(b);
                    end
                    suml = suml + sumb*HT(k+1,l);
                end
                V2_3(k) = suml;
            end
            V2_3 = V2_3(:);
            
              
            V2_1 = zeros((N-1),1);
            v2_1 = V2_1 - Pphi2_2'*Rphi;
            v2_2 = V2_2 - Pu_2'*Ru;
            v2_3 = V2_3 - Pv2_2'*Rv;
            
            Vector2 = [v2_1; v2_2 ; v2_3];

            S1 = Matrix2\Vector2;
            Sphi = S1(1:N-1,:);
            Su = S1(N:2*N-2,:);
            Sv = S1(2*N-1:end,:);
            
            
            R1 = [Rphi; Rv; Ru];
            S1 = [Sphi; Su; Sv];

            counter=counter+1;
            
            if counter == max_iterations
                
                break
            end  
        end       
  
        cumcounter = [cumcounter + counter];
        % Append the current cumcounter value to the array
        cumcounter_values = [cumcounter_values, cumcounter];
        each_count_value = [each_count_value, counter];
      
        xphi = Rphi ;
        xu= Ru;
        xv = Rv ;
        yphi = Sphi ;
        yu = Su ;
        yv = Sv ;

        %adding the known boundary values
        Xphi = [ Xphi , [ 0 ; xphi ; 0 ] ] ;
        Xu = [ Xu, xu ] ;
        Xv = [ Xv, [ 0 ; xv ; 0 ] ] ;
        Yphi = [ Yphi , [ 0 ; yphi ; 0 ] ] ;
        Yu = [ Yu , [ 0 ; yu ; 0 ] ] ;
        Yv = [ Yv, yv ] ;
        
        Q=Q+1;

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
        PHIerror= zeros(N+1) ;
        PHImesh= zeros(N+1) ;
        Umesh = zeros(N+1) ;
        Vmesh = zeros(N+1) ;
        Uerror = zeros(N+1) ;
        Verror = zeros(N+1) ;

        for k = 1 :N+1
            for l = 1 :N+1
                err = approxphi(k,l)- phi(node(k),node(l));
                Uerror(k,l) = approxu(k,l) - U(node(k),node(l));
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
               QL2phierror=QL2phierror+((w(m)*w(num)*(approxphi(m,num)-PHImesh(m,num))^2)) ;
               QL2uerror=QL2uerror+((w(m)*w(num)*(approxu(m,num)-Umesh(m,num))^2)) ;
               QL2verror=QL2verror+((w(m)*w(num)*(approxv(m,num)-Vmesh(m,num))^2)) ;
            end
        end
        QL2phierror= QL2phierror^(1/2);
        QL2uerror= QL2uerror^(1/2);
        QL2verror= QL2verror^(1/2);
         
        QL2phi(Q,number)=QL2phierror;
  
    end
    endTime = cputime;       
    
    % Calculate execution time for this iteration
    enrichmentTime = endTime - startTime;

    L2U(number)= QL2uerror;
    L2V(number)= QL2verror;
    L2phi(number)= QL2phierror;
    
    allb_cumcounter(number,:) = cumcounter_values;
    each_count_b(number,:) = each_count_value;

    Ybv = T*Yv;
    Xbu = T*Xu;
    Ybu = T*Yu;
    Xbv = T*Xv; 

    summ = 0;
    for m = 1:N+1
        sumn = 0;
        for n=1:N+1
            sumj=0;
            for j =1:N+1
                sumj = sumj + G(m,j)*Xu(j);
            end
            suml = 0;
            for l =1:N
                suml = suml + F(n,l)*Ybu(l);
            end
            sumj2 = 0;
            for j =1:N
                sumj2 = sumj2 + F(m,j)*Xbv(j);
            end
            suml2 = 0;
            for l =1:N+1
                suml2 = suml2 + G(n,l)*Yv(l);
            end
            sumj3 = 0;
            for j =1:N
                suml3 = 0;
                for l=1:N
                    suml3 = suml3 + fxy(j,l)*F(n,l);
                end
                sumj3 = sumj3 + suml3*F(m,j);
            end
        sumn = sumn + w(n)*(suml*sumj + sumj2*suml2 - sumj3);
        end
        summ = summ + sumn*w(m);
    end
    J1 = summ;
    J1= J1^2;
    
    summ = 0;
    for m = 1:N+1
        sumn = 0;
        for n=1:N+1
            sumj=0;
            for j =1:N
                sumj = sumj + F(m,j)*Xbu(j);
            end
            suml = 0;
            for l =1:N+1
                suml = suml + Yu(l);
            end
            sumj2 = 0;
            for j =1:N
                sumj2 = sumj2 + G(m,j)*Xphi(j);
            end
            suml2 = 0;
            for l =1:N+1
                suml2 = suml2 + Yphi(l);
            end
        sumn = sumn + w(n)*(suml*sumj + sumj2*suml2);   
        end
        summ = summ + w(m)*sumn;
    end
    J2 = summ;
    J2 = J2^2;
    
    summ = 0;
    for m = 1:N+1
        sumn = 0;
        for n=1:N+1
            sumj=0;
            for j =1:N+1
                sumj = sumj + Xphi(j);
            end
            suml = 0;
            for l =1:N+1
                suml = suml + G(n,l)*Yphi(l);
            end
            sumj2 = 0;
            for j =1:N+1
                sumj2 = sumj2 + Xv(j);
            end
            suml2 = 0;
            for l =1:N
                suml2 = suml2 + Ybv(l)*F(n,l);
            end
        sumn = sumn + w(n)*(suml*sumj + sumj2*suml2);   
        end
        summ = summ + w(m)*sumn;
    end
    J3 = summ;
    J3 = J3^2;
    
    
    summ = 0;
    for m = 1:N+1
        sumn = 0;
        for n=1:N+1
            sumj=0;
            for j =1:N+1
                sumj = sumj + Xv(j)*G(m,j);
            end
            suml = 0;
            for l =1:N
                suml = suml + F(n,l)*Ybv(l);
            end
            sumj2 = 0;
            for j =1:N
                sumj2 = sumj2 + Xbu(j)*F(m,j);
            end
            suml2 = 0;
            for l =1:N+1
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
 
    for i=1:N
        for j=1:N+1
            eq7_11(i,j) = Xbu(i)*Yu(j) + Yphi(j)*(Xphi(i+1)-Xphi(i));
        end
    end
    
    for i=1:N+1
        for j=1:N
            eq7_12(i,j) = Xv(i)*Ybv(j) + Xphi(i)*(Yphi(j+1)-Yphi(j));
        end
    end
    
    for i =1:N
        for j=1:N
            eq7_13(i,j) = (Xv(i+1)-Xv(i))*Ybv(j) - Xbu(i)*(Yu(j+1)- Yu(j));
            eq7_14(i,j) = (Xu(i+1)-Xu(i))*Ybu(j) + Xbv(i)*(Yv(j+1)- Yv(j)) - fxy(i,j);
            
        end
    end
    
    eq7_16_all(number) = max(max((eq7_11)));
    eq7_17_all(number) = max(max((eq7_12)));
    eq7_18_all(number) = max(max((eq7_13)));
    eq7_19_all(number) = max(max((eq7_14)));
end
 
L2Uvec = sqrt(L2U.^2 + L2V.^2);
 
fprintf('L2 errors for Example 11, i = 1, using the MLS SM PGD, after one enrichment:')
fprintf('\n');
fprintf('\n'); 
fprintf('N =       ');
fprintf('%8d ', Nvalues);
fprintf('\n-----------------------------------------------------\n');
fprintf('L2 error φ   ');
fprintf('%6.2e ', L2phi);
fprintf('\nL2 error u   ');
fprintf('%6.2e ', L2Uvec);
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
fprintf('Values obtained for expressions (7.16)-(7.19) by the MLS SM PGD for Example 11, i=1:')
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

 