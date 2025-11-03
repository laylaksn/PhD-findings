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


% Define the thresholds
thresholds = [0.001, 0.0001, 0.00001, 0.000001, 0.0000001];

results = zeros(20, length(thresholds));
 
syms x y

format short
no_enrichment=25;

% number of GLL points
Nvalues = [ 8,12,16,20];

% max iterations
stopcrit = 300;
%tolerence
epsilon = 10^(-20);

for no_runs = 1:20 
    approxphiQ=[]; 
    
    Lphi=zeros(length(Nvalues),1);
    LU=zeros(length(Nvalues),1);
    LV=zeros(length(Nvalues),1);


    for number=1:length(Nvalues)
        approxphiQ=[];
        N=Nvalues(number); 
        
        % source term
        F = ones(N+1);

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
      
        B=zeros(N1);
        for i=1:N1
            for j=1:N1
                sumB=0;
                for n=1:N1
                sumB=sumB+(w(n)*w(j)*D(n,i)*F(n,j));
                end
                B(i,j)=sumB;
            end
        end
        C=zeros(N1);
        for i=1:N1
            for j=1:N1
                sumC=0;
                for n=1:N1
                 sumC=sumC+w(i)*w(n)*F(i,n)*D(n,j);
                end
                C(i,j)=sumC;
            end
        end 
    
        A = zeros(N1); M = zeros(N1); L = zeros(N1);
        for i=1:N1
            for j=1:N1
                nn=1:N1;
                A(i,j)=sum(w(nn).*D(nn,i).*D(nn,j));           
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
    
        %Introducing (N−1)^2 matrices
        A2 = A( 2 :N, 2 :N) ;
        M2 = M( 2 :N, 2 :N) ;
        L2 = L ( 2 :N, 2 :N) ;
        B2 = B ( 2 :N, 2 :N) ;
        D2 = C( 2 :N, 2 :N) ;
        
        %Introducing (N−1) *(N+1) sized matrices
        M3 = M( 2 :N, : ) ;
        A3 = A( 2 :N, : ) ;
        L3 = L ( 2 :N, : ) ;
        B3 = B ( 2 :N, : ) ;
        D3 = C( 2 :N, : ) ;
        
        %Introducing (N+1) *(N−1) sized matrix
        L4 = L ( : , 2 :N) ;
        M4 = M( : , 2 :N) ;
        A4 = A ( : , 2 :N) ;
        B4 = B ( : , 2 :N) ;
        D4 = C ( : , 2 :N) ;
        %Start off Q
        Q=0;
        
        %Start approximations
        Xphi = zeros(N+1,Q ) ;
        Xu = zeros(N+1,Q ) ;
        Xv = zeros(N+1,Q ) ;
        Yphi = zeros(N+1,Q ) ;
        Yu = zeros(N+1,Q ) ;
        Yv = zeros(N+1,Q ) ;
    
        cumcounter =0;
        %Qth enrichments
        for term = 1:no_enrichment
        
            startTime = cputime;
            
            %Initialize the Matrices containing the known Q−1 enrichments
            Pphi = zeros(N-1) ;
            Pu = zeros(N+1,N-1) ;
            Pv = zeros(N-1,N+1) ;
            
            XphiA = zeros(N-1,Q ) ;
            XphiM= zeros(N-1,Q ) ;
            XphiL= zeros(N+1,Q ) ;
            XuA= zeros(N+1,Q ) ;
            XuM= zeros(N+1,Q ) ;
            XuL= zeros(N-1,Q ) ;
            XuLT = zeros(N-1,Q ) ;
            XvA= zeros(N-1,Q ) ;
            XvM= zeros(N-1,Q ) ;
            XvL= zeros(N+1,Q ) ;
            XvLT = zeros(N+1,Q ) ;            
            YphiA= zeros(N-1,Q ) ;
            YphiM= zeros(N-1,Q ) ;
            YphiL= zeros(N+1,Q ) ;
            YuA= zeros(N-1,Q ) ;
            YuM= zeros(N-1,Q ) ;
            YuL= zeros(N+1,Q ) ;
            YuLT = zeros(N+1,Q ) ;
            YvA= zeros(N+1,Q ) ;
            YvL= zeros(N-1,Q ) ;
            YvLT = zeros(N-1,Q ) ;
            YvM= zeros(N+1,Q ) ;
            
            for i = 1 :N-1
                for j = 1 : Q
                     for k = 1 :N+1
                         XphiA ( i , j ) = XphiA ( i , j ) + Xphi ( k , j ) *A4( k , i ) ;
                         XphiM( i , j ) = XphiM ( i , j ) + Xphi ( k , j ) *M4( k , i ) ;
                         XuL( i , j ) = XuL( i , j ) + Xu( k , j ) *L3 ( i , k ) ;
                         XuLT( i , j ) = XuLT( i , j ) + Xu( k , j ) *L4 ( k , i ) ;
                         XvA( i , j ) = XvA( i , j ) + Xv( k , j ) *A4( k , i ) ;
                         XvM( i , j ) = XvM( i , j ) + Xv( k , j ) *M4( k , i ) ;
                         YphiA ( i , j ) = YphiA ( i , j ) + Yphi ( k , j ) *A4( k , i ) ;
                         YphiM( i , j ) = YphiM ( i , j ) + Yphi ( k , j ) *M4( k , i ) ;
                         YuA( i , j ) = YuA( i , j ) + Yu( k , j ) *A4( k , i ) ;
                         YuM( i , j ) = YuM( i , j ) + Yu( k , j ) *M4( k , i ) ;
                         YvL( i , j ) = YvL( i , j ) + Yv( k , j ) *L3 ( i , k ) ;
                         YvLT( i , j ) = YvLT( i , j ) + Yv( k , j ) *L4 ( k , i ) ;
        
                     end
                end
            end
        
            for i = 1 :N+1
                 for j = 1 : Q
                     for k = 1 :N+1
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
                     end
                 end
            end
            
            
           
        
            for i = 1 :N-1
                for j = 1 : Q
                    for k = 1 :N-1 
                        Pphi(k,i) = Pphi(k,i) + (XuLT(k,j)*YuM(i,j) + XphiA(k,j)*YphiM(i,j) + XvM(k,j)*YvLT(i,j) + XphiM(k,j)*YphiA(i,j)   );
                    end
                end
            end
        
            for i = 1 :N-1
                for j = 1 : Q
                    for k = 1 :N+1
                        Pu(k,i) = Pu(k,i) + (XuA(k,j)*YuM(i,j) + XvLT(k,j)*YvL(i,j) + XuM(k,j)*YuM(i,j) + XphiL(k,j)*YphiM(i,j) - XvL(k,j)*YvLT(i,j) + XuM(k,j)*YuA(i,j)) ;
                    end
                end
            end
            
            for i = 1 :N+1
                for j = 1 : Q
                    for k = 1 :N-1 
                         Pv(k,i) = Pv(k,i) + (XvM(k,j)*YvA(i,j) + XuL(k,j)*YuLT(i,j) + XvM(k,j)*YvM(i,j) + XphiM(k,j)*YphiL(i,j) - XuLT(k,j)*YuL(i,j) + XvA(k,j)*YvM(i,j));
                    end
                end
            end
        
            %Unknown functions for enrichment
            R1 = ones( 3*N-1 ,1) ;
            S1 = ones( 3*N-1 ,1) ;
            oldR = zeros( 3*N-1 ,1) ;
            oldS = zeros( 3*N-1 ,1) ;
            Rphi = rand(N-1 ,1) ;
            Ru = rand(N+1, 1 ) ;
            Rv = rand(N-1, 1 ) ;
            
            counter = 0 ;
    
            %enrichment Stage
            while max(max( abs (R1*S1' - oldR* oldS') ) ) >epsilon
                oldR = R1 ;
                oldS = S1 ;
                
                Mr11 = (Rphi'*A2*Rphi)*M2 + (Rphi'*M2*Rphi)*A2;
                Mr12 = (Rphi'*L4'*Ru)*M2  ;
                Mr13 = (Rphi'*M2*Rv)*L4'   ;
                Mr21 = (Ru'*L4*Rphi)*M2  ;
                Mr22 = (Ru'*A*Ru) *M2 + (Ru'*M*Ru) *A2 + (Ru'*M*Ru) *M2 ;
                Mr23 = (Ru'*L3'*Rv)*L3 - (Ru'*L4*Rv)*L4';
                Mr31 = (Rv'*M2*Rphi)*L4   ;
                Mr32 = (Rv'*L3*Ru)*L3' - (Rv'*L4'*Ru)*L4 ;
                Mr33 = (Rv'*A2*Rv)*M + (Rv'*M2*Rv)*A + (Rv'*M2*Rv)*M;
    
                Matrix2 = [Mr11, Mr12, Mr13; Mr21, Mr22, Mr23; Mr31, Mr32, Mr33] ;
                Vr1= - Pphi'*Rphi;
                Vr2 = B4'*Ru-Pu'*Ru;
                Vr3 = D3'*Rv-Pv'*Rv ;
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
                
                Vs1 = -Pphi*Sphi ;
                Vs2 = B4*Su-Pu*Su ; 
                Vs3 = D3*Sv-Pv*Sv ;
                Vector1 = [ Vs1 ; Vs2 ; Vs3 ] ;
                
                R1 = Matrix1\Vector1 ;
                Rphi = R1(1:N-1);
                Ru = R1(N:2*N);
                Rv = R1(2*N+1:3*N-1);
                
                counter=counter+1; 
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
            Xphi = [ Xphi , [ 0 ; xphi ; 0 ] ] ;
            Xu = [ Xu, xu ] ;
            Xv = [ Xv, [ 0 ; xv ; 0 ] ] ;
            Yphi = [ Yphi , [ 0 ; yphi ; 0 ] ] ;
            Yu = [ Yu , [ 0 ; yu ; 0 ] ] ;
            Yv = [ Yv, yv ] ;
            Q= Q+1;    
            
            %calculating the approximations
            approxphi = zeros(N+1) ;
            approxu = zeros(N+1) ;
            approxv = zeros(N+1) ;
            approxphi_1 = zeros(N+1) ;
            approxu_1 = zeros(N+1) ;
            approxv_1 = zeros(N+1) ;
            for k = 1 :N+1
                for l = 1 :N+1
                    for j = 1 : Q
                        approxphi(k,l) = approxphi(k,l)+Yphi(k,j)*Xphi(l,j);
                        approxu ( k , l ) = approxu ( k , l ) +Yu( k , j ) *Xu( l , j ) ;
                        approxv ( k , l ) = approxv ( k , l ) + Yv( k , j ) *Xv( l , j ) ;
                    end
                end
            end
            
            for k = 1 :N+1
                for l = 1 :N+1
                    for j = Q
                        approxphi_1(k,l) = Yphi(k,j)*Xphi(l,j);
                        approxu_1 ( k , l ) = Yu( k , j ) *Xu( l , j ) ;
                        approxv_1 ( k , l ) = Yv( k , j ) *Xv( l , j ) ;
                    end
                end
            end
    
            approxphiQ(:,:, Q) = approxphi_1;

            PhiNorm=0;
            UNorm=0;
            VNorm=0;
            
            for num=1:N1
            for m=1:N1
                PhiNorm=PhiNorm+((w(m)*w(num)*(approxphi_1(m,num))^2)) ;
                UNorm=UNorm+((w(m)*w(num)*(approxu_1(m,num))^2)) ;
                VNorm=VNorm+((w(m)*w(num)*(approxv_1(m,num))^2)) ;
            end
            end
            PhiNorm= PhiNorm^(1/2);
            UNorm= UNorm^(1/2);
            VNorm= VNorm^(1/2);
            
            QLphi(Q,number)=PhiNorm;
            QLu(Q,number)=UNorm;
            QLv(Q,number)=VNorm;
    
            QBeLphi(Q,no_runs)= PhiNorm;
        end
    
    end
    
    phinorm = [QLphi]; 

                
end

average_PhiNorm = mean(QBeLphi, 2, 'omitnan');  % 25x1 average

disp('Due to small deviations in each run, the code was run 20 times and the average number of enrichments to reach each threshold is computed.')

average_enrichment = nan(length(thresholds), 1);

for i = 1:length(thresholds)
    % For decreasing PhiNorm, find when it FIRST goes below the threshold
    idx = find(average_PhiNorm <= thresholds(i), 1, 'first');
    if ~isempty(idx)
        average_enrichment(i) = idx;
        fprintf('Average number of enrichments for threshold %.e: %.f\n', thresholds(i), idx);
    else
        fprintf('Threshold %.7f was never reached in any iteration\n', thresholds(i));
    end
end


% % Calculate the average enrichment number for each threshold, ignoring NaN values
% average_enrichment = mean(results, 1);
% 
% disp('Due to small deviations in each run, the code was run 20 times and the average number of enrichments to reach each threshold is computed.')
% % Display the results
% for i = 1:length(thresholds)
%     if ~isnan(average_enrichment(i))
%         fprintf('Average number of enrichments for threshold %.e: %.f\n', thresholds(i), average_enrichment(i));
%     else
%         fprintf('Threshold %.7f was never reached in any iteration\n', thresholds(i));
%     end
% end

figure;
set(gcf, 'Color', 'white')
subplot(3,4,1);
surf(node,node,approxphiQ(:,:,1)); 
title('$\phi_1$','interpreter','latex')
subplot(3,4,2);
surf(node,node,approxphiQ(:,:,2)); 
title('$\phi_2$','interpreter','latex')
subplot(3,4,3);
surf(node,node,approxphiQ(:,:,3)); 
title('$\phi_3$','interpreter','latex')
subplot(3,4,4);
surf(node,node,approxphiQ(:,:,4)); 
title('$\phi_4$','interpreter','latex')
subplot(3,4,5);
surf(node,node,approxphiQ(:,:,5)); 
title('$\phi_5$','interpreter','latex')
subplot(3,4,6);
surf(node,node,approxphiQ(:,:,6)); 
title('$\phi_6$','interpreter','latex')
subplot(3,4,7);
surf(node,node,approxphiQ(:,:,7)); 
title('$\phi_7$','interpreter','latex')
subplot(3,4,8);
surf(node,node,approxphiQ(:,:,8)); 
title('$\phi_8$','interpreter','latex')
subplot(3,4,9);
surf(node,node,approxphiQ(:,:,9)); 
title('$\phi_9$','interpreter','latex')
subplot(3,4,10);
surf(node,node,approxphiQ(:,:,10)); 
title('$\phi_{10}$','interpreter','latex')
subplot(3,4,11);
surf(node,node,approxphiQ(:,:,11)); 
title('$\phi_{11}$','interpreter','latex')
subplot(3,4,12);
surf(node,node,approxphiQ(:,:,12)); 
title('$\phi_{12}$','interpreter','latex')

 
N_values=Nvalues';

     
figure;
set(gcf, 'Color', 'white')
set(gca,'fontsize',14);
handles = semilogy(1:no_enrichment,QLphi, 'Linewidth', 2);
 legend('$N=8$','$N=12$','$N=16$','$N=20$','interpreter','latex')
grid on
xlabel('$q$','interpreter','latex','fontsize',14)
xlim([1 no_enrichment])
%xticks(1:1:no_enrichment)
ylabel('Norm of approximations for $\phi$','interpreter','latex','fontsize',14)
ylabel('$\Vert \phi_q \Vert $','interpreter','latex','fontsize',14)