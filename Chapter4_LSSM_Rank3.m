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

syms x y
format short

% Number of enrichments
no_enrichment=5;

% N GLL points
Nvalues=[8,12,16,20] ;

% Max number of iterations
max_iterations = 100;

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
    F = F -2*(node.^2)-2*(node'.^2)+4;
    F = F + (8)*(pi^2)*sin(2*pi*node).*(sin(2*pi*node'));

    % true solutions  
    syms sx sy
    phi=@(sx,sy) sin(pi*sx)*sin(pi*sy)+(sx.^2-1)*(sy.^2-1) + sin(2*pi*sx).*sin(2*pi*sy); 
    U=@(x,y) -(2)*x.*(y.^2-1) - pi*cos(pi*x).*sin(pi*y) - 2*pi*cos(2*pi*x).*sin(2*pi*y); 
    V=@(x,y) -(2)*y.*(x.^2-1) - pi*sin(pi*x).*cos(pi*y) - 2*pi*sin(2*pi*x).*cos(2*pi*y);

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
     
    startTime = cputime;
    %Qth enrichments
    for term = 1:no_enrichment
        
        
        
        %Initialize the Matrices containing the known Q−1 enrichments
        Pphi = zeros(N-1);
        Pu = zeros(N+1,N-1);
        Pv = zeros(N-1,N+1);
        
        
        XphiA = zeros(N-1, Q);
        XphiM = zeros(N-1, Q);
        XphiL = zeros(N+1, Q);
        
        XuA   = zeros(N+1, Q);
        XuM   = zeros(N+1, Q);
        XuL   = zeros(N-1, Q);
        XuLT  = zeros(N-1, Q);
        
        XvA   = zeros(N-1, Q);
        XvM   = zeros(N-1, Q);
        XvL   = zeros(N+1, Q);
        XvLT  = zeros(N+1, Q);
        
        YphiA = zeros(N-1, Q);
        YphiM = zeros(N-1, Q);
        YphiL = zeros(N+1, Q);
        
        YuA   = zeros(N-1, Q);
        YuM   = zeros(N-1, Q);
        YuL   = zeros(N+1, Q);
        YuLT  = zeros(N+1, Q);
        
        YvA   = zeros(N+1, Q);
        YvM   = zeros(N+1, Q);
        YvL   = zeros(N-1, Q);
        YvLT  = zeros(N-1, Q);

        for i = 1:N-1
            for j = 1:Q
                for k = 1:N+1
                    XphiA(i, j) = XphiA(i, j) + Xphi(k, j) * A4(k, i);
                    XphiM(i, j) = XphiM(i, j) + Xphi(k, j) * M4(k, i);
                    XuL(i, j)   = XuL(i, j)   + Xu(k, j)   * L3(i, k);
                    XuLT(i, j)  = XuLT(i, j)  + Xu(k, j)   * L4(k, i);
                    XvA(i, j)   = XvA(i, j)   + Xv(k, j)   * A4(k, i);
                    XvM(i, j)   = XvM(i, j)   + Xv(k, j)   * M4(k, i);
                    YphiA(i, j) = YphiA(i, j) + Yphi(k, j) * A4(k, i);
                    YphiM(i, j) = YphiM(i, j) + Yphi(k, j) * M4(k, i);
                    YuA(i, j)   = YuA(i, j)   + Yu(k, j)   * A4(k, i);
                    YuM(i, j)   = YuM(i, j)   + Yu(k, j)   * M4(k, i);
                    YvL(i, j)   = YvL(i, j)   + Yv(k, j)   * L3(i, k);
                    YvLT(i, j)  = YvLT(i, j)  + Yv(k, j)   * L4(k, i);
                end
            end
        end

        for i = 1:N+1
            for j = 1:Q
                for k = 1:N+1
                    XphiL(i, j) = XphiL(i, j) + Xphi(k, j) * L(i, k);
                    XuA(i, j)   = XuA(i, j)   + Xu(k, j)   * A(k, i);
                    XuM(i, j)   = XuM(i, j)   + Xu(k, j)   * M(k, i);
                    XvL(i, j)   = XvL(i, j)   + Xv(k, j)   * L(i, k);
                    XvLT(i, j)  = XvLT(i, j)  + Xv(k, j)   * L(k, i);
                    YphiL(i, j) = YphiL(i, j) + Yphi(k, j) * L(i, k);
                    YuL(i, j)   = YuL(i, j)   + Yu(k, j)   * L(i, k);
                    YuLT(i, j)  = YuLT(i, j)  + Yu(k, j)   * L(k, i);
                    YvA(i, j)   = YvA(i, j)   + Yv(k, j)   * A(k, i);
                    YvM(i, j)   = YvM(i, j)   + Yv(k, j)   * M(k, i);
                end
            end
        end

        for i = 1 :N-1
            for j = 1 : Q
                for k = 1 :N-1 
                    Pphi(k,i) = Pphi(k,i) + (XuLT(k,j)*YuM(i,j) + XphiA(k,j)*YphiM(i,j) + XvM(k,j)*YvLT(i,j) + XphiM(k,j)*YphiA(i,j) );
                end
            end
        end
    
        for i = 1 :N-1
            for j = 1 : Q
                for k = 1 :N+1
                    Pu(k,i) = Pu(k,i) +(XuA(k,j)*YuM(i,j) + XvLT(k,j)*YvL(i,j) + XuM(k,j)*YuM(i,j) + XphiL(k,j)*YphiM(i,j) - XvL(k,j)*YvLT(i,j) + XuM(k,j)*YuA(i,j)) ;
                end
            end
        end
        
        for i = 1 :N+1
            for j = 1 : Q
                for k = 1 :N-1
                     Pv(k,i) = Pv(k,i) + ( XvM(k,j)*YvA(i,j) + XuL(k,j)*YuLT(i,j) + XvM(k,j)*YvM(i,j) + XphiM(k,j)*YphiL(i,j) - XuLT(k,j)*YuL(i,j) + XvA(k,j)*YvM(i,j));
                end
            end
        end
        
        %Unknown functions for enrichment
        R1 = ones( 3*N-1 ,1) ;
        S1 = ones( 3*N-1 ,1) ;
        oldR = zeros( 3*N-1 ,1) ;
        oldS = zeros( 3*N-1 ,1) ;
        Rphi = rand(N-1 ,1) ;        Ru = rand(N+1, 1 ) ;        Rv = rand(N-1, 1 ) ;
        Sphi = rand(N-1, 1 ) ;       Su = rand(N-1 ,1) ;         Sv = rand(N+1 ,1) ;
        counter = 0 ;
        
        if term >1
            Rphi = Xphi(2:end-1,term-1);
            Ru = Xu(:,term-1);
            Rv = Xv(2:end-1,term-1);
        else 
        end
    
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
            Vr1 = -Pphi'*Rphi;
            Vr2 = B4'*Ru-Pu'*Ru;
            Vr3 = C3'*Rv-Pv'*Rv ;
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
            Vs3 = C3*Sv-Pv*Sv ;
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
        each_count_value = [each_count_value, counter];   
        % 
        % % Calculate execution time for this iteration
        % iterationTime = endTime - startTime;
        % 
        % cumTime = cumTime + iterationTime;
       
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

    each_count_value = [each_count_value, counter];   
    endTime = cputime;       
    enrichmentTime(number) = endTime - startTime; 

    L2phi(number)=QL2phierror;
    L2U(number)=QL2uerror;
    L2V(number)=QL2verror;
end
each_count_b(:) = each_count_value;
 
N_values=Nvalues'; 

figure('Units', 'normalized', 'Position', [0.2 0.2 0.5 0.5]);;
set(gcf, 'Color', 'white')
handles = semilogy(1:no_enrichment,QL2phi, 'Linewidth', 2);
legend('$N=8$','$N=12$','$N=16$','$N=20$','interpreter','latex','FontSize', 14 )
grid on 
xlim([1 no_enrichment])
xticks(1:1:no_enrichment) 
set(gca,'fontsize',14);
xlabel('$Q$ enrichments','interpreter','latex','fontsize',20)
ylabel('$L^2$ Errors in $\phi$','interpreter','latex','fontsize',20) 

fprintf('L2 errors for Example 5, i=3, using the LS SM PGD, after three enrichments:')
fprintf('\n');
fprintf('\n'); 
fprintf('N =       ');
fprintf('%8d ', Nvalues);
fprintf('\n-----------------------------------------------------\n');
fprintf('L2 error φ   ');
fprintf('%6.2e ', QL2phi(3,:)); 
fprintf('\n');

fprintf('L2 errors for u when N=20 : ')
disp(L2U(end))
fprintf('L2 errors for v when N=20 : ')
disp(L2V(end))
fprintf('\n');
fprintf('%d, %d and %d ADFPA iterations required at the first and second enrichment stages for Example 5, i=3, using the LS SM PGD, with epsilon = %6.2e\n', ...
    each_count_b(1), each_count_b(2), each_count_b(3), epsilon);
fprintf('\n');
 
fprintf('Time taken (in seconds):')
fprintf('\n');
fprintf('N =          ');
fprintf('%8d ', Nvalues);
fprintf('\n');
fprintf('Time       ');
disp( enrichmentTime); 


figure; 
set(gcf, 'Color', 'white')
for i = 1:3
    subplot(3,1,i);
    surf(node, node, Xphi(:,i) * Yphi(:,i)');   
    title(sprintf('Approximation from q = %d', i));
    xlabel('x');
    ylabel('y');
end
