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
close all

global N
syms x y
format short

% true solution
phi=@(x,y) sin(pi*x)*sin(pi*y)+ sin(2*pi*x).*sin(2*pi*y) +(x.^2-1)*(y.^2-1); 
U=@(x,y)-pi*cos(pi*x).*sin(pi*y)  - 2*pi*cos(2*pi*x).*sin(2*pi*y) -(2)*x.*(y.^2-1); 
V=@(x,y)- pi*sin(pi*x).*cos(pi*y)   - 2*pi*sin(2*pi*x).*cos(2*pi*y) -(2)*y.*(x.^2-1);
     
% N GLL points
Nvalues=[8, 12, 16, 20] ;

% Number of enrichments
num_enrichments = 5;

% fminunc options
options = optimoptions('fminunc', ...
    'Algorithm', 'quasi-newton', ...  % Use quasi-newton since no gradient
    'MaxIterations', 1e20, ...
    'MaxFunctionEvaluations', 1e20, ...
    'FunctionTolerance', 1e-12, ...
    'StepTolerance', 1e-12, ...
    'OptimalityTolerance', 1e-12, ...
    'FiniteDifferenceType', 'central', ... 
    'FiniteDifferenceStepSize', 1e-07, ... 
    'UseParallel', false);

for number=1:length(Nvalues)
    
    N=Nvalues(number);
    N1=N+1; 
    approxphi = zeros(N+1) ;
    approxu = zeros(N+1) ;
    approxv = zeros(N+1) ;

    % Compute GLL nodes
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
     
    node_transpose = node';
    xs=node;
    ys=node_transpose;

%============= Differentiation matrix and weights ==================

    X=repmat(xx,1,N1); %Creates a square matrix with column entries x
    Xdiff=X-X'+eye(N1); %eye is for Kronecker Delta
        
    LL=repmat(P(:,N1),1,N1); %Replicates Legendre Vandemonde Matrix
    LL(1:(N1+1):N1*N1)=1;
    D=(LL./(Xdiff.*LL'));
    D(1:(N1+1):N1*N1)=0;
    D(1)=-(N1*NN)/4;
    D(N1*N1)=(N1*NN)/4; %Gives differentiation matrix
    w=2./(N*N1*(LL(:,N1).^2)); %Gives weights

%=======================================================================
% Set up matrices
%=======================================================================   

    node_transpose = node';
    % source term
    F=2*pi^2*sin(pi*node)*sin(pi*node_transpose)+ 8*pi^2*sin(2*pi*node)*sin(2*pi*node_transpose)  -2*(node.^2)-2*(node_transpose.^2)+4;

    Xu_history=[];Yu_history=[];
    Xv_history=[];Yv_history=[];
    Xphi_history=[];Yphi_history=[];
    DXu_history=[];DYu_history=[];
    DXv_history=[];DYv_history=[];
    DXphi_history=[];DYphi_history=[];

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
     
    % f = F(:)' * F(:);
    sumn = 0;
    for r=1:N1
        summ=0;
        for m=1:N1
            summ = summ + (w(m)*(F(m,r)^2));
        end
        sumn = sumn + summ*w(r);
    end
    f=sumn;
    
    Bf1 =zeros(N1); Cf1=zeros(N1); Bf2= zeros(N1); Df2=zeros(N1); Df3=zeros(N1); Cf3=zeros(N1); Bf4=zeros(N1); Cf4=zeros(N1);
        
    % Storage for DXu, Yu, Xv, DYv for each iteration
    DXu_history = [];
    Yu_history = [];
    Xv_history = [];
    DYv_history = [];
    f_1=f;
    f_2=0;
    f_3=0;
    f_4=0;

    %Now the B,C,D matrices that depend on specific fs 
    for i=1:N1
        for j=1:N1
        sumBf1=0; 
            for n=1:N1 
            sumBf1=sumBf1+(w(n)*w(j)*D(n,i)*F(n,j)); 
            end 
            Bf1(i,j)=sumBf1; 
        end
    end
  
    for i=1:N1
        for j=1:N1
            sumCf1=0; 
            for n=1:N1
             sumCf1=sumCf1+w(i)*w(n)*F(i,n)*D(n,j); 
            end
            Cf1(i,j)=sumCf1; 
        end
    end
  
    startTime(number) = cputime;
    for enrich = 1:num_enrichments
        X0 =   randn(6*N1,1);  % Initial guess
        
        X0(N1+1) = 0;
        X0(2*N1) = 0;
        X0(2*N1+1) = 0;
        X0(3*N1) = 0; 
        X0(4*N1+1) = 0;
        X0(5*N1) = 0;
        X0(5*N1+1) = 0;
        X0(6*N1) = 0;
            
        % call fminunc
        [X_opt, J_min] = fminunc(@(X) objective_function(X, A, M, L, Bf1, Cf1, Bf2, Df2, Df3, Cf3, Bf4, Cf4, f_1,f_2,f_3,f_4,N1), X0, options);

        % bcs
        X_opt(N1+1)=0; X_opt(2*N1)=0; 
        X_opt(2*N1+1)=0; X_opt(3*N1)=0;
        X_opt(4*N1+1) =0;X_opt(5*N1)=0;
        X_opt(5*N1+1)=0; X_opt(6*N1)=0;
        
        disp('Minimum Objective Function Value (J_min):');
        disp(J_min);

        % Enrichment step: extract vectors from X_opt
        Ru = X_opt(1:N1);   Su = X_opt(N1+1:2*N1);
        Rv = X_opt(2*N1+1:3*N1); Sv = X_opt(3*N1+1:4*N1);
        Rphi = X_opt(4*N1+1:5*N1); Sphi = X_opt(5*N1+1:6*N1);
        Su(1) = 0;  Su(end) = 0;
        Rv(1) = 0; Rv(end) = 0;
        Sphi(1) = 0;  Sphi(end) = 0;
        Rphi(1) = 0; Rphi(end) = 0;
     
        Xu = Ru;    Xv  = Rv;    Xphi = Rphi;
        DXu = D*Ru; DXv = D*Rv;  DXphi = D*Rphi;
        Yu  = Su;   Yv = Sv;     Yphi = Sphi;
        DYu = D*Su; DYv = D*Sv;  DYphi = D*Sphi;
    
        % Store vectors  
        Xu_history(:,enrich) = Xu;          Yu_history(:,enrich) = Yu;
        Xv_history(:,enrich) = Xv;          Yv_history(:,enrich) = Yv;
        Xphi_history(:,enrich) = Xphi;      Yphi_history(:,enrich) = Yphi;
        DXu_history(:,enrich) = DXu;        DYu_history(:,enrich) = DYu;
        DXv_history(:,enrich) = DXv;        DYv_history(:,enrich) = DYv;
        DXphi_history(:,enrich) = DXphi;    DYphi_history(:,enrich) = DYphi;
        
        % Update f1, f2, f3, f4 based on the previous enrichment
        enf1 = zeros(N1); enf2 = zeros(N1); enf3 = zeros(N1); enf4 = zeros(N1);   
     
        for x = 1:N1
            for y = 1:N1
                sumf1 = 0; sumf2 = 0; sumf3 = 0; sumf4 = 0;
                for q = 1:enrich
                    sumf1 = sumf1 + DXu_history(x,q) * Yu_history(y,q) + Xv_history(x,q) * DYv_history(y,q);
                    sumf2 = sumf2 + Xu_history(x,q) * Yu_history(y,q) + DXphi_history(x,q) * Yphi_history(y,q);
                    sumf3 = sumf3 + Xv_history(x,q) * Yv_history(y,q) + Xphi_history(x,q) * DYphi_history(y,q);
                    sumf4 = sumf4 + DXv_history(x,q) * Yv_history(y,q) - Xu_history(x,q) * DYu_history(y,q);
                end
                enf1(x,y) = sumf1;
                enf2(x,y) = sumf2;
                enf3(x,y) = sumf3;
                enf4(x,y) = sumf4;
            end
        end
        ff1 = F - enf1; 
    
        sumF1=0; sumF2=0; sumF3=0; sumF4=0;
        for m=1:N1
            prodf1=0;prodf2=0; prodf3=0; prodf4=0;
            for p=1:N1
                prodf1 = prodf1 + (w(p)*((ff1(m,p))^2));
                prodf2 = prodf2 + (w(p)*(enf2(m,p)^2));
                prodf3 = prodf3 + (w(p)*(enf3(m,p)^2));
                prodf4 = prodf4 + (w(p)*(enf4(m,p)^2));
            end
            sumF1 = sumF1 + prodf1*w(m);
            sumF2 = sumF2 + prodf2*w(m);
            sumF3 = sumF3 + prodf3*w(m);
            sumF4 = sumF4 + prodf4*w(m);
        end
        f_1=sumF1;
        f_2=sumF2 ;
        f_3=sumF3 ;
        f_4=sumF4 ;

        % Now the B,C,D matrices that depend on specific fs 
        Bf2 = zeros(N1); Bf1 = zeros(N1);Bf4 = zeros(N1); 
        for i=1:N1
            for j=1:N1
                sumBf2=0; sumBf1=0; sumBf4=0;
                for m=1:N1
                    sumBf2=sumBf2+(w(m)*w(j)*D(m,i)*enf2(m,j));
                    sumBf1=sumBf1+(w(m)*w(j)*D(m,i)*ff1(m,j));
                    sumBf4=sumBf4+(w(m)*w(j)*D(m,i)*enf4(m,j));
                end
                Bf2(i,j)=sumBf2;
                Bf1(i,j)=sumBf1;
                Bf4(i,j)=sumBf4;
            end
        end

        Cf1=zeros(N1); Cf3=zeros(N1); Cf4=zeros(N1); 
        for i=1:N1
            for j=1:N1
                sumCf1=0; sumCf3=0; sumCf4=0;
                for m=1:N1
                    sumCf1=sumCf1+w(i)*w(m)*ff1(i,m)*D(m,j);
                    sumCf3=sumCf3+w(i)*w(m)*enf3(i,m)*D(m,j);
                    sumCf4=sumCf4+w(i)*w(m)*enf4(i,m)*D(m,j);
                end
                Cf1(i,j)=sumCf1;
                Cf3(i,j)=sumCf3;
                Cf4(i,j)=sumCf4;
            end
        end
     
        Df2=zeros(N1); Df3=zeros(N1); 
        for i=1:N1
            for j=1:N1
                Df2(i,j)=w(i)*w(j)*enf2(i,j);
                Df3(i,j)= w(i)*w(j)*enf3(i,j); 
            end
        end 

        % save approximations
        XU = X_opt(1:N1);  YU = X_opt(N1+1:2*N1);
        XV = X_opt(2*N1+1:3*N1);  YV = X_opt(3*N1+1:4*N1);
        XPHI = X_opt(4*N1+1:5*N1); YPHI = X_opt(5*N1+1:6*N1);
        YU(1) = 0;  YU(end) = 0;
        XV(1) = 0; XV(end) = 0;
        XPHI(1) = 0;  XPHI(end) = 0;
        YPHI(1) = 0; YPHI(end) = 0; 
        Uerror = zeros(N+1);
        Verror = zeros(N+1);
        PHIerror = zeros(N+1);
        PHImesh = zeros(N+1);
        Umesh = zeros(N+1);
        Vmesh = zeros(N+1);
        for k = 1 :N+1
            for l = 1 :N+1
                for j = 1 
                    approxphi(k,l) = approxphi(k,l)+YPHI(l,j)*XPHI(k,j);
                    approxu ( k , l ) = approxu ( k , l ) +YU( l , j ) *XU( k , j ) ;
                    approxv ( k , l ) = approxv ( k , l ) + YV( l , j ) *XV( k , j ) ;
                end
            end
        end
        for k = 1 :N+1
            for l = 1 :N+1
                Uerror(k,l) = approxu(k,l) - U(node(k),node(l));
                Verror(k,l) = approxv(k,l) - V(node(k),node(l));
                PHIerror(k,l) = approxphi(k,l) - phi(node(k),node(l));
                PHImesh(k,l) = phi(node(k),node(l));
                Umesh(k,l) = U(node(k),node(l));
                Vmesh(k,l) = V(node(k),node(l));
            end
        end

        L2phierror=0;
        L2Uerror=0;
        L2Verror=0;
        for num=1:N1
            for m=1:N1
               L2Uerror=L2Uerror+((w(m)*w(num)*(approxu(m,num)-Umesh(m,num))^2)) ;
               L2Verror=L2Verror+((w(m)*w(num)*(approxv(m,num)-Vmesh(m,num))^2)) ;
               L2phierror=L2phierror+((w(m)*w(num)*(approxphi(m,num)-PHImesh(m,num))^2)) ;
            end
        end
        
        L2Uerror= L2Uerror^(1/2);
        L2Verror= L2Verror^(1/2);
        L2phierror= L2phierror^(1/2);
    
        L2phi(enrich,number) = L2phierror;
    end
    endTime(number) = cputime;
    enrichmentTime(number) = endTime(number)-startTime(number);
end


figure; 
set(gca,'fontsize',14);
handles = semilogy(1:num_enrichments,L2phi, 'Linewidth', 2);
set(gcf, 'Color', 'white')
legend('$N=8$','$N=12$','$N=16$','$N=20$' ,'interpreter','latex')
xlabel('$Q$ enrichments','interpreter','latex','fontsize',14)
xlim([1 num_enrichments])
xticks(1:1:num_enrichments)
ylabel('$L^{2}$ Error in $\phi$','interpreter','latex','fontsize',14)
set(gca, 'XGrid', 'on', 'YGrid', 'on');

fprintf('L2 errors for Example 11, i = 3, using the MLS SM PGD, after three enrichments:')
fprintf('\n');
fprintf('\n');
fprintf('N =       ');
fprintf('%8d   ', Nvalues);
fprintf('\n---------------------------------------------------------\n');
fprintf('L2 error φ   ');
fprintf('%6.2e   ', L2phi(end,:));
 
fprintf('\n');
fprintf('\n');

fprintf('Time taken (in seconds):')
fprintf('\n');
fprintf('N =          ');
fprintf('%8d ', Nvalues);
fprintf('\n');
fprintf('Time       ');
disp( enrichmentTime); 
function J = objective_function(X, A, M, L, Bf1, Cf1, Bf2, Df2, Df3, Cf3, Bf4, Cf4, f_1,f_2,f_3,f_4,N1)
    % Extract components from X  
    Ru = X(1:N1);  Su = X(N1+1:2*N1);
    Rv = X(2*N1+1:3*N1);  Sv = X(3*N1+1:4*N1);
    Rphi = X(4*N1+1:5*N1); Sphi = X(5*N1+1:6*N1);
    
    Su(1) = 0;  Su(end) = 0;
    Rv(1) = 0; Rv(end) = 0;
    Sphi(1) = 0;  Sphi(end) = 0;
    Rphi(1) = 0; Rphi(end) = 0;
     
    % Compute J1, J2, J3, J4
    J1 =  (Ru'*A*Ru)*(Su'*M*Su) + (Rv'*M*Rv)*(Sv'*A*Sv) + 2*(Rv'*L*Ru)*(Su'*L*Sv) - 2*(Ru'*Bf1*Su) - 2*(Rv'*Cf1*Sv) + f_1;
    J2 = (Ru'*M*Ru)*(Su'*M*Su) + (Rphi'*A*Rphi)*(Sphi'*M*Sphi) + 2*(Ru'*L*Rphi)*(Su'*M*Sphi) + f_2 + 2*(Ru'*Df2*Su) + 2*(Rphi'*Bf2*Sphi);
    J3 = (Rv'*M*Rv)*(Sv'*M*Sv) + (Rphi'*M*Rphi)*(Sphi'*A*Sphi) + 2*(Rv'*M*Rphi)*(Sv'*L*Sphi) + f_3 + 2*(Rv'*Df3*Sv) + 2*(Rphi'*Cf3*Sphi);
    J4 = (Rv'*A*Rv)*(Sv'*M*Sv) + (Ru'*M*Ru)*(Su'*A*Su) - 2*(Ru'*L*Rv)*(Sv'*L*Su) + f_4 + 2*(Rv'*Bf4*Sv) - 2*(Ru'*Cf4*Su);

    J = J1 + J2 + J3 + J4; J=abs(J);
    if rand() < 0.001  % don't print all updates
        fprintf('J = %e, J1 = %e, J2 = %e, J3 = %e, J4 = %e\n', J, J1, J2, J3, J4);
    end

end

  
  