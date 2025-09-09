close all
clear all
clc

%domain dimensions
Lx = 1;
Ly = 1;

%PGD parameters
epsilon = 1e-20; %this is the tolerence
%PGD enrichment tolerance set to zero to enforce the number of enrichments
epsilon_tilde = 0;
%Maximum number of enrichments
Max_terms = 10; 
%Maximum number of iterations in the fixed point loop
stopcrit=1000;
   
M =[8 12 16 20];
%M=[20]
for ii=1:numel(M)
  
    Nx = M(ii);
    Ny = M(ii);
    
    % Truncation + 1
    NN=Nx-1;
    N1=Nx;
 
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
x=xx;
y=xx;


     U_exact=sin(pi*(1-x.^2).*(1-y'.^2));
     F=4.*pi^2.*(x.^2.*(1-y'.^2).^2+y'.^2.*(1-x.^2).^2).*sin(pi*(1-x.^2).*(1-y'.^2))+2.*pi*((1-x.^2)+(1-y'.^2)).*cos(pi.*(1-x.^2).*(1-y'.^2));      

    X=repmat(xx,1,N1); %Creates a square matrix with column entries x
    Xdiff=X-X'+eye(N1); %eye is for Kronecker Delta
    L=repmat(P(:,N1),1,N1); %Replicates Legendre Vandemonde Matrix
    L(1:(N1+1):N1*N1)=1;
    D=(L./(Xdiff.*L'));
    D(1:(N1+1):N1*N1)=0;
    D(1)=-(N1*NN)/4;
    D(N1*N1)=(N1*NN)/4; %Gives differentiation matrix

    w=2./(NN*N1*(L(:,N1).^2)); %Computation of weights
    n=1:N1;
    for i=1:N1
        for j=1:N1
        A(i,j)=sum(w(n).*D(n,i).*D(n,j));
        end
    end 
    
   Ay=A;
   Ax=A;
   MM=diag(w(n));
   Mx=MM;
   My=MM;   
X=[];
Y=[];
    N0 = size(X,2);%size(X,dim), this function will return the size of X’s dimension, specified by the input scalar dim.
           cumTime=0;
    %main enrichment loop
    for term=1:(Max_terms)
        startTime = cputime;        
        %make an initial guess for a and b
        a = randn(Nx,1); %returns a Nx by 1 matrix (column)
        b = randn(Ny,1);
        
        %Homogeneous Dirichlet Boundary conditions for a in 0 and Lx
        %And for b in 0 and Ly 
        a(1) = 0;
        a(end)=0;
        b(1) = 0;
        b(end)=0;
     
        counter=0;
        a_old= zeros((Nx),1);
        b_old= zeros((Nx),1);
        
        %fixed point iterations
        while max(max( abs (a*b' - a_old* b_old') ) ) >epsilon
       
            %Save the old values of a & b for later comparison
            a_old = a;
            b_old = b;
            
            %Solve for a
            %LHS is:(b'M_yb)(A_xa)+(b'A_yb)(M_xa)
            %LHS coefficients
            bMyb = b'*My*b;
            bAyb = b'*Ay*b;
            
            %Construction of the RHS         
            for k=1:N1
                summ=0;
                     for i=1:N1
                        summ=summ + b(i)*w(i)*F(i,k);
                     end
                fxy_a(k)=w(k)*summ;
             end
  
  
            RHS=fxy_a';
   
            %In case this is not the first enrichment, previously computed
            %terms are added to the RHS
            if (term>1)
                %RHS coefficients
                AxX=Ax*X;
                bMyY=b'*My*Y;
                MxX=Mx*X;
                bAyY=b'*Ay*Y;
               
                RHS=RHS-(AxX*bMyY'+MxX*bAyY');
            end
            
            %Construction of the FE problem
            LHS =( bMyb*Ax + bAyb*Mx);
            
            %Solution with homogeneous boundary conditions at both ends
            a(2:end-1) = LHS(2:end-1,2:end-1)\RHS(2:end-1); %Since Aa=RHS
            
            %Solve for b
            %LHS is:(a'A_xa)(M_yb)+(a'M_xa)(A_yb)
            %LHS coefficients
            aMxa = a'*Mx*a;
            aAxa = a'*Ax*a;
            
            %Construction of the RHS
             for k=1:N1
                summ=0;
                     for i=1:N1
                        summ=summ + a(i)*w(i)*F(i,k);
                     end
                fxy_b(k)=w(k)*summ;
             end
        
            RHS=fxy_b';
            %In case this is not the first enrichment, previously computed
            %terms are added to the RHS
            if (term>1)
                %RHS coefficients
                AyY=Ay*Y;
                aMxX=a'*Mx*X;
                MyY=My*Y;
                aAxX=a'*Ax*X;
                
                RHS=RHS-(AyY*aMxX'+MyY*aAxX');
            end
            
            %construction of the FE problem
            LHS=(aMxa*Ay+aAxa*My);
            %solution with homogeneous boundary conditions at both ends
            b(2:end-1) = LHS(2:end-1,2:end-1)\RHS(2:end-1); %Since Ab=RHS
            
            %Norm of the difference between the 2 fixed point iterations
%             S_difference = sqrt(aMxa*bMyb + (a_old'*Mx*a_old)*(b_old'*My*b_old) - 2*(a'*Mx*a_old)*(b'*My*b_old));
%             %fixed point exit test
%             if(S_difference < epsilon), break; end
             counter=counter+1;
      
            if counter == stopcrit
                break
            end
            
        end
          
        term
        counter
        max(max( abs (a*b' - a_old* b_old') ) )
        
        %New normalized enrichment is added to existing
        fact_x = sqrt((aMxa)/(x(end)-x(1))^2);
        fact_y = sqrt((bMyb)/(y(end)-y(1))^2);
        fact_xy = sqrt(fact_x*fact_y);
        X = [X fact_xy*a/fact_x];
        Y = [Y fact_xy*b/fact_y];
        
        %simplified stopping criterion
        E = sqrt((aMxa)*(bMyb)) / sqrt((X(:,N0+1)'*Mx*X(:,N0+1))*(Y(:,N0+1)'*My*Y(:,N0+1)));
        
        if(E<epsilon_tilde), break; end;
        % Record end time
    endTime = cputime;       
            
    % Calculate execution time for this iteration
    iterationTime = endTime - startTime;
    
    cumTime = cumTime + iterationTime;
    % Append to the vector, cumulatively
    executionTimes(term,ii) = cumTime;
     
    end
   

%Error computation
    for j=1:Max_terms
        U_pgd = Y(:,1:j)*X(:,1:j)';
         E_MN(ii,j) = trapz(x,trapz(y,(U_exact - U_pgd).^2,1),2);
         E_MNL2(ii,j) = sqrt(trapz(x,trapz(y,(U_exact - U_pgd).^2,1),2));
    
    L2phierror=0;
    for num=1:N1
        for m=1:N1
           L2phierror=L2phierror+((w(m)*w(num)*(U_pgd(m,num)-U_exact(m,num))^2)) ;
        end
    end
    L2phierror= L2phierror^(1/2);
    E_MNL2(ii,j)=L2phierror;
    end
   
end



%Error as a function of N for different values of M
figure;
set(gca,'fontsize',15);
handles = semilogy(1:Max_terms,E_MNL2, 'linewidth', 2);
grid on
%legend(handles,cellfun(@(in) ['N=' num2str(in)],num2cell(M),'uniformoutput',false))
legend('$N=8$','$N=12$','$N=16$','$N=20$','interpreter','latex')
set(gca,'xtick',[ 1 2 3 4 5 6 7 8 9 10]);
set(gca,'xlim',[1 Max_terms]);
xlabel('$Q$ enrichments','interpreter','latex', 'fontsize',15)
ylabel('$L^{2}$ Error','interpreter','latex','fontsize',15)

% Display the execution times
disp('Execution times for each iteration:');
disp(executionTimes);

%Error for cpu time
figure;
set(gca,'fontsize',14);
handles=plot(1:Max_terms, executionTimes)
set(gca,'xtick',[0 1 2 3 4 5 6 7 8 9 10]);
set(gca,'xlim',[1 Max_terms]);
%legend(handles,cellfun(@(in) ['N=' num2str(in)],num2cell(M),'uniformoutput',false))
legend('$N=8$','$N=8$','$N=8$','$N=8$','interpreter','latex')
xlabel('$Q$ enrichments','interpreter','latex','fontsize',15)
ylabel('CPU time','interpreter','latex','fontsize',15)