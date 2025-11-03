% =====================================================================
% The scripts called in for computation of nodes, differentiation
% matrix, and weights are based on:
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
format short
global N
syms x y
format short
 
%Number of enrichments
no_enrichment=1;

%Max iterations
maxiterations = 100;
 
% number of GLL points
Nvalues=[12,16,20,24]; 

 % tolerance
 epsilon = 10^(-8);


cumTime =zeros(length(Nvalues),1);
executionTimes = zeros(no_enrichment,length(Nvalues));

for number=1:length(Nvalues)
    cumcounter =0;
    cumcounter_values =[]; 
    each_count_value=[];

    N=Nvalues(number);
    NP1=N+1;
  
    mu=1; 
    node=Nodes;
    w=Weights;
    D=DMatrix;
    Dtilde=DtildeMatrix;
    E=Ematrix; 
    
    nodex=node;
    nodey=node;

    UUU = zeros(size(nodex));
    VVV = zeros(size(nodex));
    PPP = zeros(size(nodex));
    VV1 = zeros(size(nodex));
    VV2 = zeros(size(nodex));
    VV3 = zeros(size(nodex));
    VV4 = zeros(size(nodex));
    
    %Analytical solutions
    for xx=1:numel(nodex)
        for yy=1:numel(nodey)
            UUU(xx,yy) = -sin(2*pi*nodey(yy))*(cos(2*pi*nodex(xx))-1);
            VVV(xx,yy) = sin(2*pi*nodex(xx))*(cos(2*pi*nodey(yy))-1);
            PPP(xx,yy) = sin(pi*nodex(xx))*sin(pi*nodey(yy));
            VV1(xx,yy) = 2*pi*sin(2*pi*nodex(xx))*sin(2*pi*nodey(yy));
            VV2(xx,yy) = -2*pi*(cos(2*pi*nodex(xx))-1)*cos(2*pi*nodey(yy));
            VV3(xx,yy) = 2*pi*cos(2*pi*nodex(xx))*(cos(2*pi*nodey(yy))-1);
            VV4(xx,yy) = -2*pi*sin(2*pi*nodex(xx))*sin(2*pi*nodey(yy));
         end
    end
    Umesh=UUU;
    Vmesh=VVV;
    V1mesh=VV1;
    V2mesh=VV2;
    V3mesh=VV3;
    V4mesh=VV4;
    Pmesh=PPP;
    xx=node;
    for x=1:numel(xx)
        for y=1:numel(xx)
            F1(x,y) = pi*cos(pi*xx(x))*sin(pi*xx(y))-(4*pi^2*sin(2*pi*xx(y)))*(2*cos(2*pi*xx(x))-1) ;
            F2(x,y) = pi*cos(pi*xx(y))*sin(pi*xx(x))+(4*pi^2*sin(2*pi*xx(x)))*(2*cos(2*pi*xx(y))-1) ;
        end
    end
         
    A = StiffnessMatrix;
    A_2 = A(2:N,2:N);
    A_3 = A(2:N,:);
    A_4 = A_3';

    M = MassMatrix;
    M_2 = M(2:N,2:N);
    M_3 = M(2:N,:);
    M_4 = M_3';

    L = MixMatrix';
    L_2 = L(2:N,2:N);
    L_3 = L(2:N,:);
    L_4 = L_3';
    Zrow = zeros(1,N+1);
    
    At = StiffnessMatrixt;
    At_2 = At(2:N,:);
    At_4 = At;
    At_3=At_4';
    Zcol2= zeros(N+1,1);
    At_5 = [Zcol2, At, Zcol2];
    At = [Zrow; At_3; Zrow];

    Att = StiffnessMatrixtt;
    Att_2 = Att;
    Zcol = zeros(N-1,1);
    Att_3 = [Zcol, Att, Zcol];
    Att = [Zrow; Att_3; Zrow];

    Mt = MassMatrixt;
    Mt_2 = Mt(2:N,:);
    Mt_4 = Mt;
    Mt_3 = Mt_4';
    Mt = [Zrow; Mt_3; Zrow];

    Mtt = MassMatrixtt;
    Mtt_2 = Mtt;
    Mtt_3 = [Zcol, Mtt, Zcol ]; 
    Mtt = [Zrow; Mtt_3; Zrow];

    Lt = MixMatrixt;
    Lt_2 =  Lt(2:N,:);
    Lt_3 = Lt';
    Lt = [Zrow; Lt_3; Zrow];
   

    G=Gmatrix;
    G_2=G(:,2:N);
    G_3=G; 
    G_4 = G_3';
    G = [Zrow; G; Zrow];

    %Start off Q enrichments
    Q=0;

    %Start approximations
    Xu = zeros(N+1,Q);
    Xv = zeros(N+1,Q);
    Xp = zeros(N+1,Q);
    XV1 = zeros(N+1,Q);
    XV2 = zeros(N+1,Q);
    XV3 = zeros(N+1,Q);
    XV4 = zeros(N+1,Q);
    Yu = zeros(N+1,Q);
    Yv = zeros(N+1,Q);
    Yp = zeros(N+1,Q);
    YV1 = zeros(N+1,Q);
    YV2 = zeros(N+1,Q);
    YV3 = zeros(N+1,Q);
    YV4 = zeros(N+1,Q);
    DRu = zeros(N+1,Q);
    DRv = zeros(N+1,Q); 
    DRV1 = zeros(N+1,Q);
    DRV2 = zeros(N+1,Q);
    DRV3 = zeros(N+1,Q);
    DRV4 = zeros(N+1,Q);
    DSu = zeros(N+1,Q);
    DSv = zeros(N+1,Q); 
    DSV1 = zeros(N+1,Q);
    DSV2 = zeros(N+1,Q);
    DSV3 = zeros(N+1,Q);
    DSV4 = zeros(N+1,Q);
    
    counter=0; 

    for term = 1:no_enrichment

        for m =1:NP1
            for q= 1:Q
                sumDRu = 0;
                sumDRv = 0;
                sumDRV1 = 0;
                sumDRV2 = 0;
                sumDRV3 = 0;
                sumDRV4 = 0;
                sumDSu = 0;
                sumDSv = 0;
                sumDSV1 = 0;
                sumDSV2 = 0;
                sumDSV3 = 0;
                sumDSV4 = 0;
                for k=1:NP1
                    sumDRu = sumDRu + D(m,k)*Xu(k,q);
                    sumDRv = sumDRv + D(m,k)*Xv(k,q);
                    sumDRV1 = sumDRV1 + D(m,k)*XV1(k,q);
                    sumDRV2 = sumDRV2 + D(m,k)*XV2(k,q);
                    sumDRV3 = sumDRV3 + D(m,k)*XV3(k,q);
                    sumDRV4 = sumDRV4 + D(m,k)*XV4(k,q);
                    sumDSu = sumDSu + D(m,k)*Yu(k,q);
                    sumDSv = sumDSv + D(m,k)*Yv(k,q);
                    sumDSV1 = sumDSV1 + D(m,k)*YV1(k,q);
                    sumDSV2 = sumDSV2 + D(m,k)*YV2(k,q);
                    sumDSV3 = sumDSV3 + D(m,k)*YV3(k,q);
                    sumDSV4 = sumDSV4 + D(m,k)*YV4(k,q);
                end
                DRu(m,q) = sumDRu;
                DRv(m,q) = sumDRv;
                DRV1(m,q) = sumDRV1;
                DRV2(m,q) = sumDRV2;
                DRV3(m,q) = sumDRV3;
                DRV4(m,q) = sumDRV4;
                DSu(m,q) = sumDSu;
                DSv(m,q) = sumDSv;
                DSV1(m,q) = sumDSV1;
                DSV2(m,q) = sumDSV2;
                DSV3(m,q) = sumDSV3;
                DSV4(m,q) = sumDSV4;
            end
        end
        
        for m =1:NP1 
            for q= 1:Q
            sumDRp = 0;
            sumDSp = 0;
            for k=1:N-1 
                    sumDRp = sumDRp + Dtilde(m,k+1)*Xp(k+1,q);
                    sumDSp = sumDSp + Dtilde(m,k+1)*Yp(k+1,q);
            end
                    DRp(m,q) = sumDRp;
                    DSp(m,q) = sumDSp;
            end
        end
 

        Vr1 = zeros(N-1,1);
        Vr2 = zeros(N-1,1);
        Vr3 = zeros(N-1,1);
        Vr4 = zeros(N-1,1);
        Vr5 = zeros(N+1,1);
        Vr6 = zeros(N-1,1);
        Vr7 = zeros(N+1,1);
        Vs1 = zeros(N-1,1);
        Vs2 = zeros(N-1,1);
        Vs3 = zeros(N-1,1);
        Vs4 = zeros(N+1,1);
        Vs5 = zeros(N-1,1);
        Vs6 = zeros(N+1,1);
        Vs7 = zeros(N-1,1);

        Vr1a = zeros(N-1,1);
        Vr1b = zeros(N-1,1);
        Vr1c = zeros(N-1,1);
        Vr2a = zeros(N-1,1);
        Vr2b = zeros(N-1,1);
        Vr2c = zeros(N-1,1);
        Vr3a = zeros(N-1,1);
        Vr3b = zeros(N-1,1);
        Vr4a = zeros(N-1,1);
        Vr4b = zeros(N-1,1);
        Vr4c = zeros(N-1,1);
        Vr4d = zeros(N-1,1);
        Vr4e = zeros(N-1,1);
        Vr5a = zeros(N+1,1);
        Vr5b = zeros(N+1,1);
        Vr5c = zeros(N+1,1);
        Vr6a = zeros(N-1,1);
        Vr6b = zeros(N-1,1);
        Vr6c = zeros(N-1,1);
        Vr7a = zeros(N+1,1);
        Vr7b = zeros(N+1,1);
        Vr7c = zeros(N+1,1);
        Vr7d = zeros(N+1,1);
        Vr7e = zeros(N+1,1);

        Vs1a = zeros(N-1,1);
        Vs1b = zeros(N-1,1);
        Vs1c = zeros(N-1,1);
        Vs2a = zeros(N-1,1);
        Vs2b = zeros(N-1,1);
        Vs2c = zeros(N-1,1);
        Vs3a = zeros(N-1,1);
        Vs3b = zeros(N-1,1);
        Vs4a = zeros(N+1,1);
        Vs4b = zeros(N+1,1);
        Vs4c = zeros(N+1,1);
        Vs4d = zeros(N+1,1);
        Vs4e = zeros(N+1,1);
        Vs5a = zeros(N-1,1);
        Vs5b = zeros(N-1,1);
        Vs5c = zeros(N-1,1);
        Vs6a = zeros(N+1,1);
        Vs6b = zeros(N+1,1);
        Vs6c = zeros(N+1,1);
        Vs7a = zeros(N-1,1);
        Vs7b = zeros(N-1,1);
        Vs7c = zeros(N-1,1);
        Vs7d = zeros(N-1,1);
        Vs7e = zeros(N-1,1);

        %Unknown functions for enrichment
        R1 = ones( 7*N-3 ,1) ;  
        S1 = ones( 7*N-3 ,1) ;
        oldR = zeros( 7*N-3 ,1) ;
        oldS = zeros( 7*N-3 ,1) ;
            
        Su = rand(N-1 ,1);
        Sv = rand(N-1 ,1);
        Sp = rand(N-1,1);
        SV1 = rand(N-1,1);
        SV2 = rand(N+1,1);
        SV3 = rand(N-1,1);
        SV4 = rand(N+1,1);
      
        counter=0;
    
        startTime = cputime;
        % iteration Stage
        while max(max( abs (R1*S1' - oldR* oldS') ) ) > epsilon

            oldR = R1 ;
            oldS = S1;
            
            %Entries for the first Matrix system
            Mr11 = (Su'*M_2*Su)*A_2+(Su'*M_2*Su)*A_2+(Su'*A_2*Su)*M_2;
            Mr12 = (Su'*L_2*Sv)*L_2';
            Mr21 = (Sv'*L_2'*Su)*L_2;   
            Mr13 = zeros(N-1,N-1);            
            Mr31 = zeros(N-1,N-1);
            Mr14 = -(Su'*M_2*SV1)*L_4';
            Mr41 = -(SV1'*M_2*Su)*L_4;
            Mr15 = -(Su'*L_3*SV2)*M_2; 
            Mr51 = -(SV2'*L_3'*Su)*M_2;
            Mr16 = zeros(N-1,N+1);
            Mr61 = zeros(N+1,N-1);
            Mr17 = zeros(N-1,N-1);
            Mr71 = zeros(N-1,N-1);

            Mr22 = (Sv'*A_2*Sv)*M_2+(Sv'*M_2*Sv)*A_2+(Sv'*A_2*Sv)*M_2;
            Mr23 = zeros(N-1,N-1);
            Mr32 = zeros(N-1,N-1);
            Mr24 = zeros(N-1,N+1);
            Mr42 = zeros(N+1,N-1);
            Mr25 = zeros(N-1,N-1);
            Mr52 = zeros(N-1,N-1);
            Mr26 = -(Sv'*M_2*SV3)*L_4';
            Mr62 = -(SV3'*M_2*Sv)*L_4;
            Mr27 = -(Sv'*L_4'*SV4)*M_2;
            Mr72 = -(SV4'*L_4*Sv)*M_2;
              
            H=Hmatrix;
            sum_HSp = 0;
            for j = 1:N-1
                sum_HSp = sum_HSp + H(j)*Sp(j);
            end
            J6_r33= H'*H*(sum_HSp); 
            Mr33 = (Sp'*Mtt_2*Sp)*Att_2+(Sp'*Att_2*Sp)*Mtt_2+ mu*J6_r33;
            Mr34 = -(Sp'*Mt_2'*SV1)*At_4';
            Mr43 = -(SV1'*Mt_2*Sp)*At_4;
            Mr35 = -(Sp'*G_3*SV2)*Lt_2';
            Mr53 = -(SV2'*G_3'*Sp)*Lt_2;
            Mr36 = -(Sp'*Lt_2'*SV3)*G_3;
            Mr63 = -(SV3'*Lt_2*Sp)*G_3';
            Mr37 = -(Sp'*At_4'*SV4)*Mt_2';
            Mr73 =- (SV4'*At_4*Sp)*Mt_2;
            
            Mr44 = (SV1'*M_2*SV1)*A+(SV1'*M_2*SV1)*M+(SV1'*M_2*SV1)*A+(SV1'*A_2*SV1)*M+(SV1'*A_2*SV1)*M;
            Mr45 = zeros(N+1, N-1); 
            Mr54 = zeros(N-1, N+1);
            Mr46 = zeros(N+1,N+1);
            Mr64 = zeros(N+1,N+1);
            Mr47 = (SV1'*M_3*SV4)*A_4+(SV1'*A_3*SV4)*M_4;
            Mr74 = (SV4'*M_4*SV1)*A_3+(SV4'*A_4*SV1)*M_3;

            Mr55 = (SV2'*A*SV2)*M_2+(SV2'*M*SV2)*M_2+(SV2'*M*SV2)*A_2;
            Mr56 = zeros(N-1,N+1);
            Mr65 = zeros(N+1,N-1);
            Mr57 = zeros(N-1,N-1);
            Mr75 = zeros(N-1,N-1);
            
            Mr66 = (SV3'*M_2*SV3)*A+(SV3'*M_2*SV3)*M+(SV3'*A_2*SV3)*M;
            Mr67 = zeros(N+1,N-1);
            Mr76 = zeros(N-1,N+1);
            
            Mr77 = (SV4'*A*SV4)*M_2+(SV4'*M*SV4)*M_2+(SV4'*M*SV4)*A_2+(SV4'*A*SV4)*M_2+(SV4'*M*SV4)*A_2;

            %First matrix
            Matrix1 = [Mr11, Mr12, Mr13, Mr14, Mr15, Mr16, Mr17; Mr21, Mr22, Mr23, Mr24, Mr25, Mr26, Mr27; Mr31, Mr32, Mr33, Mr34, Mr35, Mr36, Mr37; Mr41, Mr42, Mr43, Mr44, Mr45, Mr46, Mr47; Mr51, Mr52, Mr53, Mr54, Mr55, Mr56, Mr57; Mr61, Mr62, Mr63, Mr64, Mr65, Mr66, Mr67; Mr71, Mr72, Mr73, Mr74, Mr75, Mr76, Mr77];
          
            %Entries for first matrix system RHS
            Vs1=0.*Su;
            Vs2=0.*Sv;
     
            for i=1:N-1 
                summ1 = 0;
                summ2 = 0;
                    for n=1:NP1
                        outsum1=0;
                        outsum2=0;
                        insum1=0;
                        insum2=0;
                            for k=1:N-1
                                insum1=insum1+E(n,k)*Sp(k);  
                                insum2=insum2+Dtilde(n,k+1)*Sp(k);
                            end
                            for m=1:NP1
                                outsum1= outsum1+w(m)*Dtilde(m,i+1)*F1(m,n); 
                                outsum2= outsum2+w(m)*E(m,i)*F2(m,n); 
                            end
                        summ1 = summ1 + w(n)*insum1*outsum1;
                        summ2 = summ2 + w(n)*insum2*outsum2;
                    end
                    Vs3(i)=summ1+summ2;                                      
            end
            Vs3=Vs3(:);
            
            for i=1:NP1 
                outsum=0;
                for m=1:NP1
                    insum=0;
                    for n=1:N-1
                        insum=insum+w(n+1)*F1(m,n+1)*SV1(n);
                    end
                    outsum=outsum+insum*w(m)*D(m,i);
                end
                Vs4(i) = -outsum;
            end
            Vs4=Vs4(:);

            for i=1:N-1 
                outsum=0;
                for n=1:NP1
                    insum=0;
                    for k=1:NP1
                        insum=insum+SV2(k)*D(n,k);
                    end
                    outsum=outsum+insum*w(i+1)*w(n)*F1(i+1,n);
                end
                Vs5(i)=-outsum;
            end
            Vs5=Vs5(:);

            for i=1:NP1 
                outsum=0;
                for m=1:NP1
                    insum=0;
                    for n=1:N-1
                        insum=insum+w(n+1)*F2(m,n+1)*SV3(n);
                    end
                    outsum=outsum+insum*w(m)*D(m,i);
                end
                Vs6(i)=-outsum;
            end
            Vs6=Vs6(:);
             
            for i=1:N-1 
                outsum=0;
                for n=1:NP1
                    insum=0;
                    for k=1:NP1
                        insum=insum+SV4(k)*D(n,k);
                    end
                    outsum=outsum+insum*w(n)*w(i+1)*F2(i+1,n);
                end
                Vs7(i)=-outsum;
            end
            Vs7=Vs7(:);

            % NEW ENRICHMENT TERMS
            for i = 1:N-1
                outer_s1a = 0;
                for m=1:N+1
                    inner_s1a = 0;
                    for k =1:N-1
                        Enrich_s1a = 0;
                        for q=1:Q
                            Enrich_s1a = Enrich_s1a + DRu(m,q)*Yu(k+1,q) + Xv(m,q)*DSv(k+1,q);
                        end 
                        inner_s1a = inner_s1a + Enrich_s1a*w(k+1)*Su(k);   
                    end
                    outer_s1a = outer_s1a + inner_s1a*w(m)*D(m,i+1);
                end
                Vs1a(i) = outer_s1a;
            end
             
            for i = 1:N-1
                outer_s1b = 0;
                for m=1:N+1
                    inner_s1b = 0;
                    for k =1:N-1
                        Enrich_s1b = 0;
                        for q=1:Q
                            Enrich_s1b = Enrich_s1b + XV1(m,q)*YV1(k+1,q) - DRu(m,q)*Yu(k+1,q);
                        end           
                        inner_s1b = inner_s1b + Enrich_s1b*w(k+1)*Su(k);
                    end
                    outer_s1b = outer_s1b + inner_s1b*w(m)*D(m,i+1);
                end
                Vs1b(i) = - outer_s1b;
            end
             
            for i = 1:N-1
                outer_s1c = 0;
                for m=1:N+1
                        Enrich_s1c = 0;
                        for q=1:Q
                            Enrich_s1c = Enrich_s1c + XV2(i+1,q)*YV2(m,q) - Xu(i+1,q)*DSu(m,q);
                        end 
                        inner_s1c = 0;
                        for l=1:N-1
                            inner_s1c = inner_s1c + D(m,l+1)*Su(l);           
                        end
                    outer_s1c = outer_s1c + Enrich_s1c*inner_s1c*w(m);
                end
                Vs1c(i) = - outer_s1c*w(i+1);
            end
            
            Vs1a=Vs1a(:);
            Vs1b=Vs1b(:);
            Vs1c=Vs1c(:);
            Vs11= Vs1a+Vs1b+Vs1c;
                      
            for i = 1:N-1
                outer_s2a = 0;
                for m=1:N+1
                        Enrich_s2a = 0;
                        for q=1:Q
                            Enrich_s2a = Enrich_s2a + DRu(i+1,q)*Yu(m,q) + Xv(i+1,q)*DSv(m,q);
                        end
                        inner_s2a = 0;
                        for l=1:N-1
                            inner_s2a = inner_s2a + D(m,l+1)*Sv(l);
                        end
                    outer_s2a = outer_s2a + Enrich_s2a*inner_s2a*w(m);
                end
                Vs2a(i) = outer_s2a*w(i+1); 
            end
            
            for i = 1:N-1
                outer_s2b = 0;
                for m=1:N+1
                    inner_s2b = 0;
                    for k =1:N-1
                        Enrich_s2b = 0;
                        for q=1:Q
                            Enrich_s2b = Enrich_s2b + XV3(m,q)*YV3(k+1,q) - DRv(m,q)*Yv(k+1,q);
                        end
                        inner_s2b = inner_s2b + Enrich_s2b*w(k+1)*Sv(k);
                    end
                    outer_s2b = outer_s2b + inner_s2b*w(m)*D(m,i+1);
                end
                Vs2b(i) = -outer_s2b;
            end
            
            for i = 1:N-1
                outer_s2c = 0;
                for m=1:N+1
                        Enrich_s2c = 0;
                        for q=1:Q
                            Enrich_s2c = Enrich_s2c + XV4(i+1,q)*YV4(m,q) - Xv(i+1,q)*DSv(m,q);
                        end
                        inner_s2c = 0;
                        for l=1:N-1
                            inner_s2c = inner_s2c + D(m,l+1)*Sv(l);
                        end
                    outer_s2c = outer_s2c + Enrich_s2c*inner_s2c*w(m);
                end
                Vs2c(i) = -outer_s2c*w(i+1);
            end
            
            Vs2a=Vs2a(:);
            Vs2b=Vs2b(:);
            Vs2c=Vs2c(:);
            Vs22= Vs2a+Vs2b+Vs2c;
            
            for i = 1:N-1
                outer_s3a = 0;
                for m=1:N+1
                    middle_s3a = 0;
                    for n=1:N+1
                        Enrich_s3a = 0;
                        for q=1:Q
                            Enrich_s3a = Enrich_s3a +(- DRV1(m,q)*YV1(n,q) - XV2(m,q)*DSV2(n,q) + DRp(m,q)*Yp(n,q));
                        end
                        inner_s3a = 0;
                        for l=1:N-1
                            inner_s3a = inner_s3a + Sp(l)*E(n,l);                        
                        end
                        middle_s3a = middle_s3a + Enrich_s3a*inner_s3a*w(n);
                    end
                    outer_s3a = outer_s3a + middle_s3a*w(m)*Dtilde(m,i+1);
                end
                Vs3a(i) = outer_s3a;
            end
            
            for i = 1:N-1
                outer_s3b = 0;
                for m=1:N+1
                    middle_s3b = 0;
                    for n=1:N+1
                        Enrich_s3b = 0;
                        for q=1:Q
                            Enrich_s3b = Enrich_s3b +(- DRV3(m,q)*YV3(n,q) - XV4(m,q)*DSV4(n,q) + Xp(m,q)*DSp(n,q));
                        end
                        inner_s3b = 0;
                        for l=1:N-1          
                            inner_s3b = inner_s3b + Sp(l)*Dtilde(n,l+1);                 
                        end
                        middle_s3b = middle_s3b + Enrich_s3b*inner_s3b*w(n); 
                    end
                    outer_s3b = outer_s3b + middle_s3b*w(m)*E(m,i);
                end
                Vs3b(i) = outer_s3b;
            end
             
            for i = 1:N-1
                middle_s3c = 0;
                for m=1:N+1
                    inner_s3c = 0;
                    for n=1:N+1
                        Enrich_s3c = 0;
                        for q=1:Q
                            Enrich_s3c = Enrich_s3c +(Xp(m,q)*Yp(n,q));
                        end
                        inner_s3c = inner_s3c + Enrich_s3c; 
                    end
                    middle_s3c = middle_s3c + inner_s3c;
                end
                outer_s3c = 0;
                for j=1:N-1          
                   outer_s3c = outer_s3c + Sp(j)*H(j);                 
                end
                Vs3c(i) = middle_s3c*outer_s3c*H(i)*mu;
            end
            
            Vs3a=Vs3a(:);
            Vs3b=Vs3b(:);
            Vs3c=Vs3c(:);
            Vs33= Vs3a+Vs3b+Vs3c;
            
            for i = 1:N+1
                outer_s4a = 0;
                for m=1:N+1
                    inner_s4a = 0;
                    for k =1:N-1
                        Enrich_s4a = 0;
                        for q=1:Q
                            Enrich_s4a = Enrich_s4a + (- DRV1(m,q)*YV1(k+1,q) - XV2(m,q)*DSV2(k+1,q) + DRp(m,q)*Yp(k+1,q));
                        end
                        inner_s4a = inner_s4a + Enrich_s4a*w(k+1)*SV1(k);    
                    end
                    outer_s4a = outer_s4a + inner_s4a*w(m)*D(m,i);
                end
                Vs4a(i) = - outer_s4a;
            end
            
            for i = 1:N+1
                inner_s4b = 0;
                for k =1:N-1
                    Enrich_s4b = 0;
                    for q=1:Q
                        Enrich_s4b = Enrich_s4b + XV1(i,q)*YV1(k+1,q) - DRu(i,q)*Yu(k+1,q);
                    end     
                    inner_s4b = inner_s4b + Enrich_s4b*w(k+1)*SV1(k); 
                end
                Vs4b(i) = inner_s4b*w(i);
            end
            
            for i = 1:N+1
                outer_s4c = 0;
                for m=1:N+1
                    inner_s4c = 0;
                    for k =1:N-1
                        Enrich_s4c = 0;
                        for q=1:Q
                            Enrich_s4c = Enrich_s4c + DRV1(m,q)*YV1(k+1,q) + DRV4(m,q)*YV4(k+1,q);
                        end        
                        inner_s4c = inner_s4c + Enrich_s4c*w(k+1)*SV1(k); 
                    end
                    outer_s4c = outer_s4c + inner_s4c*w(m)*D(m,i);
                end
                Vs4c(i) = outer_s4c;
            end
            
            for i = 1:N+1
                outer_s4d = 0;
                for m=1:N+1
                    for k =1:N-1
                        Enrich_s4d = 0;
                        for q=1:Q
                            Enrich_s4d = Enrich_s4d + XV1(i,q)*DSV1(m,q) + XV4(i,q)*DSV4(m,q);
                        end
                        inner_s4d = 0;
                        for l =1:N-1
                            inner_s4d = inner_s4d + D(m,l+1)*SV1(l);             
                        end
                    end
                    outer_s4d = outer_s4d + Enrich_s4d*inner_s4d*w(m);
                end
                Vs4d(i) = outer_s4d*w(i);
            end
            
            for i = 1:N+1
                outer_s4e = 0;
                for m=1:N+1
                        Enrich_s4e = 0;
                        for q=1:Q
                            Enrich_s4e = Enrich_s4e + DRV2(i,q)*YV2(m,q) - XV1(i,q)*DSV1(m,q);
                        end
                        inner_s4e = 0;
                        for l =1:N-1         
                            inner_s4e = inner_s4e + D(m,l+1)*SV1(l);           
                        end
                    outer_s4e = outer_s4e + Enrich_s4e*inner_s4e*w(m);
                end
                Vs4e(i) = - outer_s4e*w(i);
            end
            Vs4a=Vs4a(:);
            Vs4b=Vs4b(:);
            Vs4c=Vs4c(:);
            Vs4d=Vs4d(:);
            Vs4e=Vs4e(:);
            Vs44= Vs4a+Vs4b+Vs4c+Vs4d+Vs4e;
            
            for i = 1:N-1
                outer_s5a = 0;
                for m=1:N+1
                        Enrich_s5a = 0;
                        for q=1:Q
                            Enrich_s5a = Enrich_s5a - DRV1(i+1,q)*YV1(m,q) - XV2(i+1,q)*DSV2(m,q) + DRp(i+1,q)*Yp(m,q);
                        end
                        inner_s5a = 0;
                        for l =1:N+1
                            inner_s5a = inner_s5a + D(m,l)*SV2(l);         
                        end
                    outer_s5a = outer_s5a + Enrich_s5a*inner_s5a*w(m);
                end
                Vs5a(i) = - outer_s5a*w(i+1);
            end
             
            for i = 1:N-1
                inner_s5b = 0;
                for k =1:N+1
                    Enrich_s5b = 0;
                    for q=1:Q
                        Enrich_s5b = Enrich_s5b + XV2(i+1,q)*YV2(k,q) - Xu(i+1,q)*DSu(k,q);
                    end
                    inner_s5b = inner_s5b + Enrich_s5b*w(k)*SV2(k); 
                end
                Vs5b(i) = inner_s5b*w(i+1);
            end
             
            for i = 1:N-1
                outer_s5c = 0;
                for m=1:N+1
                    inner_s5c = 0;
                    for k =1:N+1
                        Enrich_s5c = 0;
                        for q=1:Q
                            Enrich_s5c = Enrich_s5c + DRV2(m,q)*YV2(k,q) - XV1(m,q)*DSV1(k,q);
                        end      
                        inner_s5c = inner_s5c + Enrich_s5c*w(k)*SV2(k);
                    end
                    outer_s5c = outer_s5c + inner_s5c*w(m)*D(m,i+1);
                end
                Vs5c(i) = outer_s5c;
            end
            Vs5a=Vs5a(:);
            Vs5b=Vs5b(:);
            Vs5c=Vs5c(:);
            Vs55= Vs5a+Vs5b+Vs5c;
            
            for i = 1:N+1
                outer_s6a = 0;
                for m=1:N+1
                    inner_s6a = 0;
                    for k =1:N-1
                        Enrich_s6a = 0;
                        for q=1:Q
                            Enrich_s6a = Enrich_s6a - DRV3(m,q)*YV3(k+1,q) - XV4(m,q)*DSV4(k+1,q) + Xp(m,q)*DSp(k+1,q);
                        end
                        inner_s6a = inner_s6a + Enrich_s6a*w(k+1)*SV3(k);       
                    end
                    outer_s6a = outer_s6a + inner_s6a*w(m)*D(m,i);
                end
                Vs6a(i) = - outer_s6a;
            end
             
            for i = 1:N+1
                inner_s6b = 0;
                for k =1:N-1
                    Enrich_s6b = 0;
                    for q=1:Q
                        Enrich_s6b = Enrich_s6b + XV3(i,q)*YV3(k+1,q) - DRv(i,q)*Yv(k+1,q);
                    end        
                    inner_s6b = inner_s6b + Enrich_s6b*w(k+1)*SV3(k);       
                end
                Vs6b(i) = inner_s6b*w(i);
            end
                             
            for i = 1:N+1
                outer_s6c = 0;
                for m=1:N+1
                        Enrich_s6c = 0;
                        for q=1:Q
                            Enrich_s6c = Enrich_s6c + DRV4(i,q)*YV4(m,q) - XV3(i,q)*DSV3(m,q);
                        end
                        inner_s6c = 0;
                        for l=1:N-1
                            inner_s6c = inner_s6c + D(m,l+1)*SV3(l);
                        end  
                    outer_s6c = outer_s6c + Enrich_s6c*inner_s6c*w(m);
                end
                Vs6c(i) = - outer_s6c*w(i);
            end
            Vs6a=Vs6a(:);
            Vs6b=Vs6b(:);
            Vs6c=Vs6c(:);
            Vs66 = Vs6a+Vs6b+Vs6c;
            
            for i = 1:N-1
                outer_s7a = 0;
                for m=1:N+1
                        Enrich_s7a = 0;
                        for q=1:Q
                            Enrich_s7a = Enrich_s7a - DRV3(i+1,q)*YV3(m,q) - XV4(i+1,q)*DSV4(m,q) + Xp(i+1,q)*DSp(m,q);
                        end
                        inner_s7a = 0;
                        for l = 1:N+1
                            inner_s7a = inner_s7a + D(m,l)*SV4(l);           
                        end  
                    outer_s7a = outer_s7a + Enrich_s7a*inner_s7a*w(m);
                end
                Vs7a(i) = - outer_s7a*w(i+1);
            end
             
            for i = 1:N-1
                inner_s7b = 0;
                for k =1:N+1
                    Enrich_s7b = 0;
                    for q=1:Q
                        Enrich_s7b = Enrich_s7b + XV4(i+1,q)*YV4(k,q) - Xv(i+1,q)*DSv(k,q);
                    end
                    inner_s7b = inner_s7b + Enrich_s7b*w(k)*SV4(k);           
                end
                Vs7b(i) = inner_s7b*w(i+1);
            end
             
            for i = 1:N-1
                outer_s7c = 0;
                for m=1:N+1
                    inner_s7c = 0;
                    for k =1:N+1
                        Enrich_s7c = 0;
                        for q=1:Q
                            Enrich_s7c = Enrich_s7c + DRV1(m,q)*YV1(k,q) + DRV4(m,q)*YV4(k,q);
                        end         
                        inner_s7c = inner_s7c + Enrich_s7c*w(k)*SV4(k);           
                    end
                    outer_s7c = outer_s7c + inner_s7c*w(m)*D(m,i+1);
                end
                Vs7c(i) = outer_s7c;
            end
            
            for i = 1:N-1
                outer_s7d = 0;
                for m=1:N+1
                        Enrich_s7d = 0;
                        for q=1:Q
                            Enrich_s7d = Enrich_s7d + XV1(i+1,q)*DSV1(m,q) + XV4(i+1,q)*DSV4(m,q);
                        end
                        inner_s7d = 0;
                        for l = 1:N+1          
                        inner_s7d = inner_s7d + D(m,l)*SV4(l);           
                        end 
                    outer_s7d = outer_s7d + Enrich_s7d*inner_s7d*w(m);
                end
                Vs7d(i) = outer_s7d*w(i+1);
            end

            for i = 1:N-1
                outer_s7e = 0;
                for m=1:N+1
                    inner_s7e = 0;
                    for k =1:N+1
                        Enrich_s7e = 0;
                        for q=1:Q
                            Enrich_s7e = Enrich_s7e + DRV4(m,q)*YV4(k,q) - XV3(m,q)*DSV3(k,q);
                        end         
                        inner_s7e = inner_s7e + Enrich_s7e*w(k)*SV4(k);           
                    end
                    outer_s7e = outer_s7e + inner_s7e*w(m)*D(m,i+1);
                end
                Vs7e(i) = outer_s7e;
            end
            
            Vs7a=Vs7a(:);
            Vs7b=Vs7b(:);
            Vs7c=Vs7c(:);
            Vs7d=Vs7d(:);
            Vs7e=Vs7e(:);
            Vs77= Vs7a+Vs7b+Vs7c+Vs7d+Vs7e;

            VS1 = Vs1 - Vs11;
            VS2 = Vs2 - Vs22;
            VS3 = Vs3 - Vs33;
            VS4 = Vs4 - Vs44;
            VS5 = Vs5 - Vs55;
            VS6 = Vs6 - Vs66;
            VS7 = Vs7 - Vs77;

            Vector1 = [ VS1 ; VS2 ; VS3 ; VS4 ; VS5 ; VS6 ; VS7 ];
            
            %First solution
            R1 = Matrix1\Vector1;
        
            Ru = R1(1:N-1);          
            Rv = R1(N:2*N-2);       
            Rp = R1(2*N-1:3*N-3);   
            RV1 = R1(3*N-2:4*N-2);  
            RV2 = R1(4*N-1:5*N-3);  
            RV3 = R1(5*N-2:6*N-2);  
            RV4 = R1(6*N-1:7*N-3);  


            %Entries for the second matrix system
            Ms11 = (Ru'*A_2*Ru)*M_2+(Ru'*A_2*Ru)*M_2+(Ru'*M_2*Ru)*A_2;
            Ms12 = (Ru'*L_2'*Rv)*L_2;
            Ms21= (Rv'*L_2*Ru)*L_2';
            Ms13 = zeros(N-1,N-1); 
            Ms31 = zeros(N-1,N-1);
            Ms14= -(Ru'*L_4'*RV1)*M_2;
            Ms41 = -(RV1'*L_4*Ru)*M_2;
            Ms15= -(Ru'*M_2*RV2)*L_4';
            Ms51 = -(RV2'*M_2*Ru)*L_4;
            Ms16= zeros(N-1,N-1);
            Ms61 = zeros(N-1,N-1);
            Ms17= zeros(N-1,N+1);
            Ms71 = zeros(N+1,N-1);
            
            Ms22= (Rv'*M_2*Rv)*A_2+(Rv'*A_2*Rv)*M_2+(Rv'*M_2*Rv)*A_2;
            Ms23= zeros(N-1,N-1);
            Ms32 = zeros(N-1,N-1);
            Ms24= zeros(N-1,N-1);
            Ms42 = zeros(N-1,N-1);
            Ms25= zeros(N-1,N+1);
            Ms52 = zeros(N+1,N-1);
            Ms26= -(Rv'*L_4'*RV3)*M_2;
            Ms62 = -(RV3'*L_4*Rv)*M_2;
            Ms27= -(Rv'*M_2*RV4)*L_4';
            Ms72 = -(RV4'*M_2*Rv)*L_4;
            
            J6_s33= H'*H*(sum(H*Rp))^2;
            Ms33 = (Rp'*Att_2*Rp)*Mtt_2+(Rp'*Mtt_2*Rp)*Att_2 + mu*J6_s33;     
            Ms34 = -(Rp'*At_4'*RV1)*Mt_2';
            Ms43 = -(RV1'*At_4*Rp)*Mt_2;
            Ms35 = -(Rp'*Lt_2'*RV2)*G_3 ;
            Ms53 = -(RV2'*Lt_2*Rp)*G_3' ;
            Ms36 = -(Rp'*G_3*RV3)*Lt_2';
            Ms63 = -(RV3'*G_3'*Rp)*Lt_2;
            Ms37 = -(Rp'*Mt_2'*RV4)*At_4';
            Ms73 = -(RV4'*Mt_2*Rp)*At_4;
           
            Ms44 = (RV1'*A*RV1)*M_2+(RV1'*M*RV1)*M_2+(RV1'*A*RV1)*M_2+(RV1'*M*RV1)*A_2+(RV1'*M*RV1)*A_2;
            Ms45 = zeros(N-1,N+1);
            Ms54 = zeros(N+1,N-1); 
            Ms46 = zeros(N-1,N-1);
            Ms64 = zeros(N-1,N-1);
            Ms47 = (RV1'*A_4*RV4)*M_3+(RV1'*M_4*RV4)*A_3;
            Ms74 = (RV4'*A_3*RV1)*M_4+(RV4'*M_3*RV1)*A_4;
            
            Ms55 = (RV2'*M_2*RV2)*A+(RV2'*M_2*RV2)*M+(RV2'*A_2*RV2)*M;
            Ms56 = zeros(N+1,N-1);
            Ms65 = zeros(N-1,N+1);
            Ms57 = zeros(N+1,N+1);
            Ms75 = zeros(N+1,N+1);
            
            Ms66 = (RV3'*A*RV3)*M_2+(RV3'*M*RV3)*M_2+(RV3'*M*RV3)*A_2;
            Ms67 = zeros(N-1,N+1);
            Ms76 = zeros(N+1,N-1); 
            
            Ms77 = (RV4'*M_2*RV4)*A+(RV4'*M_2*RV4)*M+(RV4'*A_2*RV4)*M+(RV4'*M_2*RV4)*A+(RV4'*A_2*RV4)*M;
            
            %second matrix
            Matrix2 = [Ms11, Ms12, Ms13, Ms14, Ms15, Ms16, Ms17; Ms21, Ms22, Ms23, Ms24, Ms25, Ms26, Ms27; Ms31, Ms32, Ms33, Ms34, Ms35, Ms36, Ms37; Ms41, Ms42, Ms43, Ms44, Ms45, Ms46, Ms47; Ms51, Ms52, Ms53, Ms54, Ms55, Ms56, Ms57; Ms61, Ms62, Ms63, Ms64, Ms65, Ms66, Ms67; Ms71, Ms72, Ms73, Ms74, Ms75, Ms76, Ms77] ;
  
            %Entries for second RHS vectors
            Vr1=0.*Ru;
            Vr2=0.*Rv;
                    
            %Vr3
            for k=1:N-1 
                summ1=0;
                summ2=0;
                    for m=1:NP1
                        outsum1=0;
                        outsum2=0;
                        insum1=0;
                        insum2=0;
                            for i=1:N-1
                                insum1=insum1+Dtilde(m,i+1)*Rp(i);
                                insum2=insum2+E(m,i)*Rp(i);
                            end
                            for n=1:NP1
                               outsum1= outsum1+w(n)*E(n,k)*F1(m,n);
                               outsum2= outsum2+w(n)*Dtilde(n,k+1)*F2(m,n);
                            end
                        summ1 = summ1 + w(m)*outsum1*insum1;
                        summ2 = summ2 + w(m)*outsum2*insum2;
                    end 
                Vr3(k)=summ1+summ2;
            end
            Vr3=Vr3(:);
         
            %Vr4
            for k=1:N-1 
                outsum=0;
                for m=1:NP1
                    insum=0;
                    for i=1:NP1
                        insum=insum+RV1(i)*D(m,i);
                    end
                    outsum=outsum+insum*w(m)*w(k+1)*F1(m,k+1);
                end
                Vr4(k)=-outsum;
            end
            Vr4=Vr4(:);

            %Vr5
            for k=1:NP1 
                outsum=0;
                for n=1:NP1
                    insum=0;
                    for m=1:N-1
                        insum=insum+w(m+1)*F1(m+1,n)*RV2(m);
                    end
                    outsum=outsum+insum*w(n)*D(n,k);
                end
                Vr5(k)=-outsum;
            end
            Vr5=Vr5(:);

            %Vr6
            for k=1:N-1
                outsum=0;
                for m=1:NP1
                    insum=0;
                    for i=1:NP1
                        insum=insum+RV3(i)*D(m,i);
                    end
                    outsum=outsum+insum*w(k+1)*w(m)*F2(m,k+1);
                end
                Vr6(k)=-outsum;
            end
            Vr6=Vr6(:);
            
            %Vr7
            for k=1:NP1 
                outsum=0;
                for n=1:NP1 
                    insum=0;
                    for m=1:N-1
                        insum=insum+w(m+1)*F2(m+1,n)*RV4(m);
                    end
                    outsum=outsum+insum*w(n)*D(n,k);
                end
                Vr7(k)=-outsum;
            end
            Vr7=Vr7(:);
        
        % NEW ENRICHMENT TERMS
            %FOR VR1
             for k = 1:N-1
                outer_r1a = 0;
                for m=1:N+1
                        Enrich_r1a = 0;
                        for q=1:Q
                            Enrich_r1a = Enrich_r1a + DRu(m,q)*Yu(k+1,q) + Xv(m,q)*DSv(k+1,q);
                        end
                        inner_r1a = 0;
                        for l = 1:N-1
                            inner_r1a = inner_r1a + D(m,l+1)*Ru(l);
                        end 
                    outer_r1a = outer_r1a + Enrich_r1a*inner_r1a*w(m);
                end
                Vr1a(k) = outer_r1a*w(k+1);
             end
             
              for k = 1:N-1
                outer_r1b = 0;
                for m=1:N+1
                        Enrich_r1b = 0;
                        for q=1:Q
                            Enrich_r1b = Enrich_r1b + XV1(m,q)*YV1(k+1,q) - DRu(m,q)*Yu(k+1,q);
                        end
                        inner_r1b = 0;
                        for l = 1:N-1
                            inner_r1b = inner_r1b + D(m,l+1)*Ru(l); 
                        end  
                    outer_r1b = outer_r1b + Enrich_r1b*inner_r1b*w(m); %%%%%this has Enrich_s1b
                end
                Vr1b(k) = -outer_r1b*w(k+1);
              end

              for k = 1:N-1
                outer_r1c = 0;
                for m=1:N+1
                    inner_r1c = 0;
                    for i =1:N-1
                        Enrich_r1c = 0;
                        for q=1:Q
                            Enrich_r1c = Enrich_r1c + XV2(i+1,q)*YV2(m,q) - Xu(i+1,q)*DSu(m,q);
                        end                 
                        inner_r1c = inner_r1c + Enrich_r1c*w(i+1)*Ru(i);
                    end
                    outer_r1c = outer_r1c + inner_r1c*w(m)*D(m,k+1);
                end
                Vr1c(k) = -outer_r1c;
              end
               
            Vr1a=Vr1a(:);
            Vr1b=Vr1b(:);
            Vr1c=Vr1c(:);
            Vr11= Vr1a+Vr1b+Vr1c;

            %FOR VR2
             for k = 1:N-1
                outer_r2a = 0;
                for m=1:N+1
                    inner_r2a = 0;
                    for i =1:N-1
                        Enrich_r2a = 0;
                        for q=1:Q
                            Enrich_r2a = Enrich_r2a + DRu(i+1,q)*Yu(m,q) + Xv(i+1,q)*DSv(m,q);
                        end
                        inner_r2a = inner_r2a + Enrich_r2a*w(i+1)*Rv(i);  
                    end
                    outer_r2a = outer_r2a + inner_r2a*w(m)*D(m,k+1);
                end
                Vr2a(k) = outer_r2a;  
             end

             for k = 1:N-1
                outer_r2b = 0;
                for m=1:N+1
                        Enrich_r2b = 0;
                        for q=1:Q
                            Enrich_r2b = Enrich_r2b + XV3(m,q)*YV3(k+1,q) - DRv(m,q)*Yv(k+1,q);
                        end
                        inner_r2b = 0;
                        for l = 1:N-1
                            inner_r2b = inner_r2b + D(m,l+1)*Rv(l);
                        end
                    outer_r2b = outer_r2b + Enrich_r2b*inner_r2b*w(m);
                end
                Vr2b(k) = -outer_r2b*w(k+1);
             end
            
             for k = 1:N-1
                outer_r2c = 0;
                for m=1:N+1
                    inner_r2c = 0;
                    for i =1:N-1
                        Enrich_r2c = 0;
                        for q=1:Q
                            Enrich_r2c = Enrich_r2c + XV4(i+1,q)*YV4(m,q) - Xv(i+1,q)*DSv(m,q);
                        end         
                        inner_r2c = inner_r2c + Enrich_r2c*w(i+1)*Rv(i);
                    end
                    outer_r2c = outer_r2c + inner_r2c*w(m)*D(m,k+1);
                end
                Vr2c(k) = -outer_r2c;
            end
            Vr2a=Vr2a(:);
            Vr2b=Vr2b(:);
            Vr2c=Vr2c(:);
            Vr22= Vr2a+Vr2b+Vr2c;

            %FOR VR3
             for k = 1:N-1
                outer_r3a = 0;
                for n=1:N+1
                    middle_r3a = 0;
                    for m = 1:N+1
                        Enrich_r3a = 0;
                        for q=1:Q
                            Enrich_r3a = Enrich_r3a +(- DRV1(m,q)*YV1(n,q) - XV2(m,q)*DSV2(n,q) + DRp(m,q)*Yp(n,q));
                        end
                        inner_r3a = 0;
                        for l = 1:N-1
                            inner_r3a = inner_r3a + Dtilde(m,l+1)*Rp(l);     
                        end
                        middle_r3a = middle_r3a + Enrich_r3a*inner_r3a*w(m);
                    end
                    outer_r3a = outer_r3a + middle_r3a*E(n,k)*w(n);
                end
                Vr3a(k) = outer_r3a;
             end

              for k = 1:N-1
                outer_r3b = 0;
                for n=1:N+1
                    middle_r3b = 0;
                    for m = 1:N+1
                        Enrich_r3b = 0;
                        for q=1:Q
                            Enrich_r3b = Enrich_r3b +(- DRV3(m,q)*YV3(n,q) - XV4(m,q)*DSV4(n,q) + Xp(m,q)*DSp(n,q));
                        end
                        inner_r3b = 0;
                        for l = 1:N-1           
                            inner_r3b = inner_r3b + E(m,l)*Rp(l); 
                        end
                        middle_r3b = middle_r3b + Enrich_r3b*inner_r3b*w(m);
                    end
                    outer_r3b = outer_r3b + middle_r3b*Dtilde(n,k+1)*w(n);
                end
                Vr3b(k) = outer_r3b;
              end

             for i = 1:N-1
                middle_r3c = 0;
                for m=1:N+1
                    inner_r3c = 0;
                    for n=1:N+1
                        Enrich_r3c = 0;
                        for q=1:Q
                            Enrich_r3c = Enrich_r3c +(Xp(m,q)*Yp(n,q));
                        end
                        inner_r3c = inner_r3c + Enrich_r3c; 
                    end
                    middle_r3c = middle_r3c + inner_r3c;
                end
                outer_r3c = 0;
                    for j=1:N-1          
                       outer_r3c = outer_r3c + Rp(j)*H(j);                 
                    end
                Vr3c(i) = middle_r3c*outer_r3c*H(i)*mu;
             end
              
            Vr3a=Vr3a(:);
            Vr3b=Vr3b(:);
            Vr3c=Vr3c(:);
            Vr33= Vr3a+Vr3b +Vr3c;


            %FOR Vr4
            for k = 1:N-1
                outer_r4a = 0;
                for m=1:N+1
                    inner_r4a = 0;
                    for i =1:N+1
                        Enrich_r4a = 0;
                        for q=1:Q
                            Enrich_r4a = Enrich_r4a - DRV1(m,q)*YV1(k+1,q) - XV2(m,q)*DSV2(k+1,q) + DRp(m,q)*Yp(k+1,q);
                        end
                        inner_r4a = 0;
                        for l = 1:N+1
                            inner_r4a = inner_r4a + D(m,l)*RV1(l);   
                        end                    
                    end
                    outer_r4a = outer_r4a + Enrich_r4a*inner_r4a*w(m);
                end
                Vr4a(k) = - outer_r4a*w(k+1);
            end
            
            for k = 1:N-1
                outer_r4b = 0;
                    inner_r4b = 0;
                    for i =1:N+1
                        Enrich_r4b = 0;
                        for q=1:Q
                            Enrich_r4b = Enrich_r4b + XV1(i,q)*YV1(k+1,q) - DRu(i,q)*Yu(k+1,q);
                        end       
                        inner_r4b = inner_r4b + Enrich_r4b*w(i)*RV1(i);           
                    end
                Vr4b(k) = inner_r4b*w(k+1);
            end

            for k = 1:N-1
                outer_r4c = 0;
                for m=1:N+1
                        Enrich_r4c = 0;
                        for q=1:Q
                            Enrich_r4c = Enrich_r4c + DRV1(m,q)*YV1(k+1,q) + DRV4(m,q)*YV4(k+1,q);
                        end
                        inner_r4c = 0;
                        for l = 1:N+1         
                            inner_r4c = inner_r4c + D(m,l)*RV1(l); 
                        end   
                    outer_r4c = outer_r4c + Enrich_r4c*inner_r4c*w(m);
                end
                Vr4c(k) = outer_r4c*w(k+1);
            end

            for k = 1:N-1
                outer_r4d = 0;
                for m=1:N+1
                    inner_r4d = 0;
                    for i =1:N+1
                        Enrich_r4d = 0;
                        for q=1:Q
                            Enrich_r4d = Enrich_r4d + XV1(i,q)*DSV1(m,q) + XV4(i,q)*DSV4(m,q);
                        end                  
                        inner_r4d = inner_r4d + Enrich_r4d*w(i)*RV1(i);           
                    end
                    outer_r4d = outer_r4d + inner_r4d*w(m)*D(m,k+1);
                end
                Vr4d(k) = outer_r4d;
            end

            for k = 1:N-1
                outer_r4e = 0;
                for m=1:N+1
                    inner_r4e = 0;
                    for i =1:N+1
                        Enrich_r4e = 0;
                        for q=1:Q
                            Enrich_r4e = Enrich_r4e + DRV2(i,q)*YV2(m,q) - XV1(i,q)*DSV1(m,q);
                        end                 
                        inner_r4e = inner_r4e + Enrich_r4e*w(i)*RV1(i);           
                    end
                    outer_r4e = outer_r4e + inner_r4e*w(m)*D(m,k+1);
                end
                Vr4e(k) = -outer_r4e;
            end
            Vr4a=Vr4a(:);
            Vr4b=Vr4b(:);
            Vr4c=Vr4c(:);
            Vr4d=Vr4d(:);
            Vr4e=Vr4e(:);
            Vr44= Vr4a+Vr4b+Vr4c+Vr4d+Vr4e;

            %FOR Vr5
             for k = 1:N+1
                outer_r5a = 0;
                for m=1:N+1
                    inner_r5a = 0;
                    for i =1:N-1
                        Enrich_r5a = 0;
                        for q=1:Q 
                            Enrich_r5a = Enrich_r5a - DRV1(i+1,q)*YV1(m,q) - XV2(i+1,q)*DSV2(m,q) + DRp(i+1,q)*Yp(m,q);
                        end
                        inner_r5a = inner_r5a + Enrich_r5a*w(i+1)*RV2(i);       
                    end
                    outer_r5a = outer_r5a + inner_r5a*w(m)*D(m,k);
                end
                Vr5a(k) = - outer_r5a;
             end
            
             for k = 1:N+1
                outer_r5b = 0;
                    inner_r5b = 0;
                    for i =1:N-1
                        Enrich_r5b = 0;
                        for q=1:Q 
                            Enrich_r5b = Enrich_r5b + XV2(i+1,q)*YV2(k,q) - Xu(i+1,q)*DSu(k,q);
                        end        
                        inner_r5b = inner_r5b + Enrich_r5b*w(i+1)*RV2(i);
                    end
                Vr5b(k) = inner_r5b*w(k);
             end

             for k = 1:N+1
                outer_r5c = 0;
                for m=1:N+1
                        Enrich_r5c = 0;
                        for q=1:Q 
                            Enrich_r5c = Enrich_r5c + DRV2(m,q)*YV2(k,q) - XV1(m,q)*DSV1(k,q);
                        end
                        inner_r5c = 0;
                        for l = 1:N-1
                            inner_r5c = inner_r5c + D(m,l+1)*RV2(l);
                        end
                    outer_r5c = outer_r5c + Enrich_r5c*inner_r5c*w(m);
                end
                Vr5c(k) = outer_r5c*w(k);
            end
            Vr5a=Vr5a(:);
            Vr5b=Vr5b(:);
            Vr5c=Vr5c(:);
            Vr55 = Vr5a+Vr5b+Vr5c;


            %FOR Vr6
             for k = 1:N-1
                outer_r6a = 0;
                for m=1:N+1
                        Enrich_r6a = 0;
                        for q=1:Q
                            Enrich_r6a = Enrich_r6a - DRV3(m,q)*YV3(k+1,q) - XV4(m,q)*DSV4(k+1,q) + Xp(m,q)*DSp(k+1,q);
                        end
                        inner_r6a = 0;
                        for l = 1:N+1
                            inner_r6a = inner_r6a + D(m,l)*RV3(l); 
                        end   
                    outer_r6a = outer_r6a + Enrich_r6a*inner_r6a*w(m);
                end
                Vr6a(k) = -outer_r6a*w(k+1); %%
             end
            
             for k = 1:N-1
                    inner_r6b = 0;
                    for i =1:N+1
                        Enrich_r6b = 0;
                        for q=1:Q
                            Enrich_r6b = Enrich_r6b + XV3(i,q)*YV3(k+1,q) - DRv(i,q)*Yv(k+1,q);
                        end        
                        inner_r6b = inner_r6b + Enrich_r6b*w(i)*RV3(i);  
                    end
                Vr6b(k) = inner_r6b*w(k+1); 
             end
            
             for k = 1:N-1
                outer_r6c = 0;
                for m=1:N+1
                    inner_r6c = 0;
                    for i =1:N+1
                        Enrich_r6c = 0;
                        for q=1:Q
                            Enrich_r6c = Enrich_r6c + DRV4(i,q)*YV4(m,q) - XV3(i,q)*DSV3(m,q); %%%%%%%%%%
                        end                
                        inner_r6c = inner_r6c + Enrich_r6c*w(i)*RV3(i);
                    end
                    outer_r6c = outer_r6c + inner_r6c*w(m)*D(m,k+1);
                end
                Vr6c(k) = -outer_r6c;
            end
            Vr6a=Vr6a(:);
            Vr6b=Vr6b(:);
            Vr6c=Vr6c(:);
            Vr66 = Vr6a+Vr6b+Vr6c;

            %FOR Vr7
             for k = 1:N+1
                outer_r7a = 0;
                for m=1:N+1
                    inner_r7a = 0;
                    for i =1:N-1
                        Enrich_r7a = 0;
                        for q=1:Q
                            Enrich_r7a = Enrich_r7a - DRV3(i+1,q)*YV3(m,q) - XV4(i+1,q)*DSV4(m,q) + Xp(i+1,q)*DSp(m,q);
                        end
                        inner_r7a = inner_r7a + Enrich_r7a*w(i+1)*RV4(i);               
                    end
                    outer_r7a = outer_r7a + inner_r7a*w(m)*D(m,k);
                end
                Vr7a(k) = - outer_r7a;
             end

             for k = 1:N+1
                inner_r7b = 0;
                for i =1:N-1
                    Enrich_r7b = 0;
                    for q=1:Q
                        Enrich_r7b = Enrich_r7b + XV4(i+1,q)*YV4(k,q) - Xv(i+1,q)*DSv(k,q);
                    end
                    inner_r7b = inner_r7b + Enrich_r7b*w(i+1)*RV4(i);              
                end
                Vr7b(k) = inner_r7b*w(k);
             end

             for k = 1:N+1
                outer_r7c = 0;
                for m=1:N+1
                    for i =1:N-1
                        Enrich_r7c = 0;
                        for q=1:Q
                            Enrich_r7c = Enrich_r7c + DRV1(m,q)*YV1(k,q) + DRV4(m,q)*YV4(k,q);
                        end
                        inner_r7c = 0;
                        for l = 1:N-1        
                            inner_r7c = inner_r7c + D(m,l+1)*RV4(l);  
                        end            
                    end
                    outer_r7c = outer_r7c + Enrich_r7c*inner_r7c*w(m);
                end
                Vr7c(k) = outer_r7c*w(k);
             end

             for k = 1:N+1
                 outer_r7d = 0;
                 for m=1:N+1
                     inner_r7d = 0;
                     for i =1:N-1
                         Enrich_r7d = 0;
                         for q=1:Q
                             Enrich_r7d = Enrich_r7d + XV1(i+1,q)*DSV1(m,q) + XV4(i+1,q)*DSV4(m,q);
                         end        
                         inner_r7d = inner_r7d + Enrich_r7d*w(i+1)*RV4(i);              
                     end
                     outer_r7d = outer_r7d + inner_r7d*w(m)*D(m,k);
                 end
                 Vr7d(k) = outer_r7d;
             end

             for k = 1:N+1
                outer_r7e = 0;
                for m=1:N+1
                        Enrich_r7e = 0;
                        for q=1:Q
                            Enrich_r7e = Enrich_r7e + DRV4(m,q)*YV4(k,q) - XV3(m,q)*DSV3(k,q);
                        end
                        inner_r7e = 0;
                        for l = 1:N-1                
                            inner_r7e = inner_r7e + D(m,l+1)*RV4(l); 
                        end   
                    outer_r7e = outer_r7e + Enrich_r7e*inner_r7e*w(m);
                end
                Vr7e(k) = outer_r7e*w(k);
            end
            Vr7a=Vr7a(:);
            Vr7b=Vr7b(:);
            Vr7c=Vr7c(:);
            Vr7d=Vr7d(:);
            Vr7e=Vr7e(:);
            Vr77= Vr7a+Vr7b+Vr7c+Vr7d+Vr7e;

            VR1 = Vr1 - Vr11;
            VR2 = Vr2 - Vr22;
            VR3 = Vr3 - Vr33;
            VR4 = Vr4 - Vr44;
            VR5 = Vr5 - Vr55;
            VR6 = Vr6 - Vr66;
            VR7 = Vr7 - Vr77;
            
            
            Vector2 = [ VR1 ; VR2 ; VR3 ; VR4 ; VR5 ; VR6 ; VR7 ];
             
            S1 = Matrix2\Vector2;
            
            Su = S1(1:N-1);          
            Sv = S1(N:2*N-2);       
            Sp = S1(2*N-1:3*N-3);    
            SV1 = S1(3*N-2:4*N-4);   
            SV2 = S1(4*N-3:5*N-3);   
            SV3 = S1(5*N-2:6*N-4);  
            SV4 = S1(6*N-3:7*N-3);   
        
            counter=counter+1;
            if counter == maxiterations
                break
            end

        end
        
        endTime = cputime;

        cumcounter = [cumcounter + counter];
        % Append the current cumcounter value to the array
        cumcounter_values = [cumcounter_values, cumcounter];
        each_count_value = [each_count_value, counter];
      
        % Calculate execution time for this iteration
        iterationTime = endTime - startTime;

        cumTime(number) = cumTime(number) + iterationTime;
       
        %setting XQ YQ from the calculated values
        xu= Ru ;
        xv = Rv ;
        xp = Rp ;
        xV1 = RV1 ;
        xV2 = RV2 ;
        xV3 = RV3 ;
        xV4 = RV4 ;
             
        yu = Su ;
        yv = Sv ;
        yp = Sp ;
        yV1 = SV1 ;
        yV2 = SV2 ;
        yV3 = SV3 ;
        yV4 = SV4 ;
               
        %adding the known boundary values
        Xu = [Xu,[0;xu;0]];
        Xv = [Xv, [0;xv;0]];
        XV1 = [XV1, xV1];
        XV2 = [XV2, [0;xV2;0]];
        XV3 = [XV3, xV3];
        XV4 = [XV4, [0;xV4;0]];
        Yu = [Yu, [0;yu;0]];
        Yv = [Yv,[0;yv;0]];
        YV1 = [YV1, [0;yV1;0]];
        YV2 = [YV2, yV2];
        YV3 = [YV3, [0;yV3;0]];
        YV4 = [YV4, yV4];
       
        %Pressure boundary terms
        Xp1 = 0;
        for i = 1:N-1
            Xp1 = Xp1 + E(1,i)*xp(i);
        end
        XpNP1 = 0;
        for i = 1:N-1
            XpNP1 = XpNP1 + E(N+1,i)*xp(i);
        end
        Yp1 = 0;
        for i = 1:N-1
            Yp1 = Yp1 + E(1,i)*yp(i);
        end
        YpNP1 = 0;
        for i = 1:N-1
            YpNP1 = YpNP1 + E(N+1,i)*yp(i);
        end
        
        Xp= [Xp, [Xp1;xp;XpNP1]];
        Yp= [Yp, [Yp1;yp;YpNP1]];
       
        Q=Q+1;  
   
        executionTimes(Q,number) = cumTime(number);

        %calculating the approximations
        approxu = zeros(N+1) ;
        approxv = zeros(N+1) ;
        approxp = zeros(N+1) ;
        approxV1 = zeros(N+1) ;
        approxV2 = zeros(N+1) ;
        approxV3 = zeros(N+1) ;
        approxV4 = zeros(N+1) ;
        aapproxV4 = zeros(N+1); 
        
        for k = 1 :N+1
            for l = 1 :N+1
                for q = 1 : Q
                    approxu(l,k) = approxu(l,k) + Yu(k,q)*Xu(l,q) ;
                    approxv(l,k) = approxv(l,k) + Yv(k,q)*Xv(l,q) ;
                    approxp(l,k) = approxp(l,k) + Yp(k,q)*Xp(l,q) ;
                    approxV1(l,k) = approxV1(l,k) + YV1(k,q)*XV1(l,q) ;
                    approxV2(l,k) = approxV2(l,k) + YV2(k,q)*XV2(l,q) ;
                    approxV3(l,k) = approxV3(l,k) + YV3(k,q)*XV3(l,q) ;
                    approxV4(l,k) = approxV4(l,k) + YV4(k,q)*XV4(l,q) ;
                    
                end
            end
        end
          
        Uerror = zeros(N+1) ;
        Verror = zeros(N+1) ;
        Perror = zeros(N+1) ;
        V1error = zeros(N+1) ;
        V2error = zeros(N+1) ;
        V3error = zeros(N+1) ;
        V4error = zeros(N+1) ;
        
        for k = 1 :N+1
            for l = 1 :N+1
                Uerror(k,l) = approxu(k,l) - Umesh(k,l);
                Verror(k,l) = approxv(k,l) - Vmesh(k,l);
                Perror(k,l) = approxp(k,l) - Pmesh(k,l);
                V1error(k,l) = approxV1(k,l) - V1mesh(k,l);
                V2error(k,l) = approxV2(k,l) - V2mesh(k,l);
                V3error(k,l) = approxV3(k,l) - V3mesh(k,l);
                V4error(k,l) = approxV4(k,l) - V4mesh(k,l);

            end
        end 
      
        L2Uerror=0;
        L2Verror=0;
        L2Perror=0;
        L2V1error=0;
        L2V2error=0;
        L2V3error=0;
        L2V4error=0;
        for num=1:NP1
            for m=1:NP1
                L2Uerror=L2Uerror+((w(m)*w(num)*(approxu(m,num)-Umesh(m,num))^2)) ;
                L2Verror=L2Verror+((w(m)*w(num)*(approxv(m,num)-Vmesh(m,num))^2)) ;
                L2Perror=L2Perror+((w(m)*w(num)*(approxp(m,num)-Pmesh(m,num))^2)) ;
                L2V1error=L2V1error+((w(m)*w(num)*(approxV1(m,num)-V1mesh(m,num))^2)) ;
                L2V2error=L2V2error+((w(m)*w(num)*(approxV2(m,num)-V2mesh(m,num))^2)) ;
                L2V3error=L2V3error+((w(m)*w(num)*(approxV3(m,num)-V3mesh(m,num))^2)) ;
                L2V4error=L2V4error+((w(m)*w(num)*(approxV4(m,num)-V4mesh(m,num))^2)) ;
                   
            end
        end


        L2Uerror= L2Uerror^(1/2);
        L2Verror= L2Verror^(1/2);
        L2Perror= L2Perror^(1/2);
        L2V1error= L2V1error^(1/2);
        L2V2error= L2V2error^(1/2);
        L2V3error= L2V3error^(1/2);
        L2V4error= L2V4error^(1/2);
        
        L2U(Q,number)=L2Uerror;
        L2V(Q,number)=L2Verror;
        L2P(Q,number)=L2Perror;
        L2V1(Q,number)=L2V1error;
        L2V2(Q,number)=L2V2error;
        L2V3(Q,number)=L2V3error;
        L2V4(Q,number)=L2V4error;
 

    end

    allb_cumcounter(number,:) = cumcounter_values;
    each_count_b(number,:) = each_count_value;



end
 

 


N_values=Nvalues'; 
 
L2Uvec = sqrt(L2U.^2 + L2V.^2);
  
% Display the execution times
disp('Execution times for each iteration:');
disp(executionTimes);
  
fprintf('L2 errors for Example 9, using the LS SM PGD, after one enrichment:')
fprintf('\n');
fprintf('\n'); 
fprintf('N =       ');
fprintf('%8d ', Nvalues);
fprintf('\n-----------------------------------------------------\n');
fprintf('L2 error p   ');
fprintf('%6.2e ', L2P(end,:));
fprintf('\nL2 error u   ');
fprintf('%6.2e ', L2Uvec(end,:));
fprintf('\n');

each_count_bT = each_count_b';
fprintf('\n');
fprintf('Number of ADFPA iterations required at the first enrichment stage for Example 9 using the LS SM PGD, with epsilon = :')
disp(epsilon);
fprintf('\n');
fprintf('\n'); 
fprintf('N =          ');
fprintf('%8d ', Nvalues);
fprintf('\n-----------------------------------------------------\n');
fprintf('ADFPA iter   ');
fprintf('%8d ', each_count_bT); 
fprintf('\n');
 