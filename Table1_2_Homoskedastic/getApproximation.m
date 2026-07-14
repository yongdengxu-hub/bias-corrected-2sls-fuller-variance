%function a = getApproximation(y0, X, X1, Beta, Sigma, Omega, Pi, Gamma, G, g, K, k, J, j, indY, indLY, exdLY, eq, L, T, eList, rvInxL0, rvInxL1)
function a = getApproximation(y0, X, X1, Beta, Sigma, Omega, Pi, Gamma, G, g, K, k, J, j, indY, indLY, exdLY, eq, L, T, eList)
    eListLen = length(eList);
    a = zeros(eListLen, g+k+j);
    
    for eptr = 1:eListLen
        switch eList(eptr)
            case 10001;  a(eptr,:) = OLS_bias_approximation_SS(X, X1, Beta, Sigma, Omega, Pi, G, g, k, eq, T);
            case 10101;  a(eptr,:) = OLS_bias_approximation_LT(X, X1, Beta, Sigma, Omega, Pi, G, g, k, eq, T);
            %case 10201;  a(eptr,:) = OLS_con_bias_approximation(X, X1, Beta, Sigma, Omega, Pi, G, g, k, eq, T);
                
            %case 10102;  a(eptr,:) = OLS_bias_approximation(X, X1, Beta, Sigma, Omega, Pi, G, g, k, eq, T);
            case 10102;  a(eptr,:) = NaN;
            case 10202;  a(eptr,:) = OLS_bias_approximation_LT(X, X1, Beta, Sigma, Omega, Pi, G, g, k, eq, T);
                
            %case 20001;  a(eptr,:) = Mikhail_2SLS_bias_approximation(X, X1, Beta, Sigma, Omega, Pi, G, k, eq, L);
            case 20101;  a(eptr,:) = Kclass_bias_approximation(X, X1, Beta, Sigma, Pi, k, 1, eq, L, T);    
            
            
            %case 20201;  a(eptr,:) = Dynamic_2SLS_LT_bias_approximation(y0, X, X1, Beta, Sigma, Omega, Pi, Gamma, G, g, K, k, J, j, indY, indLY, exdLY, T, eq);
            %case 20201;  NaN;
            %case 20301;  a(eptr,:) = Dynamic_2SLS_SS_bias_approximation(y0, X, X1, Beta, Sigma, Omega, Pi, Gamma, G, g, K, k, J, j, indY, indLY, T, eq);
            
            case 20401;  a(eptr,:) = Kclass_bias_approximation(X, X1, Beta, Sigma, Pi, k, 1-1/T^3, eq, L, T);
            case 20501;  a(eptr,:) = Kclass_bias_approximation(X, X1, Beta, Sigma, Pi, k, 1-1/T^2, eq, L, T);    
            case 20601;  a(eptr,:) = Kclass_bias_approximation(X, X1, Beta, Sigma, Pi, k, 1-1/T^(3/2), eq, L, T);    
            case 20701;  a(eptr,:) = Kclass_bias_approximation(X, X1, Beta, Sigma, Pi, k, 1-1/T, eq, L, T);
            case 20801;  a(eptr,:) = Kclass_comb_bias_approximation(X, X1, Beta, Sigma, Pi, k, 1-1/T^3, 1-1/T, L, L-1, eq, L, T);
            case 20901;  a(eptr,:) = Kclass_comb_bias_approximation(X, X1, Beta, Sigma, Pi, k, 1-1/T, 1-(L+1)/T, 2, 1, eq, L, T);
                
            %case 20002;  a(eptr,:) = Mikhail_2SLS_bias_approximation(X, X1, Beta, Sigma, Omega, Pi, G, k, eq, L);
            case 20102;  a(eptr,:) = NaN;
            %case 20102;  a(eptr,:) = Kclass_bias_approximation(X, X1, Beta, Sigma, Pi, k, 1, eq, L, T);    
            
            case 20402;  a(eptr,:) = NaN;
            case 20502;  a(eptr,:) = NaN;
            case 20602;  a(eptr,:) = NaN;
            case 20702;  a(eptr,:) = NaN;
            case 20802;  a(eptr,:) = NaN;
            case 20902;  a(eptr,:) = NaN;
                
            case 21001;  a(eptr,:) = NaN;
            case 21101;  a(eptr,:) = NaN;
                
            %case 20602;  a(eptr,:) = Kclass_comb_bias_approximation(X, X1, Beta, Sigma, Pi, k, 1-1/T^3, 1-1/T, L, L-1, eq, L, T);
            %case 20702;  a(eptr,:) = Kclass_comb_bias_approximation(X, X1, Beta, Sigma, Pi, k, 1-1/T, 1-(L+1)/T, 2, 1, eq, L, T);
            
            case 30101;  a(eptr,:) = NaN;
            case 30201;  a(eptr,:) = NaN;
            case 30301;  a(eptr,:) = NaN;
            case 30401;  a(eptr,:) = NaN;
            case 30501;  a(eptr,:) = NaN;
                
            case 30102;  a(eptr,:) = NaN;
            case 30202;  a(eptr,:) = NaN;
            case 30302;  a(eptr,:) = NaN;
            case 30402;  a(eptr,:) = NaN;
            case 30502;  a(eptr,:) = NaN;
           
            case 40101;  a(eptr,:) = NaN;
            case 40201;  a(eptr,:) = NaN;
            case 40301;  a(eptr,:) = NaN;
            case 40401;  a(eptr,:) = NaN;
            case 40501;  a(eptr,:) = NaN;
                
            case 40102;  a(eptr,:) = NaN;
            case 40202;  a(eptr,:) = NaN;
            case 40302;  a(eptr,:) = NaN;
            case 40402;  a(eptr,:) = NaN;
            case 40502;  a(eptr,:) = NaN;
            
            case 50101;  a(eptr,:) = NaN;
            case 50201;  a(eptr,:) = NaN;
            case 50301;  a(eptr,:) = NaN;
            case 50401;  a(eptr,:) = NaN;
            case 50501;  a(eptr,:) = NaN;
                
            case 50102;  a(eptr,:) = NaN;
            case 50202;  a(eptr,:) = NaN;
            case 50302;  a(eptr,:) = NaN;
            case 50402;  a(eptr,:) = NaN;
            case 50502;  a(eptr,:) = NaN;
            
            case 60101;  a(eptr,:) = NaN;
            case 60201;  a(eptr,:) = NaN;
            case 60301;  a(eptr,:) = NaN;
            case 60401;  a(eptr,:) = NaN;
            case 60501;  a(eptr,:) = NaN;
                
            case 60102;  a(eptr,:) = NaN;
            case 60202;  a(eptr,:) = NaN;
            case 60302;  a(eptr,:) = NaN;
            case 60402;  a(eptr,:) = NaN;
            case 60502;  a(eptr,:) = NaN;
            
            case 70101;  a(eptr,:) = NaN;
            case 70201;  a(eptr,:) = NaN;
            case 70301;  a(eptr,:) = NaN;
                
            case 70102;  a(eptr,:) = NaN;
            case 70202;  a(eptr,:) = NaN;
            case 70302;  a(eptr,:) = NaN;
            
            case 80101;  a(eptr,:) = Kclass_bias_approximation(X, X1, Beta, Sigma, Pi, k, 1, eq, L, T);    
            
            case 80102;  a(eptr,:) = NaN;
            case 80201;  a(eptr,:) = NaN;
            case 80202;  a(eptr,:) = NaN;              
                
          
            %{
            case 30101;  a(eptr,:) = RV1_bias_approximation(X, X1, Beta, Sigma, Pi, g, K, k, 1, eq, T, rvInxL0);
            case 30201;  a(eptr,:) = RV1_bias_approximation(X, X1, Beta, Sigma, Pi, g, K, k, 1-1/T^3, eq, T, rvInxL0);
            case 30301;  a(eptr,:) = RV1_bias_approximation(X, X1, Beta, Sigma, Pi, g, K, k, 1-1/T^2, eq, T, rvInxL0);
            case 30401;  a(eptr,:) = RV1_bias_approximation(X, X1, Beta, Sigma, Pi, g, K, k, 1-1/T^(3/2), eq, T, rvInxL0);    
            case 30501;  a(eptr,:) = RV1_bias_approximation(X, X1, Beta, Sigma, Pi, g, K, k, 1-1/T, eq, T, rvInxL0);
                
            case 30102;  a(eptr,:) = RV1_bias_approximation(X, X1, Beta, Sigma, Pi, g, K, k, 1, eq, T, rvInxL0);
            case 30202;  a(eptr,:) = RV1_bias_approximation(X, X1, Beta, Sigma, Pi, g, K, k, 1-1/T^3, eq, T, rvInxL0);
            case 30302;  a(eptr,:) = RV1_bias_approximation(X, X1, Beta, Sigma, Pi, g, K, k, 1-1/T^2, eq, T, rvInxL0);
            case 30402;  a(eptr,:) = RV1_bias_approximation(X, X1, Beta, Sigma, Pi, g, K, k, 1-1/T^(3/2), eq, T, rvInxL0);    
            case 30502;  a(eptr,:) = RV1_bias_approximation(X, X1, Beta, Sigma, Pi, g, K, k, 1-1/T, eq, T, rvInxL0);    
           
            case 40101;  a(eptr,:) = RV1_bias_approximation(X, X1, Beta, Sigma, Pi, g, K, k, 1, eq, T, rvInxL1);
            case 40201;  a(eptr,:) = RV1_bias_approximation(X, X1, Beta, Sigma, Pi, g, K, k, 1-1/T^3, eq, T, rvInxL1);    
            case 40301;  a(eptr,:) = RV1_bias_approximation(X, X1, Beta, Sigma, Pi, g, K, k, 1-1/T^2, eq, T, rvInxL1);
            case 40401;  a(eptr,:) = RV1_bias_approximation(X, X1, Beta, Sigma, Pi, g, K, k, 1-1/T^(3/2), eq, T, rvInxL1);    
            case 40501;  a(eptr,:) = RV1_bias_approximation(X, X1, Beta, Sigma, Pi, g, K, k, 1-1/T, eq, T, rvInxL1);
                
            case 40102;  a(eptr,:) = RV1_bias_approximation(X, X1, Beta, Sigma, Pi, g, K, k, 1, eq, T, rvInxL1);
            case 40202;  a(eptr,:) = RV1_bias_approximation(X, X1, Beta, Sigma, Pi, g, K, k, 1-1/T^3, eq, T, rvInxL1);    
            case 40302;  a(eptr,:) = RV1_bias_approximation(X, X1, Beta, Sigma, Pi, g, K, k, 1-1/T^2, eq, T, rvInxL1);
            case 40402;  a(eptr,:) = RV1_bias_approximation(X, X1, Beta, Sigma, Pi, g, K, k, 1-1/T^(3/2), eq, T, rvInxL1);    
            case 40502;  a(eptr,:) = RV1_bias_approximation(X, X1, Beta, Sigma, Pi, g, K, k, 1-1/T, eq, T, rvInxL1);
            
            case 50101;  a(eptr,:) = RV2_bias_approximation(X, X1, Beta, Sigma, Pi, g, k, 1, eq, T, rvInxL0);          
            case 50201;  a(eptr,:) = RV2_bias_approximation(X, X1, Beta, Sigma, Pi, g, k, 1-1/T^3, eq, T, rvInxL0);    
            case 50301;  a(eptr,:) = RV2_bias_approximation(X, X1, Beta, Sigma, Pi, g, k, 1-1/T^2, eq, T, rvInxL0);        
            case 50401;  a(eptr,:) = RV2_bias_approximation(X, X1, Beta, Sigma, Pi, g, k, 1-1/T^(3/2), eq, T, rvInxL0);    
            case 50501;  a(eptr,:) = RV2_bias_approximation(X, X1, Beta, Sigma, Pi, g, k, 1-1/T, eq, T, rvInxL0);    
                
            case 50102;  a(eptr,:) = RV2_bias_approximation(X, X1, Beta, Sigma, Pi, g, k, 1, eq, T, rvInxL0);          
            case 50202;  a(eptr,:) = RV2_bias_approximation(X, X1, Beta, Sigma, Pi, g, k, 1-1/T^3, eq, T, rvInxL0);    
            case 50302;  a(eptr,:) = RV2_bias_approximation(X, X1, Beta, Sigma, Pi, g, k, 1-1/T^2, eq, T, rvInxL0);        
            case 50402;  a(eptr,:) = RV2_bias_approximation(X, X1, Beta, Sigma, Pi, g, k, 1-1/T^(3/2), eq, T, rvInxL0);    
            case 50502;  a(eptr,:) = RV2_bias_approximation(X, X1, Beta, Sigma, Pi, g, k, 1-1/T, eq, T, rvInxL0);    
            
            case 60101;  a(eptr,:) = RV2_bias_approximation(X, X1, Beta, Sigma, Pi, g, k, 1, eq, T, rvInxL1);          
            case 60201;  a(eptr,:) = RV2_bias_approximation(X, X1, Beta, Sigma, Pi, g, k, 1-1/T^3, eq, T, rvInxL1);    
            case 60301;  a(eptr,:) = RV2_bias_approximation(X, X1, Beta, Sigma, Pi, g, k, 1-1/T^2, eq, T, rvInxL1);        
            case 60401;  a(eptr,:) = RV2_bias_approximation(X, X1, Beta, Sigma, Pi, g, k, 1-1/T^(3/2), eq, T, rvInxL1);    
            case 60501;  a(eptr,:) = RV2_bias_approximation(X, X1, Beta, Sigma, Pi, g, k, 1-1/T, eq, T, rvInxL1);    
                
            case 60102;  a(eptr,:) = RV2_bias_approximation(X, X1, Beta, Sigma, Pi, g, k, 1, eq, T, rvInxL1);          
            case 60202;  a(eptr,:) = RV2_bias_approximation(X, X1, Beta, Sigma, Pi, g, k, 1-1/T^3, eq, T, rvInxL1);    
            case 60302;  a(eptr,:) = RV2_bias_approximation(X, X1, Beta, Sigma, Pi, g, k, 1-1/T^2, eq, T, rvInxL1);        
            case 60402;  a(eptr,:) = RV2_bias_approximation(X, X1, Beta, Sigma, Pi, g, k, 1-1/T^(3/2), eq, T, rvInxL1);    
            case 60502;  a(eptr,:) = RV2_bias_approximation(X, X1, Beta, Sigma, Pi, g, k, 1-1/T, eq, T, rvInxL1);     
            
            case 70101;  a(eptr,:) = NaN;
            case 70201;  a(eptr,:) = NaN;
            case 70301;  a(eptr,:) = NaN;
                
            case 70102;  a(eptr,:) = NaN;
            case 70202;  a(eptr,:) = NaN;
            case 70302;  a(eptr,:) = NaN;
            %}
        end
    end
end

function a = OLS_bias_approximation_LT(X, X1, Beta, Sigma, Omega, Pi, G, g, k, eq, T)
    C               = zeros( G-1+k, G-1+k );
    Omega(eq,:)     = [];
    Omega(:,eq)     = [];
    C(1:G-1,1:G-1)  = Omega;
            
    Pi2         = Pi;
    Pi2(:,eq)   = [];
    
    E               = zeros( G-1+k );
    E(1:G-1,1:G-1)  = T*Omega;
    
    Yb  = X*Pi2;
    Q   = [Yb'*Yb Yb'*X1; X1'*Yb X1'*X1] + E;
            
    Beta2       = inv(Beta);
    Beta2(:,eq) = [];
    
    S1 = Sigma(:,eq);
    qb = [Beta2'*S1; zeros(k,1)];
            
    
    T1 = (T-(g+k+1))*(Q\qb);
    T2 = T*trace(Q\C)*(Q\qb);
    T3 = (T*(g+k+2))*(Q\C)*(Q\qb);
    T4 = -1*T^2*trace(Q\C)*(Q\C)*(Q\qb);
    T5 = -1*T^2*(Q\C)*(Q\C)*(Q\qb);
    
    a = [T1  T2  T3  T4  T5];
    a = a';
    a = sum(a);
end

function a = OLS_bias_approximation_SS(X, X1, Beta, Sigma, Omega, Pi, G, g, k, eq, T)
    %%C               = zeros( G-1+k, G-1+k );
    %%Omega(eq,:)     = [];
    %%Omega(:,eq)     = [];
    %%C(1:G-1,1:G-1)  = Omega;
    
    Pi2         = Pi;
    Pi2(:,eq)   = [];
    
    %%E               = zeros( G-1+k );
    %%E(1:G-1,1:G-1)  = T*Omega;
    
    %Yb  = X*Pi2;
    %Q   = [Yb'*Yb Yb'*X1; X1'*Yb X1'*X1] + E;
    Zb = [X*Pi2 X1];
    Q = inv(Zb'*Zb);
            
    Beta2       = inv(Beta);
    Beta2(:,eq) = [];
    
    S1 = Sigma(:,eq);
    qb = [Beta2'*S1; zeros(k,1)];

    
    a = (T-(g+k+1))*Q*qb;
    
    %%T1 = (T-(g+k+1))*(Q\q);
    %%T2 = T*trace(Q\C)*(Q\q);
    %%T3 = (T*(g+k+2))*(Q\C)*(Q\q);
    %%T4 = -1*T^2*trace(Q\C)*(Q\C)*(Q\q);
    %%T5 = -1*T^2*(Q\C)*(Q\C)*(Q\q);
    
    %%a = [T1  T2  T3  T4  T5];
    %%a = a';
    %%a = sum(a);
end

function a = OLS_con_bias_approximation(X, X1, Beta, Sigma, Omega, Pi, G, g, k, eq, T)
    x1 = T/(2*(g+k+1));
    c1 = 2*x1 - 1;
    c2 = 2*x1 - 2;

    a1 = OLS_bias_approximation(X(1:T/2,:), X1(1:T/2,:), Beta, Sigma, Omega, Pi, G, g, k, eq, T/2);
    a2 = OLS_bias_approximation(X(T/2+1:T,:), X1(T/2+1:T,:), Beta, Sigma, Omega, Pi, G, g, k, eq, T/2);
    a3 = OLS_bias_approximation(X, X1, Beta, Sigma, Omega, Pi, G, g, k, eq, T);
     a = c1*((a1+a2)/2) - c2*a3;
end

function a = Kclass_bias_approximation(X, X1, Beta, Sigma, Pi, k, kc, eq, L, T)
    Pi2       = Pi;
    Pi2(:,eq) = [];
    
    Yb = X*Pi2;
    Q = [Yb'*Yb Yb'*X1; X1'*Yb X1'*X1];
             
    Beta2 = inv(Beta);
    Beta2(:,eq) = [];
    
    S1 = Sigma(:,eq);
    q = [Beta2'*S1; zeros(k,1)];
    
    theta = (kc-1)*T; 
    %a = (L-theta-1);
    ccc = (Q\q);
    a = (L-theta-1)*(Q\q);
end

%{
function a = Kclass_comb1_bias_approximation(X, X1, Beta, Sigma, Pi, k, eq, L, T)
    a1 = Kclass_bias_approximation(X, X1, Beta, Sigma, Pi, k, 1-1/T^3, eq, L, T);
    %a1 = Kclass_bias_approximation(X, X1, Beta, Sigma, Pi, k, 1, eq, L, T);
    a2 = Kclass_bias_approximation(X, X1, Beta, Sigma, Pi, k, 1-1/T, eq, L, T);
    a = L*a1 - (L-1)*a2;
end

function a = Kclass_comb2_bias_approximation(X, X1, Beta, Sigma, Pi, k, eq, L, T)
    a1 = Kclass_bias_approximation(X, X1, Beta, Sigma, Pi, k, 1-1/T, eq, L, T);
    a2 = Kclass_bias_approximation(X, X1, Beta, Sigma, Pi, k, 1-(L+1)/T, eq, L, T);
    a = 2*a1 - a2;
end
%}

function a = Kclass_comb_bias_approximation(X, X1, Beta, Sigma, Pi, k, kc1, kc2, c1, c2, eq, L, T)
    a1 = Kclass_bias_approximation(X, X1, Beta, Sigma, Pi, k, kc1, eq, L, T);
    a2 = Kclass_bias_approximation(X, X1, Beta, Sigma, Pi, k, kc2, eq, L, T);
    a = c1*a1 - c2*a2;
end

function a = Mikhail_2SLS_bias_approximation(X, X1, Beta, Sigma, Omega, Pi, G, k, eq, L)
    C               = zeros( G-1+k, G-1+k );
    Omega(eq,:)     = [];
    Omega(:,eq)     = [];
    C(1:G-1,1:G-1)  = Omega;
            
    Pi2 = Pi;
    Pi2(eq,:) = [];
    Yb = X*Pi2;
    Q = [Yb'*Yb Yb'*X1; X1'*Yb X1'*X1];
    
    Beta2 = inv(Beta);
    Beta2(:,eq) = [];
    S1 = Sigma(:,eq);
    q = [Beta2'*S1; zeros(k,1)];
    
    I = eye( G-1+k, G-1+k );
    a = (L-1)*(I + trace(Q\C)*I - (L-2)*(Q\C))*(Q\q);
end

function a = Dynamic_2SLS_LT_bias_approximation(y0, X, X1, Beta, Sigma, Omega, Pi, Gamma, G, g, K, k, J, j, indY,indLY, exdLY, T, eq)
    Gamma2 = Gamma(:,indY);

    Qw = zeros(G,G);
    for t = 1:T-1
        Qw = Qw + (T-t)*(Gamma^(t-1))'*Omega*(Gamma^(t-1));
        %ELWLW = ELWLW + (Gamma^(t-1))'*Omega*(Gamma^(t-1));
    end
    
    %%%%%%%
    %{
    Qwww = zeros(G,G);
    for t = 0:T-1
        Qwww = Qwww + (T-t)*(Gamma^t)'*Omega*(Gamma^t);
        %ELWLW = ELWLW + (Gamma^(t-1))'*Omega*(Gamma^(t-1));
    end 
    %}
    %%%%%%%%
    
    Yb = zeros(T,G);
    Yb(1,:) = y0*Gamma + X(1,:)*Pi;
    for t = 2:T
            Yb(t,:) = Yb(t-1,:)*Gamma + X(t,:)*Pi;
    end
    Y2b = Yb(:,indY);
    
    if isempty(indLY)
        Rb = [Y2b X1];
        Zb = X;
        I2s = zeros(K,G);
        
        %Qz = Zb'*Zb + [Qw(indLY,indLY) zeros(indLY,K); zeros(K, indLY) zeros(K, K)];
        Qz = inv(Zb'*Zb);
        Qzs = I2s'*Qz*I2s;
    else
        LYb = [y0; Yb(1:T-1,:)];
        LY1b = LYb(:,indLY);
        Rb = [Y2b LY1b X1];
        %Zb = [LYb X];
        %I2s = [eye(G); zeros(K,G)];
        %Qz = Zb'*Zb + [Qw zeros(G,K); zeros(K, G) zeros(K, K)];
        %Qz = inv(Qz);
        %Qzs = I2s'*Qz*I2s;
        
        Zb = [LYb(:,sort([indLY, exdLY])) X];
        zvector = zeros(1,G);
        zvector(:,sort([indLY exdLY])) = 1;
        Qwz = Qw;
        dropOff = setdiff(1:1:G,sort([indLY, exdLY]));
        Qwz(:,dropOff) = [];
        Qwz(dropOff,:) = [];
        I2s = [diag(zvector); zeros(K,G)];
        I2s(dropOff,:) = [];
        Qz = Zb'*Zb + [Qwz zeros(J,K); zeros(K, J) zeros(K, K)];
        %Qz = Zb'*Zb + [Qw(:,[indLY exdLY]) zeros(G,K); zeros(K, G) zeros(K, K)];
        Qz = inv(Qz);
        Qzs = I2s'*Qz*I2s;
        
        %Zb = [LYb X];
        %zvector = zeros(1,G);
        %zvector(:,[indLY exdLY]) = 1;
        %I2s = [diag(zvector); zeros(K,G)];
        %Qz = Zb'*Zb + [Qw zeros(G,K); zeros(K, G) zeros(K, K)];
        %Qz = inv(Qz);
        %Qzs = I2s'*Qz*I2s;
    end
        
    %Zb = [LYb(:,indLY) X];
    
    
    Ed1d1 = [Gamma2'*Qw*Gamma2 Gamma2'*Qw(:,indLY) zeros(g,k); ...
             Qw(indLY,:)*Gamma2 Qw(indLY,indLY) zeros(j,k); ...
             zeros(k,g) zeros(k,j) zeros(k,k)];
    %Ed1d1 = [Gamma2'*Qwww*Gamma2 0 zeros(g,k); ...
    %         0 0 zeros(j,k); ...
    %         zeros(k,g) zeros(k,j) zeros(k,k)];

    Qs = Rb'*Rb + Ed1d1;
    Qs = inv(Qs);
    
%{
    Zb = [LYb X];
    Qz = Zb'*Zb + [Qw zeros(G,K); zeros(K, G+K)];
    Qz = inv(Qz);
    I2s = [eye(G); zeros(K,G)];
    Qzs = I2s'*Qz*I2s;
             %}
    
    %Zb = [LYb(:,indLY) X];
    %Qz = Zb'*Zb + [Qw(indLY,indLY) zeros(indLY,K); zeros(K, indLY) zeros(K, K)];
    %Qz = inv(Qz);
    %Qzs = I2s'*Qz*I2s;

    D = diag(ones(1,T-1),-1);
    
    
    
    %B = [eye(g+j) zeros(g+j,k)];
    I = eye(G);
    B = [I(:,1:g) zeros(G,j+k)];
    
    %psi = B'*vphi;
    phi = (inv(Beta')*Sigma(:,eq));
    
    
    %Beta2 = inv(Beta);
    Beta2 = inv(Beta);
    Beta2 = Beta2(:,indY); 
    
    %psi = (Beta2'*Sigma(:,eq)) / Sigma(eq,eq);
    psi = (Beta2'*Sigma(:,eq));
    
    %psi = [psi zeros(g,j+k)]';
    psi = [psi; zeros(j+k,1)];
    
    I1 = [eye(j); zeros(G-j,j)];
    A = [Gamma2 I1 zeros(G,k)];
    %A = [Gamma2 [1;1] zeros(G,k)]; 
   
    I = eye(g+j+k);
    I2 = eye(G);
    %%F = zeros(g+j+k,1);
    
    T1 = zeros(g+j+k,1);
    T2 = T1; T3 = T1; T4 = T1; T5 = T1; T6 = T1; T7 = T1; T8 = T1; T9 = T1; T10 = T1; T11 = T1; T12 = T1;
    T13 = T1;
    T14 = T1;
    
    
    %T21 = T2;
    %T22 = T2;
    T9 = zeros(T,1);
    T11 = zeros(G,1);
    
    T1 = (Rb'*Zb*Qz*Zb'*Rb*Qs + (trace(Zb*Qz*Zb'*Rb*Qs*Rb')*I))*psi;
    T5 = ( (trace(Qw*I2s'*Qz*Zb'*Rb*Qs*A')*I) + Rb'*Zb*Qz*I2s*Qw*A*Qs + A'*Qw*I2s'*Qz*Zb'*Rb*Qs )*psi;
    
    T8 = (A'*Qw*Qzs*Qw*A*Qs + (trace(Qw*Qzs*Qw*A*Qs*A')*I))*psi;
    
    
    for t = 1:T-1
        T2 = T2 + (Rb'*(D^t)*Rb*Qs + (trace(Rb'*(D^t)*Rb*Qs)*I))*A'*(Gamma^(t-1))'* phi;
        %T21 = T21 + (Rb'*D^t*Rb*Qs*A'*(Gamma^(t-1))')*phi; 
        %22 = T22 + (A'*(Gamma^(t-1))'*(trace(Rb'*D^t*Rb*Qs)*I2))*phi;
        for r = 1:T-1
            T3  = T3  + (Rb'*(D^t)*(D^r)'*Rb*Qs + (trace((D^t)*(D^r)'*Rb*Qs*Rb')*I))*(trace(Omega*(Gamma^(r-1))*Qzs*(Gamma^(t-1))')*I)*psi;
            T4  = T4  + (trace((D^t)*(D^r)'*Zb*Qz*I2s*(Gamma^(t-1))'*Omega*(Gamma^(r-1))*A*Qs*Rb')*I)*psi;
            T6  = T6  + (trace(Zb*Qz*Zb'*(D^t)*(D^r)')*I) * ((trace(Omega*(Gamma^(t-1))*A*Qs*A'*(Gamma^(r-1))')*I) + A'*(Gamma^(t-1))'*Omega*(Gamma^(r-1))*A*Qs)*psi;
            T7  = T7  + (A'*(Gamma^(t-1))'*Omega*(Gamma^(r-1))*I2s'*Qz*Zb'*(D^t)*(D^r)'*Rb*Qs*psi);
            %T9  = T9  + ( Rb'*(D^t)*(D^r)*Rb*Qs*B'*Omega*(Gamma^(t-1))*Qzs*(Gamma^(r-1))'*phi);
            T9  = T9  + ((D^t)*(D^r)*Rb*Qs*B'*Omega*(Gamma^(t-1))*Qzs*(Gamma^(r-1))'*phi);
          
            T10 = T10 + ( A'*(Gamma^(t-1))'* (Omega*B*Qs*Rb'*((D^r)'*(D^t)'+(D^r)'*(D^t))*Zb*Qz*I2s*(Gamma^(r-1))' + ...
                (trace((D^t)'*(D^r)*Zb*Qz*I2s*(Gamma^(r-1))'*Omega*B*Qs*Rb')*eye(G)))*phi);
            %T11 = T11 + ( A'*(Gamma^(t-1))'* (trace(Omega*B*Qs*A'*(Gamma^(r-1))')*trace((D^t)'*Zb*Qz*Zb'*(D^r)')*eye(G)) * phi);
            T11 = T11 + ( (Gamma^(t-1))'* (trace(Omega*B*Qs*A'*(Gamma^(r-1))')*trace((D^t)'*Zb*Qz*Zb'*(D^r)')*eye(G)) * phi);
            T12 = T12 + ( B'*Omega*(Gamma^(r-1))* (Qzs*(Gamma^(t-1))'*(trace((D^t)*Rb*Qs*Rb'*(D^r))*eye(G)) + ...
                A*Qs*A'*(Gamma^(t-1))'*(trace(Zb*Qz*Zb'*(D^t)*(D^r))*eye(G)) + A*Qs*Rb'*(D^t)'*(D^r)'*Zb*Qz*I2s*(Gamma^(t-1))')*phi);
            for s = 1:T-1
                if t == r+s 
                    T13 = T13 + ( A'*(Gamma^(t-1))'*Omega*Gamma^(t-r-1)*A*Qs*A'*(Gamma^(r-1))' * trace((D^t)'*(D^r)*(D^(t-r)))*eye(G))*phi;
                elseif r == t+s
                    T14 = T14 + ( A'*(Gamma^(t-1))'* trace(Gamma^(r-t-1)*A*Qs*A'*(Gamma^(r-1))'*Omega) * trace((D^t)'*(D^r)*(D^(r-t))')*eye(G))*phi;
                else
                end
            end
        end
    end
    
    T1 =  - Qs*T1;
    T2 =  + Qs*(trace(Zb*Qz*Zb')*I)*psi - Qs*T2;
    T3 =  - Qs*T3;
    T4  = - Qs*T4;
    T5  = - Qs*T5;
    T6  = - Qs*T6;
    T7  = - Qs*T7;
    T8  = - Qs*T8;
    %T9  = - Qs*T9;
    T9  = - Qs*Rb'*T9;
    T10 = - Qs*T10;
    %T11 = - Qs*T11;
    T11 = - Qs*A'*T11;
    T12 = - Qs*T12;
    T13 = - Qs*T13;
    T14 = - Qs*T14;
    
    a = T1 + T2 + T3 + T4 + T5 + T6 + T7 + T8 + T9 + T10 + T11 + T12 + T13 + T14;
    %aa = [T1 T2 T3 T4 T5 T6 T7 T8 T9 T10 T11 T12 T13 T14]
    
    
    %psi = Beta'*vphi;
    %
    %jjjj = Rb'*Zb*Qz*Zb'*Rb*Qs;
    %kk = Zb*Qz*Zb'*Rb*Qs*Rb';
    %I = eye(size(A,2));
    %T1 = Qs * (Rb'*Zb*Qz*Zb'*Rb*Qs + (trace(Zb*Qz*Zb'*Rb*Qs*Rb')*I)) * vphi;
end

function a = Dynamic_2SLS_SS_bias_approximation(y0, X, X1, Beta, Sigma, Omega, Pi, Gamma, G, g, K, k, J, j, indY, indLY, T, eq)
    
    %Gamma2 = (Gamma(indLY,:))';
    Gamma2 = Gamma(:,indY);
    %Gamma2 = [0;0];

    %Sigma = iSigma/iSigma(eq,eq);
    %Omega = iOmega/iOmega(eq,eq);
    %sigma = Sigma(eq,eq);
    
    Yb = zeros(T,G);
    Yb(1,:) = y0*Gamma + X(1,:)*Pi;
    for t = 2:T
            Yb(t,:) = Yb(t-1,:)*Gamma + X(t,:)*Pi;
    end
    %Y1b = Yb(:,indLY);
    Y2b = Yb(:,indY);
    LYb = [y0; Yb(1:T-1,:)];
    LY1b = LYb(:,indLY);
    Rb = [Y2b LY1b X1];
    
    
    %phi = (inv(Beta')*Sigma(:,eq)) / Sigma(eq,eq);
    phi = (inv(Beta')*Sigma(:,eq));
    
    
    %Beta2 = inv(Beta);
    Beta2 = inv(Beta);
    Beta2 = Beta2(:,indY); 
    
    %psi = (Beta2'*Sigma(:,eq)) / Sigma(eq,eq);
    psi = (Beta2'*Sigma(:,eq));
    
    psi = [psi zeros(g,j+k)]';
    
    
    I1 = [eye(j); zeros(G-j,j)];
    A = [Gamma2 I1 zeros(G,k)];

    D = diag(ones(1,T-1),-1);
    
    T1 = zeros(g+j+k,G);
    I = eye(g+j+k);
    for t = 1:T
        T1 = T1 + (Rb'*(D^t)*Rb*inv(Rb'*Rb) + trace(Rb'*(D^t)*Rb*inv(Rb'*Rb))*I)*A'*(Gamma^(t-1))';
    end
    
    %a = sigma*(G+K-(g+j+k)-1)*inv(Rb'*Rb)*psi - sigma*inv(Rb'*Rb)*T1*phi;
    %a = sigma*(G+K-2*g-k-2)*inv(Rb'*Rb)*psi - sigma*inv(Rb'*Rb)*T1*phi;
    a = (G+K-(g+j+k)-1)*inv(Rb'*Rb)*psi - inv(Rb'*Rb)*T1*phi;
    %a = (G+K-2*g-k-2)*inv(Rb'*Rb)*psi - inv(Rb'*Rb)*T1*phi;
    
    %{
    
    
    
    Q = inv(Rb'*Rb);
    
    D = diag(ones(1,T-1),-1);
    
    I1 = eye(g+1);
    I1 = I1(:,indLY);
    A = [Gamma2 I1 zeros(G,k)];
    
    sigma = Sigma(eq,eq);
    %{
    EWu1 = zeros(g,1);
    Sigma1 = Sigma(:,eq);
    sigma = Sigma(eq,eq);
    for t = 0:T-1
        EWu1 = EWu1 + (T-t)*(Gamma^t)'*inv(Beta')*Sigma1;
    end
    %}
    
    Beta2 = inv(Beta);
    Beta2(:,eq) = []; 
    EV2u1 = Beta2'*Sigma(:,eq);
    vphi2 = (1/sigma)*EV2u1;
    vphi = [EV2u1(1:eq-1,:); 0; EV2u1(eq:end,:)'];
    psi = [vphi2 zeros(1,j+k)]';
       
    
    T1 = zeros(g+j+k,1);
    
    I = eye(g+j+k);
    for t = 1:T-1
        T1 = T1 + (Rb'*(D^t)*Rb*Q + trace(Rb'*(D^t)*Rb*Q)*I)*(Gamma^(t-1)*A)'*vphi;
    end
    
    a = (J+K-(g+j+k)-1)*Q*psi - Q*T1;
    %}
end

function a = RV1_bias_approximation(X, X1, Beta, Sigma, Pi, g, K, k, kc, eq, T, rvInx)
    coefLen = g + k;
    X1 = [X1 X(:,rvInx)];
    
    k = k + length(rvInx);
    L = K - k - g;
    a = Kclass_bias_approximation(X, X1, Beta, Sigma, Pi, k, kc, eq, L, T);
    a = a(1:coefLen,:);
end

function a = RV2_bias_approximation(X, X1, Beta, Sigma, Pi, g, k, kc, eq, T, rvInx)
    X(:,rvInx) = [];
    Pi(:, rvInx) = [];
    K = size(X,2);
    L = K - k - g;
    a = Kclass_bias_approximation(X, X1, Beta, Sigma, Pi, k, kc, eq, L, T);
end

%{
function [a Qq] = OLS_bias_approsimation(X, X1, eq, Beta, Gamma, Sigma,Omega,Pi,G,K,G,k,T,eq)
    Omega = inv(Beta) * Sigma *(inv(Beta))';
    Pi = -1 * inv(Beta) * Gamma;
    G = size(Beta,2);
    k = length(find(Gamma(eq,:) ~= 0)) + length(find(sum(Gamma)==0));
    g = length(find(Beta(eq,:) ~= 0)) - 1;
        
    C = zeros( G-1+k, G-1+k );
    C(1:G-1,1:G-1) = Omega(2:end,2:end);
            
    Pi2 = Pi';
    Pi2 = Pi2(:,2:end);
    E = zeros( G-1+k );
    E(1:G-1,1:G-1) = T*Omega(2:end,2:end);
    Q = [Pi2'*(X'*X)*Pi2 Pi2'*X'*X1; X1'*X*Pi2 X1'*X1] + E;
    Q = inv(Q);
            
    Beta2 = inv(Beta');
    Beta2 = Beta2(:,2:end);
    S1 = Sigma(:,1);
    q = [Beta2'*S1; zeros(k,1)];
            
    %g = length(indY)-1;
    %k = length(indX);
            
    T1 = (T-(g+k+1))*Q*q;
    T2 = T*trace(Q*C)*Q*q;
    T3 = (T*(g+k+2))*Q*C*Q*q;
    T4 = -1*T^2*trace(Q*C)*Q*C*Q*q;
    T5 = -1*T^2*Q*C*Q*C*Q*q;
    a = [T1  T2  T3  T4  T5];
    Qq = Q;
end
%}

