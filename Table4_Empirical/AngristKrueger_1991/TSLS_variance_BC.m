function [b, Rho, S11, Var_AS,Var_BC] = TSLS_variance_BC(y, Y2, X, X1, k, T, g,L, G,K)          
%  	Bias corrected variance estimator
    
    [b, R] = kClass(y, Y2, X, X1, 1, T); 
    
    Yb = X*(X\Y2);
    Q = [Yb'*Yb Yb'*X1; X1'*Yb X1'*X1];
    Z1 = [Y2 X1];
    e =     y - Z1*b;
    I = eye( G-1+k, G-1+k );
    S11 =    e'*e/(T- g - k);
    Var_AS = diag(S11*(I/Q)); 
    
    V2 = Y2 - Yb;
    v1=y-X*(X\y);
    Tau = V2'*(v1-V2*b(1:g))/(T-K);
    C               = zeros( G-1+k, G-1+k );
    C(1:G-1,1:G-1)  = V2'*V2/(T-K);
    C1               = zeros( G-1+k, G-1+k );
    C1(1:G-1,1:G-1)  =(Tau*Tau')/S11;
   
    Var_bias = S11*((L+1)*(Q\C1/Q) + trace(Q\C)*(I/Q));
    
    Var_BC = Var_AS - diag(Var_bias);
    if sum(Var_BC<=0)>=1  | sum(Var_BC>=Var_AS)>=1
       Var_BC = Var_AS;
       Var_BC = TSLS_variance_boot(y, Y2, X, X1, T, g, k, K, 999)' ;
    end
    Rho = corr(e, V2);
end