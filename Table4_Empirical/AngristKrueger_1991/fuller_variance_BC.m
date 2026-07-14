
function [b, Rho,S11,Var_AS,Var_BC] = fuller_variance_BC(y, Y2, X, X1,k, T, g,L, G,K)          
%  	Bias corrected variance estimator
    
    Yd = [y Y2];
    YdtYd = Yd'*Yd;
    Wsdd = YdtYd - ( Yd' * X1 * (X1\Yd) );
    Wdd = YdtYd - ( Yd' * X * (X\Yd) );
    lambda = min(eig(Wdd\Wsdd));
    lambda = lambda - 1/(T-K);
    
    [b, R] = kClass(y, Y2, X, X1, lambda, T);  
                                
    V2 = Y2 - X*(X\Y2);
    UL = Y2'*Y2-lambda*(V2'*V2);
    UR = Y2'*X1;
    LL = UR';
    LR = X1'*X1;
    Q1= ([UL UR; LL LR]); 
        
    Z1 = [Y2 X1];
    e =     y - Z1*b;
    S11 =   e'*e/(T- g - k);
    I = eye(size(Q1));
    Var_AS = diag(S11*(I/Q1));
    
    Yb= X*(X\Y2);
    Q = [Yb'*Yb Yb'*X1; X1'*Yb X1'*X1];
    v1=y-X*(X\y);
    Tau = V2'*(v1-V2*b(1:g))/(T);
    C               = zeros( G-1+k, G-1+k );
    C(1:G-1,1:G-1)  = V2'*V2/(T);
    C1               = zeros( G-1+k, G-1+k );
    C1(1:G-1,1:G-1)  =(Tau*Tau')/S11;
    C2              = C - C1;
   
    Var_bias =  S11*((Q\C/Q) +5*(Q\C1/Q) + 2*trace(Q\C)*(I/Q)- trace(Q\C1)*(I/Q) - trace(Q\C2)*(I/Q) - L*(Q\C2/Q));% (L*(T-K+L-2)+2*L)/(T-K-2)*(Q\C2/Q));
    %Var_bias =  S11*(6*(Q\C1/Q) + trace(Q\C)*(I/Q) - (L-1)*(Q\C2/Q));% (L*(T-K+L-2)+2*L)/(T-K-2)*(Q\C2/Q));
    Var_BC = Var_AS - diag(Var_bias);
    if sum(Var_BC<=0)>=1
       %Var_BC =  Var_AS;
       Var_BC = fuller_variance_boot(y, Y2, X, X1, T, g, k, K, 999)' ;
    end
    Rho = corr(e, V2);
end