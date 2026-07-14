function [b, I] = kClass(y, Y2, X, X1, k, T)
       %V = (eye(T) - X * inv(X'*X) * X') * Y2;
       %X = [X(:,1) 10*X(:,2:end)];
       %X1 = [X1(:,1) 10*X1(:,2:end)];
       %X = 10*X;
       V = Y2 - X*(X\Y2);
    
        UL = Y2'*Y2-k*(V'*V);
        UR = Y2'*X1;
       %LL = X1'*Y2;
        LL = UR';
        LR = X1'*X1;
    
       %b = inv([UL UR; LL LR])*[(Y2- k*V)'*y; X1'*y];
        b = ([UL UR; LL LR])\([(Y2- k*V)'*y; X1'*y]);
        %b(5) = b(5)*100;
       %R1 = 1 - (V'*V) / (Y2'*Y2 - size(Y2,2)* (mean(Y2))^2);
        %M = eye(T) - ones(T,T)/T;
        %R2 = diag(eye(size(Y2,2)) - (V'*V)./(Y2'*M*Y2));
        %I = diag(eye(size(Y2,2)) - ((V'*V)./(T-size(X,2)))./((Y2'*M*Y2)./(T-1)));
        %I(2) = 1-((T-1)/(T-size(X,2)))*(1-R2);
        %I(1) = Y2'*M*Y2;
        %I(2) = V'*V;
        I = NaN; 
        
    %concentration parameter
    %I(2) = (1/o22) * p22'* X2' * (X2-X1*(X1\X2)) * p22;
end