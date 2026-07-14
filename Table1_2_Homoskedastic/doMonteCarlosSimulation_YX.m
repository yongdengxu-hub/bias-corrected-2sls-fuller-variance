function doMonteCarlosSimulation_YX (config)
    replications    = config.replications;  % number of replications
    eq              = config.eq;            % esitmated equation
    bst             = config.bst;           % number of bootstrapping
    constT          = config.constT;        % with/without constant term
    
    mListLen = length(config.model);
    tListLen = length(config.tList);
    dListLen = length(config.dList);
    eListLen = length(config.eList);
    
    eList    = config.eList;
    C       = config.model(1).C;
    k1=nnz(C(:,1))+1;
    O = zeros(mListLen,eListLen,tListLen,dListLen,k1,13);
    Tstats=NaN(mListLen,k1,eListLen,k1,replications);                 
    FirstF=NaN(mListLen,replications);

    V_All = zeros(mListLen,k1,eListLen,replications);
    S11_All = NaN(mListLen,k1,eListLen,replications);
    Var_AS_All = NaN(mListLen,k1,eListLen,replications);
    Var_BC_All = NaN(mListLen,k1,eListLen,replications);


    %** LAYER 1: Model*****************************************************
    for mptr = 1:mListLen
        Beta    = config.model(mptr).Beta;
        C       = config.model(mptr).C;
        Lambda  = config.model(mptr).Lambda;
        Sigma   = config.model(mptr).Sigma;
        Gamma   = -1 * Lambda/Beta;
        Pi      = -1 * C/Beta;
        Omega   = Beta'\Sigma/Beta;
        EigVal  = eig(Gamma);
        rho     = chol(Omega,'lower');
        
        % Get included (endogenous/lagged endogenous/exgenous) index
        indX        = find(C(:,eq) ~= 0);
        indY        = find(Beta(:,eq) ~= 0);
        indY(eq,:)  = [];
        indLY       = find(Lambda(:,eq) ~= 0);
        
        % Ignore the index of (Pi) with row elements all equal to zero
        inxX0   = find(sum(abs(C),2) ~= 0);
        exdX    = find(C(:,eq) == 0);
        exdX    = intersect(exdX, inxX0);
        
        % Ignore the index of (Beta) with row elements all equal to zero
        inxYO   = find(sum(abs(Beta),2) ~= 0);
        exdY    = find(Beta(:,eq) == 0);
        exdY    = intersect(exdY, inxYO);
        
        % Ignore the index of (Lambda) with row elements all equal to zero
        inxLYO  = find(sum(abs(Lambda)) ~= 0);
        exdLY   = find(Lambda(:,eq) == 0);
        exdLY   = intersect(exdLY, inxLYO);
        
        % Coefficients(True value) of estimated equation
        delta1      = -1*[Beta(indY, eq); Lambda(indLY, eq); C(indX, eq)];
        delta1Len   = length(delta1);
        
        K = length([indX; exdX]);
        G = length([indY; exdY]) + 1;
        J = length([indLY; exdLY]);
        k = length(indX);
        g = length(indY);
        j = length(indLY);
        L = (K + J) - (k + j) - g;
        
      
        %** LAYER 2: Sample Size ******************************************  

        tptr = tListLen;
        T = config.tList(tptr);
        rw = 0;             % runway length
        if constT == 1
            cX = config.X(rw+1:T+rw,1:K);
        else
            cX = config.X(rw+1:T+rw,2:K+1);
        end
        X = cX;
        X1 = cX(:,indX);
        X2 = cX(:,exdX);
        y0 = zeros(1,G);

        %** LAYER 3: Distributions ************************************
        dptr = dListLen;
        dist = config.dList(dptr);
        %** LAYER 4: Simulation ***************************************
        V = zeros(delta1Len,eListLen,replications);
        R = zeros(2,eListLen,replications);
        S11 = NaN(delta1Len,eListLen,replications);
        Var_AS = NaN(delta1Len,eListLen,replications);
        Var_BC = NaN(delta1Len,eListLen,replications);
        tic
        for rptr = 1:replications
            [Y LY Vb]      = genY(y0, X, rho, Gamma, Pi, T, G, dist);
            y           = Y(:,eq);
            Y2          = Y;
            Y2(:,eq)    = [];
            LY1         = LY(:,indLY);
            LY1X1       = [LY1 X1];
            LYX         = [LY(:,[indLY exdLY]) X];
            [FirstF(mptr,rptr), Wald, CD] = FirstStageTest(Y2, X2);
            %rejection(mptr,rptr) = boots_test(y, Y2, X, X1, k, T, bst, delta1);  
            for eptr = 1:eListLen
                switch eList(eptr)
                    case 90101;  [V(:,eptr,rptr), R(:,eptr,rptr),S11(1,eptr,rptr), Var_AS(:,eptr,rptr), Var_BC(:,eptr,rptr)] = TSLS_variance_BC(y, Y2, X, X1, k, T, g,L, G,K)   ;
                        %Var_boot2(:,eptr,rptr)          =      TSLS_variance_boot(y, Y2, X, X1, T, g, k, K, bst)  ;
                    case 90102;  [V(:,eptr,rptr), R(:,eptr,rptr),S11(1,eptr,rptr), Var_AS(:,eptr,rptr),Var_BC(:,eptr,rptr)] = kClass_combined_variance_BC(y, Y2, X, X1, 1-1/T^3, 1-1/T, L, L-1,k, T, g,L, G,K);
                        %Var_boot2(:,eptr,rptr)          =      kClass_combined_variance_boot(y, Y2, X, X1, 1-1/T^3, 1-1/T, L, L-1, T, g, k, K, bst);
                    case 90103;  [V(:,eptr,rptr), R(:,eptr,rptr),S11(:,eptr,rptr), Var_AS(:,eptr,rptr), Var_BC(:,eptr,rptr)] = kClass_combined_variance_BC(y, Y2, X, X1, 1-1/T, 1-(L+1)/T, 2, 1,k, T, g,L, G,K);
                        %Var_boot2(:,eptr,rptr)          =      kClass_combined_variance_boot(y, Y2, X, X1, 1-1/T, 1-(L+1)/T, 2, 1, T, g, k, K, bst);
                    case 90104;  [V(:,eptr,rptr), R(:,eptr,rptr),S11(1,eptr,rptr), Var_AS(:,eptr,rptr), Var_BC(:,eptr,rptr)] = kClass_combined_variance_BC(y, Y2, X, X1, 1-1/T^3, 1-(L-1)/T, 2, 1,k, T, g,L, G,K);
                        %Var_boot2(:,eptr,rptr)          =      kClass_combined_variance_boot(y, Y2, X, X1, 1-1/T^3, 1-(L-1)/T, 2, 1, T, g, k, K, bst);
                    case 90105;  [V(:,eptr,rptr), R(:,eptr,rptr),S11(1,eptr,rptr), Var_AS(:,eptr,rptr), Var_BC(:,eptr,rptr)] = fuller_variance_BC(y, Y2, X, X1,k, T, g,L, G,K)   ;
                        %Var_boot2(:,eptr,rptr)          =      fuller_variance_boot(y, Y2, X, X1, T, g, k, K, bst)  ;
                   otherwise
                end
            end
        end

        O(mptr,:,tptr,dptr,:,:) = getInfo(V,R,delta1);
        a = getApproximation(y0, X, X1, Beta, Sigma, Omega, Pi, Gamma, G, g, K, k, J, j, indY, indLY, exdLY, eq, L, T, eList);
        O(mptr,:,tptr,dptr,:,2) = a;
        ahat = reshape(O(mptr,:,tptr,dptr,:,4),eListLen,delta1Len);
        O(mptr,:,tptr,dptr,:,3) = ((a - ahat)./abs(ahat))*100;

        M=NaN(delta1Len,eListLen,5,4);
        ConvP=NaN(delta1Len,eListLen,5);               %convergency probability
        PowerT=NaN(delta1Len,eListLen,5);              % power of t test
        for eptr = 1:eListLen         % For case  {90101 90102 90103 90104 90105}  variance of variance only
            M(:,eptr,:,1)=getstat(Var_AS(:,eptr,:));
            M(:,eptr,:,2)=getstat(Var_BC(:,eptr,:));
            M(:,eptr,:,4)=getstat(S11(:,eptr,:));
            B = reshape(O(mptr,:,tptr,dptr,:,:),eListLen,delta1Len,13);
            [PowerT(:,eptr,1),ConvP(:,eptr,1),Tstats(mptr,:,eptr,1,:)] = conv_power(V(:,eptr,:),Var_AS(:,eptr,:),delta1',T);
            [PowerT(:,eptr,2),ConvP(:,eptr,2),Tstats(mptr,:,eptr,2,:)] = conv_power(V(:,eptr,:),Var_BC(:,eptr,:),delta1',T);
        end
        toc
        printResult_var(mptr, eq, T, replications, bst, constT, dist, L, Beta, Lambda, C, Gamma, Pi, Sigma, Omega, EigVal, B, M, eList,ConvP,PowerT,mean(FirstF(mptr,:)))
        V_All(mptr,:,:,:)=V;
        S11_All(mptr,:,:,:)=S11;
        Var_AS_All(mptr,:,:,:)=Var_AS;
        Var_BC_All(mptr,:,:,:)=Var_BC;
    end
    date_string = datestr(now(), 'yyyymmdd_HHMMSS');
    savefile = ['output' ' (' date_string ').mat'];
    save(savefile, 'O','PowerT','FirstF','V_All','Var_AS_All','Var_BC_All','S11_All','Tstats');
end



function A = getInfo(V,R,delta1)
    eListLen = size(V(:,:,1),2);
    cListLen = size(V(:,:,1),1);
    
    % Total 13 things to be reported.
    A = zeros(eListLen,cListLen, 13); 
    
    for eptr = 1:eListLen
        for cptr = 1:cListLen
            tmp = V(cptr,eptr,:);
            A(eptr,cptr,1)    = delta1(cptr);
            A(eptr,cptr,4)    = mean(tmp) - delta1(cptr);
            A(eptr,cptr,5)    = ((mean(tmp) - delta1(cptr))/abs(delta1(cptr)))*100;
            A(eptr,cptr,6)    = std(tmp);
            A(eptr,cptr,7)    = max(tmp);
            A(eptr,cptr,8)    = min(tmp);
            A(eptr,cptr,9)    = median(tmp) - delta1(cptr);
            A(eptr,cptr,10)    = iqr(tmp);
            A(eptr,cptr,11)   = mse(tmp, delta1(cptr));
            A(eptr,cptr,12)   = mean(R(1,eptr,:));
            A(eptr,cptr,13)   = mean(R(2,eptr,:));
        end
    end
end


function [b, I, V] = OLS(y, Y2, X1, T)
    Y1X1 = [Y2 X1];
    b = Y1X1\y;   
    V = y - Y1X1*b;
    M = eye(T) - ones(T,T)/T;
    I(1) = 1 - (V'*V)/(y'*M*y);
end

function [b, I] = OLS_boot(y, Y2, X, X1, T, bst)
    Y1X1 = [Y2 X1];
    Y1_hat = X*(X\Y2);
    v = Y2 - Y1_hat;
    
    b = Y1X1\y;
    u = y - Y1X1*b;
    
    b_star = zeros(bst, size(b,1));
    I_star = zeros(bst, 2);
    for p = 1:bst
        inx = randi(T,T,1);
        Y1_star = Y1_hat + v(inx, :);
        y_star = [Y1_star X1]*b + u(inx, :);        
        %b_star(p,:) = [Y1_star X1]\y_star;
        [b_star(p,:) I_star(p,:)] = OLS(y_star, Y1_star, X1, T);
    end
    b = 2*b' - mean(b_star);
    I(1) = mean(I_star(1,:));
end

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
        M = eye(T) - ones(T,T)/T;
        %R2 = diag(eye(size(Y2,2)) - (V'*V)./(Y2'*M*Y2));
        I = diag(eye(size(Y2,2)) - ((V'*V)./(T-size(X,2)))./((Y2'*M*Y2)./(T-1)));
        %I(2) = 1-((T-1)/(T-size(X,2)))*(1-R2);
        %I(1) = Y2'*M*Y2;
        %I(2) = V'*V;
        %I(1) = NaN;1 - 
        
    %concentration parameter
    %I(2) = (1/o22) * p22'* X2' * (X2-X1*(X1\X2)) * p22;
end


function [b, I] = kClass_boot(y, Y2, X, X1, k, T, bst)
    [b, dummy] = kClass(y, Y2, X, X1, k, T);
    
    V = Y2 - X*(X\Y2);
    Y1_hat = Y2 - V;
    U = y - [Y1_hat X1]*b;
    
    b_star = zeros(bst, size(b,1));
    I_star = zeros(bst, 2);
    for p = 1:bst
        inx = randi(T,T,1);
        Y1_star = Y1_hat + V(inx, :);
        y_star = [Y1_star X1]*b + U(inx, :);
        [b_star(p,:) I_star(p,:)] = kClass(y_star, Y1_star, X, X1, k, T);
    end
    b = 2*b' - mean(b_star);
    I = mean(I_star(1,:));
    
    %M = eye(T) - ones(T,T)/T;
    %I = diag(eye(size(Y2,2)) - (V'*V)./(Y2'*M*Y2));
       %I(2) = mean(I_star(2,:));
end
        
        
%% Fuller
function [b, I] = fuller(y, Y2, X, X1, T, adj)
    Yd = [y Y2];
    YdtYd = Yd'*Yd;
    %Wsdd = Yd' * Yd - Yd' * X1 * inv(X1' * X1) * X1' * Yd;
    Wsdd = YdtYd - ( Yd' * X1 * (X1\Yd) );
    
    %Wdd = Yd' * Yd - Yd' * X * inv(X' * X) * X' * Yd;
    Wdd = YdtYd - ( Yd' * X * (X\Yd) );
    
    %lambda = min(eig(inv(Wdd) * Wsdd));
    lambda = min(eig(Wdd\Wsdd));
    
    lambda = lambda - adj;
    [b, I] = kClass(y, Y2, X, X1, lambda, T);
end
        
function [b, I] = fuller_boot(y, Y2, X, X1, T, adj, bst)
    [b, dummy] = fuller(y, Y2, X, X1, T, adj);
    
    V = Y2 - X*(X\Y2);
    Y1_hat = Y2 - V;
    U = y - [Y1_hat X1]*b;
    
    b_star = zeros(bst, size(b,1));
    I_star = zeros(bst, 1);
    for p = 1:bst
        inx = randi(T,T,1);
        Y1_star = Y1_hat + V(inx, :);
        y_star = [Y1_star X1]*b + U(inx, :);
        %Yd_star = [y_star Y1_star];
        [b_star(p,:) I_star(p,:)] = fuller(y_star, Y1_star, X, X1, T, adj);
    end
    b = 2*b' - mean(b_star);
    I(1) = mean(I_star(1,:));
end

function rejection = boots_test_2SLS(y, Y2, X, X1, T, bst, delta1)        
   % Boostrap t test
   [b, dummy] = kClass(y, Y2, X, X1, 1, T);
    
    V = Y2 - X*(X\Y2);
    Y1_hat = Y2 - V;
    U = y - [Y1_hat X1]*b;
    
    b_star = zeros(bst, size(b,1));
    I_star = zeros(bst, 2);
    for p = 1:bst
        inx = randi(T,T,1);
        Y1_star = Y1_hat + V(inx, :);
        y_star = [Y1_star X1]*b + U(inx, :);
        [b_star(p,:) I_star(p,:)] = kClass(y_star, Y1_star, X, X1, 1, T);
    end

    if ((b-delta1)< quantile(b_star, 0.025) | (b-delta1)>quantile(b_star, 0.975))
        rejection=1;
    else
         rejection=0;
    end
end

function rejection = boots_test_Fuller(y, Y2, X, X1, T, bst, delta1,adj)        
   % Boostrap t test
   [b, dummy] = fuller(y, Y2, X, X1, T, adj);
    
    V = Y2 - X*(X\Y2);
    Y1_hat = Y2 - V;
    U = y - [Y1_hat X1]*b;
    
    b_star = zeros(bst, size(b,1));
    I_star = zeros(bst, 2);
    for p = 1:bst
        inx = randi(T,T,1);
        Y1_star = Y1_hat + V(inx, :);
        y_star = [Y1_star X1]*b + U(inx, :);
        [b_star(p,:) I_star(p,:)] = fuller(y_star, Y1_star, X, X1, T,adj);
    end

    if ((b-delta1)< quantile(b_star, 0.025) | (b-delta1)>quantile(b_star, 0.975))
        rejection=1;
    else
         rejection=0;
    end
end


function printResult_var(m, eq, T, replications, bst, const, dist, L, Beta, Lambda, C, Gamma, Pi, Sigma, Omega, EigVal, B, M,eList,ConvP,PowerT,FirstF)
    eListLen = length(B(:,1,1));        % model number
    cptrLen = length(B(1,:,1));         % coefficient number
    mList=length(M(1,1,1,:));
    
    date_string = datestr(now(), 'yyyy-mm-dd HHMMSS');
    
    fname = ['M',num2str(m), '_L', num2str(L), '_N',num2str(T), '_P',num2str(Sigma(1,2)), '_','_F',num2str(round(FirstF,1)), ' (',date_string, ')', '.txt'];
    fileID = fopen(fname,'w');
    %fileID = 1;
        
    fprintf(fileID, [date_string, '\r\n']);
    fprintf(fileID,'\r\n');

    fprintf(fileID, 'Elapsed time is %.4f seconds. \r\n',    toc );
    fprintf(fileID,'\r\n');
    
    fprintf(fileID, 'Equation(%d)  N(%d)  R(%d)  Boot(%d)  ConstT(%d)  Dist(%d)  L(%d) \r\n', eq, T, replications, bst, const, dist, L);
    fprintf(fileID, '%s\r\n', '');
    
    fprintf(fileID,'endogenous (structural) = \r\n');
    ftmp = [repmat('%+2.4f   ', 1, size(Beta',2)), '\r\n'];
    fprintf(fileID, ftmp, transpose(Beta'));
    fprintf(fileID,'\r\n');
    
    fprintf(fileID,'lagged endogenous|exgenous (structural) = \r\n');
    ftmp = [repmat('%+2.4f   ', 1, size([Lambda; C]',2)), '\r\n'];
    fprintf(fileID, ftmp, transpose([Lambda; C]'));
    fprintf(fileID,'\r\n');
    
    fprintf(fileID,'lagged endogenous|exgenous (reduced) = \r\n');
    ftmp = [repmat('%+2.4f   ', 1, size([Gamma; Pi]',2)), '\r\n'];
    fprintf(fileID, ftmp, transpose([Gamma; Pi]'));
    fprintf(fileID,'\r\n');
    
    fprintf(fileID,'Sigma = \r\n');
    ftmp = [repmat('%+2.4f   ', 1, size(Sigma,2)), '\r\n'];
    fprintf(fileID, ftmp, transpose(Sigma));
    fprintf(fileID,'\r\n');
    
    fprintf(fileID,'Omega = \r\n');
    ftmp = [repmat('%+10.4f   ', 1, size(Omega,2)), '\r\n'];
    fprintf(fileID, ftmp, transpose(Omega));
    fprintf(fileID,'\r\n');
    
    fprintf(fileID,'Eign roots = \r\n');
    ftmp = [repmat('%+2.4f   ', 1, size(EigVal,2)), '\r\n'];
    fprintf(fileID, ftmp, transpose(EigVal));
    fprintf(fileID,'\r\n');
   
    
    for p = 1:cptrLen
        fprintf(fileID,'%-11s', 'Coefficient');
        fprintf(fileID,'(%1.f)', p);
        fprintf(fileID,'%5s', '');
        %fprintf(fileID,'%-14s', '');
        fprintf(fileID,'%8s', 'True');
        %fprintf(fileID,'%17s', 'Approx.');
        fprintf(fileID,'%18s', 'Bias');
        fprintf(fileID,'%12s', 'Var');
        fprintf(fileID,'%18s', 'Var_Mean');
        fprintf(fileID,'%15s', 'Var_std');
        fprintf(fileID,'%15s', 'Var_Max');
        fprintf(fileID,'%15s', 'Var_Min');
        %fprintf(fileID,'%15s', 'Var_Medium');
        fprintf(fileID,'%17s', 'Covg_Prob');
        fprintf(fileID,'%17s', 'Size_pTest');
        fprintf(fileID,'%17s', 'FistF_test');
        fprintf(fileID,'%8s', 'R2');
        fprintf(fileID,'\r\n');
        fprintf(fileID,'================================================================');
        fprintf(fileID,'================================================================');
        fprintf(fileID,'==========================================================\r\n');
        
       
        meanInx = abs(B(:,p,4));
        [dummy, meanInx] = sort(meanInx);
        [dummy, meanInx] = sort(meanInx);
        
        mseInx = abs(B(:,p,11));
        [dummy, mseInx] = sort(mseInx);
        [dummy, mseInx] = sort(mseInx);
        %}
         
        for q = 1:eListLen
            fprintf(fileID,'%12s', getName(eList(q)));
            fprintf(fileID,'\r\n');
           for m=1:mList 
            fprintf(fileID,'%18s', getNamevar(m));
            fprintf(fileID, '%6.4f', B(q,p,1));
            %fprintf(fileID, '%12.4f',B(q,p,2));
            %fprintf(fileID, ' (%+4.0f%%)',B(q,p,3));
            fprintf(fileID, '     ');
            %fprintf(fileID, '(%2.0f) %6.4f', meanInx(q), B(q,p,4));
            %fprintf(fileID, '%12.4f (%2.0f)', B(q,p,4), meanInx(q));
            fprintf(fileID, '%6.4f',B(q,p,4));
            fprintf(fileID, ' (%+2.0f%%)',B(q,p,5));
            
            fprintf(fileID, '%12.4f', B(q,p,6).^2);
            %fprintf(fileID, '%12.4f', B(q,p,7));
            %fprintf(fileID, '(%2.2f)',(B(q,p,7)/(B(q,p,6).^2)));                     
                fprintf(fileID, '%12.4f', M(p,q,1,m));
                fprintf(fileID, '(%2.2f)',(M(p,q,1,m)/(B(q,p,6).^2)));
                fprintf(fileID, '%14.4f', M(p,q,2,m));
                fprintf(fileID, '%14.4f', M(p,q,3,m));
                fprintf(fileID, '%14.4f', M(p,q,4,m));
                %fprintf(fileID, '%17.4f', M(p,q,5,m));
                
            fprintf(fileID, '%16.4f', ConvP(p,q,m));
            fprintf(fileID, '%16.4f', PowerT(p,q,m));
            fprintf(fileID, '%16.4f', FirstF);
            fprintf(fileID, '%11.4f', B(q,p,13));
            
            fprintf(fileID, '\r\n');
            end;
         end;
                         
        fprintf(fileID, '\r\n');
    end
    fclose(fileID);      
end

function estName = getNamevar(est)
    switch est
        case 1; estName = 'Var_Est            : ';
        case 2; estName = 'Var_BC             : ';
        %case 3; estName = 'Var_Boot1          : ';
        case 3; estName = 'Var_boot2          : ';
        case 4; estName = 'S11                : ';
        otherwise;      
    end
end

function estName = getName(est)
    switch est
            
        case 90101; estName = 'TSLS_Var_Est: ';
        case 90105; estName = 'Fuller1_Var_Est: ';

         otherwise;  
    end
end

function [Y LY Vb] = genY(y0, X, rho, Gamma, Pi, T, G, dist)
    switch dist
        case 1; e = randn(T, G);
        case 2; e = -sqrt(12)/2 + sqrt(12).*rand(T, G);
        otherwise;
    end
    
    Vb = e*rho';
    %Y  = zeros(T, G);
    Y  = zeros(T, G);
    
    if ~any(any(Gamma))
        Y = X*Pi + Vb;
        LY = [];
    else
        %Y(1,:) = y0*Gamma + X(1,:)*Pi + Vb(1,:);
        XV = X*Pi + Vb;
        Y(1,:) = y0*Gamma + XV(1,:);
        
        for t = 2:T
            %Y(t,:) = Y(t-1,:)*Gamma + X(t,:)*Pi + Vb(t,:);
            Y(t,:) = Y(t-1,:)*Gamma + XV(t,:);
        end
        %LY = [y0; Y(1:T-1,:)];
        LY = [y0; Y(1:T-1,:)];
    end
end

function [Yb LYb] = genYb(y0, X, Gamma, Pi, T, G)
    Yb      = zeros(T, G);
    Yb(1,:) = y0*Gamma + X(1,:)*Pi;
    
    for t = 2:T
        Yb(t,:) = Yb(t-1,:)*Gamma + X(t,:)*Pi;
    end
    LYb = [y0; Yb(1:T-1,:)];
end

function b = iqr(X)
    XS  = sort(X);
    N   = length(X);

    q1      = (N+1) / 4;
    q1L     = floor(q1);
    q1R     = q1L + 1;
    qDiff   = q1 - q1L;
    Q1      = XS(q1L) + qDiff * ( XS(q1R) - XS(q1L) );

    q3      = 3*(N+1) / 4;
    q3L     = floor(q3);  
    q3R     = q3L + 1;
    qDiff   = q3 - q3L;
    Q3      = XS(q3L) + qDiff * ( XS(q3R) - XS(q3L) );

    b = Q3 - Q1;
end

function b = mse(X, trueValue)
    tmp = squeeze(X);
    n = length(tmp);
    b = ((tmp - trueValue)' * (tmp - trueValue))/n;
end



function [b, R, Var_AS] = TSLS_variance_AS(y, Y2, X, X1, k, T, g)                
% Estimated Variance 
   [b, R] = kClass(y, Y2, X, X1, 1, T);
    
    Yb = X*(X\Y2);
    Q = [Yb'*Yb Yb'*X1; X1'*Yb X1'*X1];
    Z1 = [Y2 X1];
    e =     y - Z1*b;
    S11 =   e'*e/(T- g - k);
    I = eye(size(Q));
    Var_AS = diag(S11*(I/Q));
 end



function [b, R, S11, Var_AS,Var_BC] = TSLS_variance_BC(y, Y2, X, X1, k, T, g,L, G,K)          
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
    % Positivity safeguard: revert just the offending coordinate to Var_AS, matching the
    % per-coordinate revert used in the heteroskedastic 2SLS case (TSLS_variance_BC_HC).
    nonpos = (Var_BC <= 0);
    Var_BC(nonpos) = Var_AS(nonpos);
end

function Var_boot =  TSLS_variance_boot(y, Y2, X, X1, T, g, k, K, bst)   
        
    [b, R] = kClass(y, Y2, X, X1, 1, T); 
        
    V = Y2 - X*(X\Y2);
    Y2_hat = Y2 - V;
    U = y - [Y2 X1]*b;
    V=(T/(T-K))^0.5*V;
    U=(T/(T-k-g))^0.5*U;
         
    b_star = zeros(bst, size(b,1));
    Var_star = zeros(bst, size(b,1));
    I_star = zeros(bst, 2);
    for p = 1:bst
        inx = randi(T,T,1);
        Y2_star = Y2_hat + V(inx, :);
        y_star = [Y2_star X1]*b + U(inx, :);
        [b_star(p,:), I_star(p,:)]= kClass(y_star, Y1_star, X, X1, 1, T);
        Var_star(p,:)= diag((b_star(p,:)-b')'*( b_star(p,:)-b'));
    end
    Var_boot = sum(Var_star)/(bst);
    
end 
 
%%    fuller_vairance Variance estimator   See Phillips and Xu (2014)    %%  


function [b, R, Var_AS] = fuller_variance_AS(y, Y2, X, X1, k, T, g,K)                
% Estimated Variance 

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
end


function [b, R,S11,Var_AS,Var_BC] = fuller_variance_BC(y, Y2, X, X1,k, T, g,L, G,K)          
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
    Tau = V2'*(v1-V2*b(1:g))/(T);               % e=v1-V2*b(1:g)
    C               = zeros( G-1+k, G-1+k );
    C(1:G-1,1:G-1)  = V2'*V2/(T);
    C1               = zeros( G-1+k, G-1+k );
    C1(1:G-1,1:G-1)  =(Tau*Tau')/S11;
    C2              = C - C1;
   
    Var_bias =  S11*(6*(Q\C1/Q) + trace(Q\C)*(I/Q) - (L-1)*(Q\C2/Q));
    Var_BC = Var_AS - diag(Var_bias);
    % Positivity safeguard: revert just the offending coordinate to Var_AS, matching the
    % per-coordinate revert used in the heteroskedastic 2SLS case (TSLS_variance_BC_HC).
    nonpos = (Var_BC <= 0);
    Var_BC(nonpos) = Var_AS(nonpos);
    
end

function Var_boot = fuller_variance_boot(y, Y2, X, X1, T, g, k, K, bst)   
        
    [b, R] = fuller(y, Y2, X, X1, T, 1/(T-K));   
    V = Y2 - X*(X\Y2);
    Y1_hat = Y2 - V;
    U = y - [Y2 X1]*b;
    V=(T/(T-K))^0.5*V;
    U=(T/(T-k-g))^0.5*U;
         
    b_star = zeros(bst, size(b,1));
    Var_star = zeros(bst, size(b,1));
    I_star = zeros(bst, 2);
    for p = 1:bst
        inx = randi(T,T,1);
        Y1_star = Y1_hat + V(inx, :);
        y_star = [Y1_star X1]*b + U(inx, :);
        [b_star(p,:), I_star(p,:)]= fuller(y_star, Y1_star, X, X1, T, 1/(T-K));
        Var_star(p,:)= diag((b_star(p,:)-b')'*( b_star(p,:)-b'));
    end
    Var_boot = sum(Var_star)/(bst);
    
end