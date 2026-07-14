function doMonteCarlosSimulation_YX_Hetro(config)
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
    Fc_All=NaN(mListLen,replications);    % classical first-stage F
    Fhc_All=NaN(mListLen,replications);   % robust first-stage F (F_HC)

    V_All = zeros(mListLen,k1,eListLen,replications);
    S11_All = NaN(mListLen,k1,eListLen,replications);
    Var_AS_All = NaN(mListLen,k1,eListLen,replications);
    Var_BC_All = NaN(mListLen,k1,eListLen,replications);
    Var_BCsr_All = NaN(mListLen,k1,eListLen,replications);   % scalar (CP*,rho*) BC, robust first-stage F
    Var_BCsc_All = NaN(mListLen,k1,eListLen,replications);   % scalar (CP,rho) BC, classical first-stage F
    Var_BCf_All  = NaN(mListLen,k1,eListLen,replications);   % FULL scalar BC (delta + leverage add-back)
    Var_BCfQ_All = NaN(mListLen,k1,eListLen,replications);   % diagnostic: full scalar with Q-based scale


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
        Var_BCsr = NaN(delta1Len,eListLen,replications);
        Var_BCsc = NaN(delta1Len,eListLen,replications);
        Var_BCf  = NaN(delta1Len,eListLen,replications);
        Var_BCfQ = NaN(delta1Len,eListLen,replications);
        tic
        for rptr = 1:replications
            [Y LY Vb]      = genY_Hetro(y0, X, rho, Gamma, Pi, T, G, dist);
            y           = Y(:,eq);
            Y2          = Y;
            Y2(:,eq)    = [];
            LY1         = LY(:,indLY);
            LY1X1       = [LY1 X1];
            LYX         = [LY(:,[indLY exdLY]) X];
            [FirstF1, Wald, CD, FirstF(mptr,rptr)] = FirstStageTest(Y2, X2);
            Fc_All(mptr,rptr)=FirstF1; Fhc_All(mptr,rptr)=Wald;
            %rejection(mptr,rptr) = boots_test(y, Y2, X, X1, k, T, bst, delta1);  
            for eptr = 1:eListLen
                switch eList(eptr)
                    %case 90101;  [V(:,eptr,rptr), R(:,eptr,rptr),S11(1,eptr,rptr), Var_AS(:,eptr,rptr), Var_BC(:,eptr,rptr)] = TSLS_variance_BC(y, Y2, X, X1, k, T, g,L, G,K)   ;
                        %Var_boot2(:,eptr,rptr)          =
                        %TSLS_variance_boot(y, Y2, X, X1, T, g, k, K, bst)  ;   
                   %case 90105;  [V(:,eptr,rptr), R(:,eptr,rptr),S11(1,eptr,rptr), Var_AS(:,eptr,rptr), Var_BC(:,eptr,rptr)] = fuller_variance_BC(y, Y2, X, X1,k, T, g,L, G,K)   ;
                        %Var_boot2(:,eptr,rptr)          =
                        %fuller_variance_boot(y, Y2, X, X1, T, g, k, K, bst)  ;  
                   case 90101;  [V(:,eptr,rptr), R(:,eptr,rptr),S11(1,eptr,rptr), Var_AS(:,eptr,rptr), Var_BC(:,eptr,rptr), Var_BCsr(:,eptr,rptr), Var_BCsc(:,eptr,rptr), Var_BCf(:,eptr,rptr), Var_BCfQ(:,eptr,rptr)] = TSLS_variance_BC_HC(y, Y2, X, X1, k, T, g, L, G, K, FirstF1, Wald) ;
                   case 90105;  [V(:,eptr,rptr), R(:,eptr,rptr),S11(1,eptr,rptr), Var_AS(:,eptr,rptr), Var_BC(:,eptr,rptr)] = fuller_variance_BC_HC(y, Y2, X, X1,k, T, g,L, G,K, Wald)   ;

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
       % printResult_var(mptr, eq, T, replications, bst, constT, dist, L, Beta, Lambda, C, Gamma, Pi, Sigma, Omega, EigVal, B, M, eList,ConvP,PowerT,mean(FirstF(mptr,:)))
        V_All(mptr,:,:,:)=V;
        S11_All(mptr,:,:,:)=S11;
        Var_AS_All(mptr,:,:,:)=Var_AS;
        Var_BC_All(mptr,:,:,:)=Var_BC;
        Var_BCsr_All(mptr,:,:,:)=Var_BCsr;
        Var_BCsc_All(mptr,:,:,:)=Var_BCsc;
        Var_BCf_All(mptr,:,:,:)=Var_BCf;
        Var_BCfQ_All(mptr,:,:,:)=Var_BCfQ;
    end
    date_string = datestr(now(), 'yyyymmdd_HHMMSS');
    savefile = ['output' ' (' date_string ').mat'];
    save(savefile, 'O','PowerT','FirstF','Fc_All','Fhc_All','V_All','Var_AS_All','Var_BC_All','Var_BCsr_All','Var_BCsc_All','Var_BCf_All','Var_BCfQ_All','S11_All','Tstats');
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
        case 90106; estName = 'TSLE_Var_HC: ';
        case 90107; estName = 'Fuller1_Var_HC: ';    

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
    Zb=[Yb X1];
    H = (Zb' * Zb) \Zb';
    %Var_AS = ((T-1)/(T-k-g)) * diag(H * diag(e.^2) * H');


    V2 = Y2 - Yb;
    v1=y-X*(X\y);
    Tau = V2'*(v1-V2*b(1:g))/(T-K);
    C               = zeros( G-1+k, G-1+k );
    C(1:G-1,1:G-1)  = V2'*V2/(T-K);
    C1               = zeros( G-1+k, G-1+k );
    C1(1:G-1,1:G-1)  =(Tau*Tau')/S11;
   
    Var_bias = S11*((L+1)*(Q\C1/Q) + trace(Q\C)*(I/Q));
    
    Var_BC = Var_AS - diag(Var_bias);
    if sum(Var_BC<=0)>=1  %| sum(Var_BC>=Var_AS)>=1
       Var_BC = Var_AS;
    end
end



%%    fuller_vairance Variance estimator   See Phillips and Xu (2014)    %%  

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
   
    Zb=[Y2-sqrt(lambda)*V2 X1];
    %Zb=[X*(X\Y2) X1];
    H = Q1 \Zb';
    %Var_HC = ((T-1)/(T-k-g)) * diag(H * diag(e.^2) * H');

    %Var_AS =  ((T-1)/(T-k-g))*diag(inv(Zb'*Zb)* Zb' * diag(e.^2) * Zb * inv(Zb'*Zb));


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
    if sum(Var_BC<=0)>=1  %| sum(Var_BC>=Var_AS)>=1
       Var_BC = Var_AS;
    end
    
end

function [b, R, S11, Var_HC, Var_BC, Var_BCsr, Var_BCsc, Var_BCf, Var_BCfQ] = TSLS_variance_BC_HC(y, Y2, X, X1, k, T, g, L, G, K, Fc, Fhc)
% TSLS_variance_BC_HC: Heteroskedasticity-robust bias corrected variance estimator
% Based on Theorem 3 of the provided theory.
    global DISABLE_LEV

    % 1. Standard 2SLS Estimation (k=1)
    [b, R] = kClass(y, Y2, X, X1, 1, T); 
    
    % 2. Construct Q_inv (The (Z'Z)^-1 term)
    % Note: In your notation Q is the inverse. In code, 'Q' is usually the matrix itself.
    % We construct the projected instruments Z_hat = [Yb X1]
    Yb = X * (X \ Y2);
    Z_hat = [Yb X1];
    Q_mat = Z_hat' * Z_hat;  % This is \hat{Z}_1' \hat{Z}_1
    Q = inv(Q_mat);          % This corresponds to Q in the theorem
    
    % 3. Residuals and Heteroskedasticity Statistics
    Z1 = [Y2 X1];
    e = y - Z1 * b;          % 2SLS Residuals
    e2 = e.^2;
    sig2_bar = mean(e2);     % Average variance (T^-1 tr(Omega))
    m4 = mean(e.^4);         % Fourth moment (Proxy for T^-1 tr(Omega^2))
    
    % Calculate delta (Assuming normality for the kurtosis adjustment factor of 3)
    % delta = (mean(sigma^2))^2 / mean(sigma^4). 
    % Under normality assumption: mean(u^4) approx 3 * mean(sigma^4)
    delta_hat = (3 * sig2_bar^2) / m4; 
    
    % Cap delta at 1 as per definition
    if delta_hat > 1; delta_hat = 1; end

    % 4. Standard HC Variance Estimator (Sandwich)
    % Var_HC = Q * (Z_hat' * Omega * Z_hat) * Q
    Omega = diag(e2);
    Meat = Z_hat' * (e2 .* Z_hat); % Efficient calculation of Z' * diag(e^2) * Z
    Var_HC = Q * Meat * Q;
    
    % 5. Construct Reduced Form Components for Bias Correction
    % V2: Reduced form residuals of Y2
    V2 = Y2 - Yb; 
    
    % pi_hat: Scaled correlation between V2 and u
    % formula: (1/T * V2' * u) / sig2_bar
    pi_hat = (V2' * e) / T / sig2_bar;
    
    % Construct C_tilde matrices
    % C_tilde_1 (Upper left block G-1 x G-1)
    % Formula: (m4 / 3) * pi * pi'
    C1_upper = (m4 / 3) * (pi_hat * pi_hat');
    
    % C_tilde_2 (Upper left block)
    % Formula: T^-1 * (V2' * Omega * V2) - C_tilde_1
    V2_Omega_V2 = V2' * (e2 .* V2); % V2' * Omega * V2
    C2_upper = (V2_Omega_V2 / T) - C1_upper;
    
    % Embed into full (G-1+k) matrices
    C_tilde_1 = zeros(G-1+k, G-1+k);
    C_tilde_1(1:G-1, 1:G-1) = C1_upper;
    
    C_tilde_2 = zeros(G-1+k, G-1+k);
    C_tilde_2(1:G-1, 1:G-1) = C2_upper;
    
    C_tilde = C_tilde_1 + C_tilde_2;
    
    % 6. Robust HC0 bias (corrected eq:HC_bias, delta-bookkeeping fixed 2026-07-03;
    %    the old Coeff2 = delta*(L-1)^2-(L^2-3L) dropped delta from the pi-pi' channel):
    %    tr(Q*C_tilde_2)*Q + delta*tr(Q*C_tilde_1)*Q + delta*(L+1)*Q*C_tilde_1*Q
    Term1  = trace(Q * C_tilde_2) * Q + delta_hat * trace(Q * C_tilde_1) * Q;
    Coeff2 = delta_hat * (L + 1);
    Term2 = Coeff2 * (Q * C_tilde_1 * Q);
    Var_bias_conv = Term1 + Term2;

    % 6b. Leverage / residual-estimation correction Delta_lev  (v8: Lemma 4 / App. D.7).
    %    E[V_HC] = (conventional expectation) - Delta_lev; corrected bias = (conv bias) - Delta_lev,
    %    so  V_HC^BC = V_HC - (Var_bias_conv - Delta_lev).  Estimable HC2-flavoured per-observation
    %    form (exact under homoskedasticity, validated in the semi-weak regime mu^2 ~ 15-25):
    %       Delta_lev = Q * Z_hat' * diag( e_t^2 * h_t/(1-h_t) ) * Z_hat * Q,   h_t = zhat_t' Q zhat_t.
    h_t       = sum((Z_hat * Q) .* Z_hat, 2);  % leverage h_t = zhat_t' Q zhat_t (diag only, fast)
    w_lev     = e2 .* h_t ./ (1 - h_t);        % per-observation leverage weight
    Delta_lev = Q * (Z_hat' * (w_lev .* Z_hat)) * Q;
    if ~isempty(DISABLE_LEV) && DISABLE_LEV; Delta_lev(:) = 0; end

    % 7. Corrected bias-corrected estimator:  V_HC - (conventional bias - Delta_lev)
    Var_bias = Var_bias_conv - Delta_lev;
    dHC      = diag(Var_HC);
    Var_BC   = dHC - diag(Var_bias);

    % Positivity safeguard: at L=0 the 2SLS variance has no finite moments, so a
    % corrected coordinate can turn non-positive; revert just that coordinate to HC.
    nonpos         = (Var_BC <= 0);
    Var_BC(nonpos) = dHC(nonpos);

    % --- New scalar (CP,rho) / (CP*,rho*) bias correction for beta (1st coefficient) ---
    %   V_HC^BC(beta) = V_HC(beta) - sig2_bar*(1+(L+1)*rho^2)/CP^2.
    %   CP  = kz*F (classical first-stage F)  -> Var_BCsc  (classical)
    %   CP* = kz*F_HC (robust first-stage F)  -> Var_BCsr  (robust, the new formula)
    %   rho = corr(structural resid e, first-stage resid V2)  [= rho* under common het].
    kz       = L + 1;
    rho_hat  = (V2' * e) / sqrt( (V2'*V2) * (e'*e) );
    CP_class = kz * max(Fc  - 1, 1e-6);          % E(F) ~ 1 + CP/kz  =>  CP = kz(F-1)
    CP_star  = kz * max(Fhc - 1, 1e-6);          % robust:  CP* = kz(F* - 1)
    corr_sc  = sig2_bar * (1 + (L+1)*rho_hat^2) / CP_class^2;
    corr_sr  = sig2_bar * (1 + (L+1)*rho_hat^2) / CP_star^2;
    Var_BCsc = dHC;  Var_BCsc(1) = dHC(1) - corr_sc;
    Var_BCsr = dHC;  Var_BCsr(1) = dHC(1) - corr_sr;
    if Var_BCsc(1) <= 0; Var_BCsc(1) = dHC(1); end
    if Var_BCsr(1) <= 0; Var_BCsr(1) = dHC(1); end

    % --- FULL scalar BC (closest to the matrix form): scalar het conventional
    %     coefficient over CP*^2, PLUS the actual leverage add-back.
    %     Corrected 2026-07-03: substituting the fixed matrix-form bias through the
    %     (rho^2/delta,CP^2) dictionary makes delta cancel exactly, giving the clean
    %     delta-free 1+(L+1)rho^2 -- matching Var_BCsc/Var_BCsr below, which were
    %     already using this correct form.
    conv_coef = 1 + (L+1)*rho_hat^2;
    % FULL scalar BC in the CONCENTRATION PARAMETER: uses the robust first-stage estimate
    % CP* = kz(F*-1). Exact given the concentration; CP* is its practical robust estimate.
    corr_f    = sig2_bar * conv_coef / CP_star^2 - Delta_lev(1,1);
    Var_BCf   = dHC;  Var_BCf(1) = dHC(1) - corr_f;
    if ~isfinite(Var_BCf(1)) || Var_BCf(1) <= 0 || Var_BCf(1) > 3*dHC(1)
        Var_BCf(1) = dHC(1);   % revert to HC when the formula is unstable (per-rep CP*->0 / delta->0)
    end

    % Reference only: same correction with the EXACT realized concentration 1/Q(1,1)
    % (equals the matrix-form BC; used to isolate the CP*-estimation gap).
    corr_fQ   = sig2_bar * conv_coef * Q(1,1)^2 - Delta_lev(1,1);
    Var_BCfQ  = dHC;  Var_BCfQ(1) = dHC(1) - corr_fQ;
    if ~isfinite(Var_BCfQ(1)) || Var_BCfQ(1) <= 0 || Var_BCfQ(1) > 3*dHC(1)
        Var_BCfQ(1) = dHC(1);
    end

    % Output as column vectors (compatible with previous signature)
    Var_HC = dHC;
    S11    = sig2_bar;                          % average residual variance
end

function [b, R, S11, Var_HC, Var_BC] = fuller_variance_BC_HC(y, Y2, X, X1, k, T, g, L, G, K, Fhc)
% fuller_variance_BC_HC: Heteroskedasticity-robust bias corrected variance estimator for Fuller
% Based on Theorem 4 of the provided theory.
    global DISABLE_LEV DIAG_FULLER

    % 1. Fuller Estimation
    Yd = [y Y2];
    YdtYd = Yd'*Yd;
    Wsdd = YdtYd - ( Yd' * X1 * (X1\Yd) );
    Wdd = YdtYd - ( Yd' * X * (X\Yd) );
    lambda = min(eig(Wdd\Wsdd));
    
    % Fuller constant c=1
    k_fuller = lambda - 1/(T-K); 
    
    [b, R] = kClass(y, Y2, X, X1, k_fuller, T);  
    
    % 2. Residuals and Heteroskedasticity Statistics
    Z1 = [Y2 X1];
    e = y - Z1 * b;          % Fuller Residuals
    e2 = e.^2;
    sig2_bar = mean(e2);
    m4 = mean(e.^4);
    
    % Calculate delta (with normality assumption factor 3)
    delta_hat = (3 * sig2_bar^2) / m4;
    if delta_hat > 1; delta_hat = 1; end
    
    % 3. Construct Q matrix (2SLS projection, used for the CORRECTION terms only;
    %    the theory's Q = (Z_bar'Z_bar)^-1 is approximated by the projected Z).
    Yb = X * (X \ Y2);
    Z_hat = [Yb X1];
    Q_mat = Z_hat' * Z_hat;
    Q = inv(Q_mat);

    % 4. HC Variance Estimator: PROPER k-class sandwich (fixed 2026-07-03).
    %    The paper's eq:HC_est with k = k_F and the Supplement's eqG:identity define
    %    V_HC(alpha_F) = A_F^-1 * Z1F' diag(uF^2) Z1F * A_F^-1, with pseudo-instrument
    %    Z1F = Z1 - k_F*M_X*Z1 = [Yhat2 + (1-k_F)*Vhat2 : X1] and bread A_F = Z1F'Z1.
    %    The previous code used a hybrid (2SLS bread Q with meat on Z_hat), which is a
    %    DIFFERENT estimator: its finite-sample bias is negative at L>=4 and Theorem 4
    %    does not apply to it (verified by MC, c4_fuller_hybrid.m). Theorem 4's own-
    %    observation cancellation lives in the (1-k_F)*Vhat2 part of the meat below.
    V2   = Y2 - Yb;                    % first-stage residuals Vhat2 = M_X*Y2
    Z1F  = [Yb + (1 - k_fuller) * V2, X1];
    A_F  = Z1F' * Z1;                  % = Z1'(I - k_F M_X) Z1, symmetric
    MeatF  = Z1F' * (e2 .* Z1F);
    AFinv  = inv(A_F);
    Var_HC = AFinv * MeatF * AFinv';

    % 5. Construct Reduced Form Components (V2 already formed above)
    pi_hat = (V2' * e) / T / sig2_bar;

    C1_upper = (m4 / 3) * (pi_hat * pi_hat');
    V2_Omega_V2 = V2' * (e2 .* V2);
    C2_upper = (V2_Omega_V2 / T) - C1_upper;

    C_tilde_1 = zeros(G-1+k, G-1+k);
    C_tilde_1(1:G-1, 1:G-1) = C1_upper;

    C_tilde_2 = zeros(G-1+k, G-1+k);
    C_tilde_2(1:G-1, 1:G-1) = C2_upper;

    C_tilde = C_tilde_1 + C_tilde_2;

    % 6. Bias-corrected Fuller variance: CONVENTIONAL k-class form (decision 2026-07-03,
    %    fuller_decomp.m). In the semi-weak band (CP 8-25) the proper HC0 sandwich is
    %    dominated by bread-Jensen inflation from A_F^-1 * (.) * A_F^-1 (0.6-1.8 x Var,
    %    a distributional effect no O(T^-2) meat correction can remove; the own/cross
    %    pseudo-instrument meat terms are only 1-7% of Var). The conventional k-class
    %    variance sig2F * A_F^-1, with one power of A_F^-1, is near-unbiased there
    %    (ratios 0.97-1.17 at L=4), and its verified conventional bias formula
    %    (Supplement S.3.3, MC-checked) is
    %       bias_conv = tr(Q*C2)Q + delta*tr(Q*C1)Q + 6*delta*Q*C1*Q - (L-1)*Q*C2*Q,
    %    valid under Assumption 1 (balance). The robust Theorem-4 sandwich correction
    %    remains valid at strong instruments; Var_HC below still reports the proper
    %    sandwich for the V_HC/Var column.
    p_dim    = G - 1 + k;
    sig2F    = (e' * e) / (T - p_dim);
    Var_conv = sig2F * AFinv;

    bias_conv = trace(Q * C_tilde_2) * Q + delta_hat * trace(Q * C_tilde_1) * Q ...
              + 6 * delta_hat * (Q * C_tilde_1 * Q) - (L - 1) * (Q * C_tilde_2 * Q);

    % Delta_lev kept for diagnostics only (not part of the conventional correction)
    h_t       = sum((Z_hat * Q) .* Z_hat, 2);
    w_lev     = e2 .* h_t ./ (1 - h_t);
    Delta_lev = Q * (Z_hat' * (w_lev .* Z_hat)) * Q;
    if ~isempty(DISABLE_LEV) && DISABLE_LEV; Delta_lev(:) = 0; end

    % 7. Two-branch gated multiplicative correction (cpr:posBC form; decision 2026-07-04,
    %    fuller_L0_fix.m). The branches follow the structural fact E[1-lambda] = -L/(T-K):
    %    L = 0 (just-identified): the k-class bread equals the 2SLS bread up to mean-zero
    %      noise, so the sandwich has NO bread-Jensen inflation, and its meat-bread tail
    %      cancellation is essential precisely because the single-dof bread is heavy-tailed
    %      (the conventional form fails here: ratios 1.9/1.6/1.3 at CP=8/10/15). Use the
    %      SANDWICH with the Theorem-4 bias and a PROPORTIONAL gate w = min(CP*/10,1),
    %      which engages from CP*=0; MC: BC/Var = 0.998/0.954/0.933/0.937 at CP=8/10/15/25.
    %    L >= 1: bread-Jensen inflation appears (sandwich raw 1.59 already at L=1, CP=8) and
    %      the conventional k-class form with the standard gate is the right tool
    %      (fuller_decomp.m; L=4: BC/Var = 1.18/1.16/1.11/1.02/1.02).
    kz       = L + 1;
    CP_star  = kz * max(Fhc - 1, 1e-6);
    if L == 0
        h_tF   = sum((Z_hat * Q) .* Z_hat, 2);
        bias_thm4 = trace(Q * C_tilde_2) * Q + delta_hat * trace(Q * C_tilde_1) * Q ...
                  + delta_hat * (L + 5) * (Q * C_tilde_1 * Q) - Delta_lev;
        w = min(max(CP_star / 10, 0), 1);
        dHC0      = diag(Var_HC);
        bias_diag = diag(bias_thm4);
        Var_BC    = dHC0 ./ (1 + w * bias_diag ./ dHC0);
        dConv = diag(Var_conv);   % kept for diagnostics
    else
        % L >= 1 (decision 2026-07-06): additive correction with per-coordinate
        % revert-if-negative, the SAME rule as the paper's other three cases (homo 2SLS,
        % homo Fuller, hetero 2SLS). Same-seed comparison against the gated multiplicative
        % form showed the two are close on net (L=2 improves: beta CP=8 1.36->1.11;
        % L=4 worsens: beta CP=8 1.18->1.38; zero negative draws in 20k reps at L>=1),
        % so consistency wins here. The gate remains ONLY at L=0, where the Theorem-4
        % bias estimate is wild-tailed (b/V up to 4.7e5) and every revert-style
        % alternative undercorrects (0.80-0.96 vs 0.93-1.00 gated).
        w = NaN;   % kept for the DIAG_FULLER row layout
        dConv     = diag(Var_conv);
        bias_diag = diag(bias_conv);
        Var_BC    = dConv - bias_diag;
        nonpos    = (Var_BC <= 0);
        Var_BC(nonpos) = dConv(nonpos);
    end

    if ~isempty(DIAG_FULLER) && DIAG_FULLER
        row = [delta_hat, CP_star, w, k_fuller, dConv(1), bias_diag(1), Delta_lev(1,1), Var_BC(1), Var_HC(1,1)];
        try
            acc = evalin('base','diag_fuller_acc');
        catch
            acc = zeros(0,9);
        end
        acc(end+1,:) = row;
        assignin('base','diag_fuller_acc', acc);
    end

    Var_HC = diag(Var_HC);
    S11 = sig2_bar;

end