function CI_AS_All = doMonteCarlosSimulation_CI(config)
% Monte Carlo simulation with CI coverage analysis (90% and 95% CI coverage)
%
% KEY ADDITIONS vs doMonteCarlosSimulation_YX:
%   - 90% and 95% CI coverage probabilities
%   - Average, median CI length
%   - Separate reporting for conventional vs bias-corrected variance estimators
%
% OUTPUTS:
%   CI_AS_All: struct array with CI coverage for each model
%   Also saves output_CI_YYYYMMDD_HHMMSS.mat
%
% USAGE:
%   cd('C:/Users/xuyon/IV_Variance_Correction/matlab/New_codes');
%   main_CI  % sets up config then calls this

%% Extract config
replications = config.replications;
eq           = config.eq;
bst          = config.bst;
constT       = config.constT;

mListLen = length(config.model);
tListLen = length(config.tList);
dListLen = length(config.dList);
eListLen = length(config.eList);

eList = config.eList;
C      = config.model(1).C;
k1 = nnz(C(:, 1)) + 1;

FirstF = NaN(mListLen, replications);

V_All      = zeros(mListLen, k1, eListLen, replications);
S11_All    = NaN(mListLen, k1, eListLen, replications);
Var_AS_All = NaN(mListLen, k1, eListLen, replications);
Var_BC_All = NaN(mListLen, k1, eListLen, replications);

CI_AS_All = struct('ConvP90', {}, 'ConvP95', {}, 'AvgLen90', {}, 'AvgLen95', {}, ...
    'MedLen90', {}, 'MedLen95', {}, 'MinLen90', {}, 'MinLen95', {}, ...
    'MaxLen90', {}, 'MaxLen95', {});
CI_BC_All = CI_AS_All;

for mptr = 1:mListLen
    Beta   = config.model(mptr).Beta;
    C      = config.model(mptr).C;
    Lambda = config.model(mptr).Lambda;
    Sigma  = config.model(mptr).Sigma;
    Gamma  = -1 * Lambda / Beta;
    Pi     = -1 * C / Beta;
    Omega  = Beta' \ Sigma / Beta;
    EigVal = eig(Gamma);
    rho    = chol(Omega, 'lower');

    indX    = find(C(:, eq) ~= 0);
    indY    = find(Beta(:, eq) ~= 0);
    indY(eq, :) = [];
    indLY   = find(Lambda(:, eq) ~= 0);

    inxX0   = find(sum(abs(C), 2) ~= 0);
    exdX    = find(C(:, eq) == 0);
    exdX    = intersect(exdX, inxX0);

    inxYO   = find(sum(abs(Beta), 2) ~= 0);
    exdY    = find(Beta(:, eq) == 0);
    exdY    = intersect(exdY, inxYO);

    inxLYO  = find(sum(abs(Lambda)) ~= 0);
    exdLY   = find(Lambda(:, eq) == 0);
    exdLY   = intersect(exdLY, inxLYO);

    delta1    = -1 * [Beta(indY, eq); Lambda(indLY, eq); C(indX, eq)];
    delta1Len = length(delta1);

    K = length([indX; exdX]);
    G = length([indY; exdY]) + 1;
    J = length([indLY; exdLY]);
    k = length(indX);
    g = length(indY);
    j = length(indLY);
    L = (K + J) - (k + j) - g;

    tptr = tListLen;
    T    = config.tList(tptr);
    rw   = 0;
    if constT == 1
        cX = config.X(rw+1:T+rw, 1:K);
    else
        cX = config.X(rw+1:T+rw, 2:K+1);
    end
    X  = cX;
    X1 = cX(:, indX);
    X2 = cX(:, exdX);
    y0 = zeros(1, G);

    dptr = dListLen;
    dist = config.dList(dptr);

    V      = zeros(delta1Len, eListLen, replications);
    R      = zeros(2, eListLen, replications);
    S11    = NaN(delta1Len, eListLen, replications);
    Var_AS = NaN(delta1Len, eListLen, replications);
    Var_BC = NaN(delta1Len, eListLen, replications);

    tic;
    for rptr = 1:replications
        [Y, LY, Vb] = genY(y0, X, rho, Gamma, Pi, T, G, dist);
        y   = Y(:, eq);
        Y2  = Y;
        Y2(:, eq) = [];
        LY1   = LY(:, indLY);
        LY1X1 = [LY1 X1];
        LYX   = [LY(:, [indLY exdLY]) X];
        [FirstF(mptr, rptr)] = FirstStageTest(Y2, X2);

        for eptr = 1:eListLen
            switch eList(eptr)
                case 90101
                    [V(:, eptr, rptr), R(:, eptr, rptr), S11(1, eptr, rptr), ...
                        Var_AS(:, eptr, rptr), Var_BC(:, eptr, rptr)] = ...
                        TSLS_variance_BC(y, Y2, X, X1, k, T, g, L, G, K);
                case 90105
                    [V(:, eptr, rptr), R(:, eptr, rptr), S11(1, eptr, rptr), ...
                        Var_AS(:, eptr, rptr), Var_BC(:, eptr, rptr)] = ...
                        fuller_variance_BC(y, Y2, X, X1, k, T, g, L, G, K);
            end
        end
    end
    toc;

    %% ---- CI Coverage Analysis (NEW) ----
    % Squeeze to 2D for CI functions: (delta1Len x replications)
    V_squeeze      = squeeze(V(:, 1, :));
    Var_AS_squeeze = squeeze(Var_AS(:, 1, :));
    Var_BC_squeeze = squeeze(Var_BC(:, 1, :));

    CI_AS = computeCIstats(V_squeeze, Var_AS_squeeze, delta1', T);
    CI_BC = computeCIstats(V_squeeze, Var_BC_squeeze, delta1', T);

    CI_AS_All(mptr) = CI_AS;
    CI_BC_All(mptr) = CI_BC;

    %% ---- Standard statistics (from original) ----
    O(mptr, :, tptr, dptr, :, :) = getInfo(V, R, delta1);
    a = getApproximation(y0, X, X1, Beta, Sigma, Omega, Pi, Gamma, G, g, K, k, J, j, ...
        indY, indLY, exdLY, eq, L, T, eList);
    O(mptr, :, tptr, dptr, :, 2) = a;
    ahat = reshape(O(mptr, :, tptr, dptr, :, 4), eListLen, delta1Len);
    O(mptr, :, tptr, dptr, :, 3) = ((a - ahat) ./ abs(ahat)) * 100;

    M = NaN(delta1Len, eListLen, 5, 4);
    ConvP = NaN(delta1Len, eListLen, 5);
    PowerT = NaN(delta1Len, eListLen, 5);

    for eptr = 1:eListLen
        M(:, eptr, :, 1) = getstat(Var_AS(:, eptr, :));
        M(:, eptr, :, 2) = getstat(Var_BC(:, eptr, :));
        M(:, eptr, :, 4) = getstat(S11(:, eptr, :));
        B = reshape(O(mptr, :, tptr, dptr, :, :), eListLen, delta1Len, 13);
        [PowerT(:, eptr, 1)] = conv_power(V(:, eptr, :), Var_AS(:, eptr, :), delta1', T);
        [PowerT(:, eptr, 2)] = conv_power(V(:, eptr, :), Var_BC(:, eptr, :), delta1', T);
    end

    %% ---- Print enhanced results ----
    printResult_CI(mptr, eq, T, replications, bst, constT, dist, L, Beta, Lambda, C, ...
        Gamma, Pi, Sigma, Omega, EigVal, B, M, eList, ConvP, PowerT, ...
        mean(FirstF(mptr, :)), CI_AS, CI_BC);

    %% ---- Store for saving ----
    V_All(mptr, :, :, :)      = V;
    S11_All(mptr, :, :, :)    = S11;
    Var_AS_All(mptr, :, :, :) = Var_AS;
    Var_BC_All(mptr, :, :, :) = Var_BC;
end

date_string = datestr(now(), 'yyyymmdd_HHMMSS');
savefile = ['output_CI_' date_string '.mat'];
save(savefile, 'O', 'FirstF', 'V_All', 'Var_AS_All', 'Var_BC_All', ...
    'S11_All', 'CI_AS_All', 'CI_BC_All');

    %% ============================================================
    %  LOCAL FUNCTIONS (nested inside doMonteCarlosSimulation_CI)
    % ============================================================

    function S = computeCIstats(V, VAR, B1, T)
    % Compute CI coverage stats: 90% and 95% coverage, CI lengths
        [m, n] = size(V);

        df = T - m;
        CV95_u = tinv(0.975, df);   CV95_l = tinv(0.025, df);
        CV90_u = tinv(0.95,   df);   CV90_l = tinv(0.05,   df);

        count90 = false(m, n);  count95 = false(m, n);
        len90 = zeros(m, n);    len95 = zeros(m, n);

        for rptr = 1:n
            se = sqrt(VAR(:, rptr));
            b = V(:, rptr);

            % 95% CI
            u95 = b + CV95_u * se;  l95 = b - CV95_u * se;
            count95(:, rptr) = (B1' >= l95) & (B1' <= u95);
            len95(:, rptr) = u95 - l95;

            % 90% CI
            u90 = b + CV90_u * se;  l90 = b - CV90_u * se;
            count90(:, rptr) = (B1' >= l90) & (B1' <= u90);
            len90(:, rptr) = u90 - l90;
        end

        S.ConvP90  = mean(count90, 2);
        S.ConvP95  = mean(count95, 2);
        S.AvgLen90 = mean(len90, 2);
        S.AvgLen95 = mean(len95, 2);
        S.MedLen90 = median(len90, 2);
        S.MedLen95 = median(len95, 2);
        S.MinLen90 = min(len90, [], 2);
        S.MinLen95 = min(len95, [], 2);
        S.MaxLen90 = max(len90, [], 2);
        S.MaxLen95 = max(len95, [], 2);
    end  % computeCIstats


    function printResult_CI(m, eq, T, replications, bst, constT, dist, L, Beta, Lambda, C, ...
        Gamma, Pi, Sigma, Omega, EigVal, B, M, eList, ConvP, PowerT, FirstF, CI_AS, CI_BC)
    % Print simulation results including CI coverage

        eListLen = length(B(:, 1, 1));
        cptrLen  = length(B(1, :, 1));

        date_str = datestr(now(), 'yyyy-mm-dd HHMMSS');
        fname = ['CI_M', num2str(m), '_L', num2str(L), '_N', num2str(T), ...
            '_P', num2str(Sigma(1, 2)), '_F', num2str(round(FirstF, 1)), ...
            ' (', date_str, ')', '.txt'];
        fid = fopen(fname, 'w');

        fprintf(fid, '=== Monte Carlo with CI Coverage Analysis ===\r\n');
        fprintf(fid, 'Date: %s\r\n', date_str);
        fprintf(fid, 'N=%d  R=%d  L=%d  F=%.1f\r\n\r\n', ...
            T, replications, L, FirstF);

        % Header
        fprintf(fid, '----------------------------------------------------------\r\n');
        fprintf(fid, '95%% CONFIDENCE INTERVAL COVERAGE\r\n');
        fprintf(fid, '----------------------------------------------------------\r\n');
        fprintf(fid, '%-6s  %-8s  %-12s  %-12s  %-12s  %-12s\r\n', ...
            'Coeff', 'True', 'Cov(AS)', 'Cov(BC)', 'AvgLen(AS)', 'AvgLen(BC)');
        fprintf(fid, '%-6s  %-8s  %-12s  %-12s  %-12s  %-12s\r\n', ...
            '------', '------', '-------', '-------', '---------', '---------');
        for p = 1:cptrLen
            fprintf(fid, '%-6d  %+8.4f  %12.4f  %12.4f  %12.4f  %12.4f\r\n', ...
                p, B(1,p,1), CI_AS.ConvP95(p), CI_BC.ConvP95(p), ...
                CI_AS.AvgLen95(p), CI_BC.AvgLen95(p));
        end
        fprintf(fid, '\r\n');

        fprintf(fid, '----------------------------------------------------------\r\n');
        fprintf(fid, '90%% CONFIDENCE INTERVAL COVERAGE\r\n');
        fprintf(fid, '----------------------------------------------------------\r\n');
        fprintf(fid, '%-6s  %-8s  %-12s  %-12s  %-12s  %-12s\r\n', ...
            'Coeff', 'True', 'Cov(AS)', 'Cov(BC)', 'AvgLen(AS)', 'AvgLen(BC)');
        for p = 1:cptrLen
            fprintf(fid, '%-6d  %+8.4f  %12.4f  %12.4f  %12.4f  %12.4f\r\n', ...
                p, B(1,p,1), CI_AS.ConvP90(p), CI_BC.ConvP90(p), ...
                CI_AS.AvgLen90(p), CI_BC.AvgLen90(p));
        end
        fprintf(fid, '\r\n');

        fprintf(fid, '----------------------------------------------------------\r\n');
        fprintf(fid, 'MEDIAN CI LENGTH COMPARISON\r\n');
        fprintf(fid, '----------------------------------------------------------\r\n');
        fprintf(fid, '%-6s  %-8s  %-12s  %-12s  %-12s  %-12s\r\n', ...
            'Coeff', 'True', 'Med95(AS)', 'Med95(BC)', 'Med90(AS)', 'Med90(BC)');
        for p = 1:cptrLen
            fprintf(fid, '%-6d  %+8.4f  %12.4f  %12.4f  %12.4f  %12.4f\r\n', ...
                p, B(1,p,1), CI_AS.MedLen95(p), CI_BC.MedLen95(p), ...
                CI_AS.MedLen90(p), CI_BC.MedLen90(p));
        end
        fprintf(fid, '\r\n');

        fprintf(fid, '----------------------------------------------------------\r\n');
        fprintf(fid, 'VARIANCE RATIO (Var_Est / True Variance)\r\n');
        fprintf(fid, '----------------------------------------------------------\r\n');
        fprintf(fid, '%-6s  %-8s  %-12s  %-12s  %-12s\r\n', ...
            'Coeff', 'True', 'Var_AS/Var', 'Var_BC/Var', 'Bias');
        for p = 1:cptrLen
            var_true = B(1,p,6)^2;
            var_as = M(p,1,1,1);  var_bc = M(p,1,2,1);
            if var_true > 0
                fprintf(fid, '%-6d  %+8.4f  %12.4f  %12.4f  %+12.4f\r\n', ...
                    p, B(1,p,1), var_as/var_true, var_bc/var_true, B(1,p,4));
            end
        end
        fprintf(fid, '\r\n');

        fprintf(fid, '----------------------------------------------------------\r\n');
        fprintf(fid, 't-TEST REJECTION RATES (5%% level)\r\n');
        fprintf(fid, '----------------------------------------------------------\r\n');
        fprintf(fid, '%-6s  %-8s  %-12s  %-12s  %-12s\r\n', ...
            'Coeff', 'True', 't-rej(AS)', 't-rej(BC)', 'Bias');
        for p = 1:cptrLen
            fprintf(fid, '%-6d  %+8.4f  %12.4f  %12.4f  %+12.4f\r\n', ...
                p, B(1,p,1), PowerT(p,1,1), PowerT(p,1,2), B(1,p,4));
        end

        fclose(fid);
        fprintf('Results -> %s\n', fname);
    end  % printResult_CI
end  % doMonteCarlosSimulation_CI


%% ================================================================
%  SUPPORT FUNCTIONS (from doMonteCarlosSimulation_YX.m)
%  Copied here to avoid "Incorrect number or types of inputs" errors
% ================================================================

function A = getInfo(V, R, delta1)
    eListLen = size(V(:,:,1),2);
    cListLen = size(V(:,:,1),1);
    A = zeros(eListLen, cListLen, 13);
    for eptr = 1:eListLen
        for cptr = 1:cListLen
            tmp = V(cptr, eptr, :);
            A(eptr, cptr, 1)  = delta1(cptr);
            A(eptr, cptr, 4)  = mean(tmp) - delta1(cptr);
            A(eptr, cptr, 5)  = ((mean(tmp) - delta1(cptr)) / abs(delta1(cptr))) * 100;
            A(eptr, cptr, 6)  = std(tmp);
            A(eptr, cptr, 7)  = max(tmp);
            A(eptr, cptr, 8)  = min(tmp);
            A(eptr, cptr, 9)  = median(tmp) - delta1(cptr);
            A(eptr, cptr, 10) = iqr(tmp);
            A(eptr, cptr, 11) = mse(tmp, delta1(cptr));
            A(eptr, cptr, 12) = mean(R(1, eptr, :));
            A(eptr, cptr, 13) = mean(R(2, eptr, :));
        end
    end
end

function S = getstat(V)
    V = squeeze(V);
    S(:, 1) = mean(V, 2);
    S(:, 2) = std(V, 0, 2);
    S(:, 3) = max(V, [], 2);
    S(:, 4) = min(V, [], 2);
    S(:, 5) = median(V, 2);
end

function [PowerT, ConvP] = conv_power(V1, VAR, B1, T)
    [m, n] = size(V1);
    CV = tinv(0.975, T - m);
    CV_u = tinv(0.975, T - m);
    CV_l = tinv(0.025, T - m);
    count2 = zeros(m, n);
    t = zeros(m, n);
    if isempty(VAR)
        VAR = var(squeeze(V1)')';
        VAR = VAR' * ones(1, n);
    end
    for rptr = 1:n
        se = sqrt(VAR(:, rptr));
        t(:, rptr) = (V1(:, rptr) - B1') ./ se;
        count2(:, rptr) = t(:, rptr) > CV_u;
    end
    ConvP = mean(count2, 2);
    PowerT = mean(t < CV_l'*ones(1,n), 2) + mean(t > CV_u'*ones(1,n), 2);
end

function [Y LY Vb] = genY(y0, X, rho, Gamma, Pi, T, G, dist)
    switch dist
        case 1; e = randn(T, G);
        case 2; e = -sqrt(12)/2 + sqrt(12).*rand(T, G);
        case 3                                   % standardized t_5 (fat tails), base-MATLAB only
            W = sum(randn(T, G, 5).^2, 3);
            e = (randn(T, G) ./ sqrt(W/5)) / sqrt(5/3);
        case 4                                   % standardized chi^2_4 (skewed), base-MATLAB only
            e = (sum(randn(T, G, 4).^2, 3) - 4) / sqrt(8);
    end
    Vb = e * rho';
    Y = zeros(T, G);
    if ~any(any(Gamma))
        Y = X * Pi + Vb;
        LY = [];
    else
        XV = X * Pi + Vb;
        Y(1, :) = y0 * Gamma + XV(1, :);
        for t = 2:T
            Y(t, :) = Y(t-1, :) * Gamma + XV(t, :);
        end
        LY = [y0; Y(1:T-1, :)];
    end
end

function b = iqr(X)
    XS = sort(X);  N = length(X);
    q1 = (N+1)/4;  q1L = floor(q1);  q1R = q1L + 1;
    q1Diff = q1 - q1L;
    Q1 = XS(q1L) + q1Diff * (XS(q1R) - XS(q1L));
    q3 = 3*(N+1)/4;  q3L = floor(q3);  q3R = q3L + 1;
    q3Diff = q3 - q3L;
    Q3 = XS(q3L) + q3Diff * (XS(q3R) - XS(q3L));
    b = Q3 - Q1;
end

function b = mse(X, trueValue)
    tmp = squeeze(X);
    n = length(tmp);
    b = ((tmp - trueValue)' * (tmp - trueValue)) / n;
end

function a = getApproximation(y0, X, X1, Beta, Sigma, Omega, Pi, Gamma, G, g, K, k, J, j, indY, indLY, exdLY, eq, L, T, eList)
    % getApproximation returns NaN for eList [90101, 90105] (2SLS and Fuller)
    % Since config.eList = [90101] for our CI simulation, this returns NaN
    eListLen = length(eList);
    a = zeros(eListLen, g + k + j);
    for eptr = 1:eListLen
        a(eptr, :) = NaN;  % eList 90101/90105 have no approximation formula here
    end
end

function [FirstF] = FirstStageTest(Y2, X2)
    [n, g1] = size(Y2);
    [n2, k2] = size(X2);
    if n ~= n2; error('Sample size mismatch'); end
    B = Y2' * X2 * pinv(X2' * X2);
    E = Y2 - X2 * B';
    SSE = E' * E;
    SSE = (SSE + SSE') / 2;
    if g1 == 1
        r2 = 1 - SSE / (Y2' * (eye(n) - ones(n,n)/n) * Y2);
        FirstF = r2 * (n - k2) / k2 / (1 - r2);
    else
        invSSE = pinv(SSE);
        M = invSSE / (invSSE + 1);
        F = (B' * X2' * (eye(n) - ones(n,n)/n) * Y2 * M) / k2;
        F = (F + F') / 2;
        FirstF = n * mean(eig(F));
    end
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

%Var_BC = Var_AS - diag(Var_bias);
%if sum(Var_BC<=0)>=1  %| sum(Var_BC>=Var_AS)>=1
%    Var_BC = Var_AS;
%end

Correction = diag(Var_bias);
% 7. Final Bias Corrected Estimator
if L==0
    threshold = 0.90 * Var_AS;  
    mask = (Correction > 0) & (Correction > threshold);  Correction(mask) = threshold(mask);
elseif L>=2
    threshold = 0.50 * Var_AS; 
    mask = (Correction > 0) & (Correction > threshold);  Correction(mask) = threshold(mask);
end
Var_BC = Var_AS - Correction;


end