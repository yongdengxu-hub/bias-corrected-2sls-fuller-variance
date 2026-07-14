
function [F, F_HC, CD, F_Eff] = FirstStageTest(Y2, X2)
    [T,k2] = size(X2);  
    pi = X2\Y2;  %Z=X2;
    v2 = Y2-X2*pi;
    sigmaZ= v2'*v2/(T); 

    F  =(X2*pi)'*inv(sigmaZ)*(X2*pi)/k2;
    Wald = pi'*((sigmaZ)^(-1).*(X2'*X2))*pi/k2;
    CD = (sigmaZ)^(-0.5)*(X2*pi)'*X2*pi* (sigmaZ)^(-0.5)/k2;
    %lm=fitlm(X2,Y2,'Intercept',false);
    %tStats = lm.Coefficients.tStat;
    % Compute standard heteroskedasticity-robust F-statistics
    HC = X2' * diag(v2.^2) * X2;
    F_r= Y2'*X2*inv(HC)*X2'*Y2/k2;
    F_HC= (X2*pi)'*X2*inv(HC)*X2'*(X2*pi)/k2;

   % Compute effective F statistics of Montiel Olea and Pflueger (2013) 
    F0=HC*inv(X2'*X2);    
    % HC is W_2 in  Lewis and Mertens(2022) and Windmeijer (2003) or
    % Sigma_vv in Andrews et al. (2019)

    F_Eff=pi'*(X2'*X2)*pi*inv(trace(F0));
end