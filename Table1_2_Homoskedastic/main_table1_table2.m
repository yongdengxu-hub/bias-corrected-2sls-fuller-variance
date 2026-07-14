% Driver for Table 1 (Table_Est_b) and Table 2 (tab:size_main) of the main paper:
% homoskedastic 2SLS/Fuller Monte Carlo, T=200, rho=0.279, beta=0.2, gamma=0.6,
% CP={8,10,15,25,70}, L={0,2,4}. Uses the current (2026-07-06) per-coordinate-revert
% bias correction in doMonteCarlosSimulation_YX.m. Fixed seed for exact reproducibility
% (the original development runs used rng('shuffle'); this replication driver fixes it).
clear;
rng(20260628);
load X2.mat;

config.eq           = 1;
config.replications = 20000;
config.bst          = 0;
config.dList        = 1;
config.constT       = 1;
config.tList        = 200;
config.X            = X2;
config.X(:,6:8)      = randn(1000, 3);
config.eList        = [90101 90105];

beta1 = -0.2; rho = 0.279;
CP    = [8 10 15 25 70];

for k2 = [1 3 5]
    pii = sqrt(CP/(config.tList*k2));
    for i = 1:length(pii)
        config.model(i).Beta    = [1.000 beta1; 0.00 1.000]';
        config.model(i).Lambda  = [0.000 0.000; 0.000 0.000]';
        config.model(i).C       = [-1.000 -0.60 zeros(1,k2); 0 0 pii(i)*ones(1,k2)]';
        config.model(i).Sigma   = [1.000 rho; rho 1.000];
    end
    doMonteCarlosSimulation_YX(config);
    d = dir('output (*).mat'); [~,ix] = max([d.datenum]);
    movefile(d(ix).name, sprintf('table1_2_L%d.mat', k2-1));
    fprintf('done L=%d\n', k2-1);
end
disp('TABLE1_TABLE2_RUN_DONE');
