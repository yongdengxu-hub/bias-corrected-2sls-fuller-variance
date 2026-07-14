% run_size_rho.m
% Generate finite-sample size MC outputs for the v9 supplement size tables.
% Runs BOTH estimators (2SLS=90101, Fuller=90105) at rho in {0.279, 0.8},
% for the beta design (endogenous coef tested, col 1) and the gamma design
% (exogenous coef tested, col 3). Homoskedastic DGP (genY_Hetro default).
% Outputs saved to Size_<design>_rho<tag>/ subfolders (timestamped .mat).
clear;
base = 'C:/Users/xuyon/IV_Variance_Correction/matlab/t_test_Hetro';
addpath(base);
S = load([base '/X2.mat']); X2 = S.X2;

CP  = [8 10 15 20 25 70]-1;     % concentration grid (labels 8,10,15,(20),25,70)
k2s = [1 3 5];                  % L = 0, 2, 4
R   = 20000;                    % replications (uniform across the table)

designs = struct('name',  {'beta','beta','gama','gama'}, ...
                 'rho',   {0.279, 0.8,  0.279, 0.8}, ...
                 'isbeta',{true,  true, false, false});

for di = 1:numel(designs)
    D = designs(di);
    rhotag = strrep(sprintf('%.3g', D.rho), '.', 'p');     % 0p279 / 0p8
    outdir = sprintf('%s/Size_%s_rho%s', base, D.name, rhotag);
    if ~exist(outdir,'dir'); mkdir(outdir); end

    cfg = struct();
    cfg.eq = 1; cfg.replications = R; cfg.bst = 0; cfg.dList = 1; cfg.constT = 1;
    cfg.tList = [200]; cfg.X = X2; cfg.eList = [90101 90105];

    if D.isbeta
        b2 = -1e-10;  x1coef = -0.60;   % test beta (col1, true ~0); gamma fixed at -0.6
    else
        b2 = -0.200;  x1coef = -1e-10;  % test gamma (col3, true ~0); beta fixed at -0.2
    end

    cd(outdir);
    fprintf('=== design %s  rho=%.3f  R=%d  -> %s ===\n', D.name, D.rho, R, outdir);
    for kk = 1:3
        k2 = k2s(kk);
        pic = sqrt(CP/(cfg.tList*k2));
        for j = 1:numel(pic)
            cfg.model(1).Beta   = [1 b2; 0 1]';
            cfg.model(1).Lambda = [0 0; 0 0]';
            cfg.model(1).C      = [-1 x1coef zeros(1,k2); 0 0 pic(j)*ones(1,k2)]';
            cfg.model(1).Sigma  = [1 D.rho; D.rho 1];
            doMonteCarlosSimulation_YX_Hetro(cfg);
        end
    end
end
fprintf('ALL SIZE RUNS COMPLETE\n');
