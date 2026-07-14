% run_homo_all.m
% Homoskedastic 2SLS+Fuller runs for the v9 supplement: confidence-interval
% coverage, variance ratios (incl. rho=0.8), and many-instruments (L=20).
% REQUIRES doMonteCarlosSimulation_YX_Hetro.m switch in HOMO mode
% (cases 90101/90105 calling TSLS_variance_BC / fuller_variance_BC).
clear;
base = 'C:/Users/xuyon/IV_Variance_Correction/matlab/t_test_Hetro';
addpath(base);
S = load([base '/X2.mat']); X2 = S.X2;
CP  = [8 10 15 20 25 70];
k2s = [1 3 5 21];     % L = 0, 2, 4, 20
R   = 20000;

cfg0.eq=1; cfg0.replications=R; cfg0.bst=0; cfg0.dList=1; cfg0.constT=1;
cfg0.tList=[200]; cfg0.X=X2; cfg0.eList=[90101 90105];

% ---- Nonzero design (variance ratios + CI coverage): beta=0.2, gamma=0.6 ----
for rho = [0.279 0.8]
    rt = strrep(sprintf('%.3g',rho),'.','p');
    outdir = sprintf('%s/Homo_nz_rho%s', base, rt);
    if ~exist(outdir,'dir'); mkdir(outdir); end
    cd(outdir);
    fprintf('=== nonzero design rho=%.3f -> %s ===\n', rho, outdir);
    for kk = 1:numel(k2s)
        k2 = k2s(kk); pic = sqrt(CP/(200*k2));
        for j = 1:numel(pic)
            cfg = cfg0;
            cfg.model(1).Beta   = [1 -0.200; 0 1]';
            cfg.model(1).Lambda = [0 0; 0 0]';
            cfg.model(1).C      = [-1 -0.60 zeros(1,k2); 0 0 pic(j)*ones(1,k2)]';
            cfg.model(1).Sigma  = [1 rho; rho 1];
            doMonteCarlosSimulation_YX_Hetro(cfg);
        end
    end
end

% ---- Size designs at L=20 (k2=21), rho=0.279: beta-design and gamma-design ----
rho = 0.279; k2 = 21; pic = sqrt(CP/(200*k2));
dn = {'B','G'}; db = [-1e-10 -0.200]; dx = [-0.60 -1e-10];
for di = 1:2
    outdir = sprintf('%s/Homo_size%s_L20_rho0p279', base, dn{di});
    if ~exist(outdir,'dir'); mkdir(outdir); end
    cd(outdir);
    fprintf('=== size design %s L=20 rho=0.279 -> %s ===\n', dn{di}, outdir);
    for j = 1:numel(pic)
        cfg = cfg0;
        cfg.model(1).Beta   = [1 db(di); 0 1]';
        cfg.model(1).Lambda = [0 0; 0 0]';
        cfg.model(1).C      = [-1 dx(di) zeros(1,k2); 0 0 pic(j)*ones(1,k2)]';
        cfg.model(1).Sigma  = [1 rho; rho 1];
        doMonteCarlosSimulation_YX_Hetro(cfg);
    end
end
fprintf('HOMO ALL RUNS COMPLETE\n');
