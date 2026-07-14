function run_CI_nn(distIdx, Lidx, R, outroot)
% Non-normal-error robustness: homoskedastic size/CI experiment under fat-tailed
% (t_5) and skewed (chi^2_4) errors, for 2SLS (90101) and Fuller (90105). One
% SLURM array task = one (dist, L) block; each task loops over rho in {0.279,0.8}
% and the CP grid (keeps the array <= the htc_genoa 10-job submit cap).
%   distIdx : 1 -> dist=3 (standardized t_5), 2 -> dist=4 (standardized chi^2_4)
%   Lidx    : 1..3 -> L = 0,2,4  (k2 = 1,3,5)
%   R       : replications (default 200000); outroot default '.'
% Requires on path: doMonteCarlosSimulation_CI.m (genY dist cases 3,4),
% fuller_variance_BC.m, kClass.m, X2.mat. Deterministic cell_CP<NN>.mat names.
    if nargin<3||isempty(R), R=200000; end
    if nargin<4||isempty(outroot), outroot='.'; end
    if ischar(distIdx)||isstring(distIdx), distIdx=str2double(distIdx); end
    if ischar(Lidx)||isstring(Lidx), Lidx=str2double(Lidx); end
    if ischar(R)||isstring(R), R=str2double(R); end

    here = fileparts(mfilename('fullpath')); addpath(here);
    Sd = load(fullfile(here,'X2.mat')); X2 = Sd.X2;

    dists = [3 4]; rhos = [0.279 0.8]; k2s = [1 3 5]; CP = [8 10 15 20 25 70];
    dist = dists(distIdx); k2 = k2s(Lidx); Lval = k2-1; k1 = k2+2;
    t0 = tic;
    for ri = 1:numel(rhos)
        rho = rhos(ri); rt = strrep(sprintf('%.3g',rho),'.','p');
        rng(20260627 + 10000*distIdx + 100*ri + Lidx, 'twister');
        outdir = fullfile(outroot,'raw_nn', sprintf('dist%d_rho%s_L%d', dist, rt, Lval));
        if ~exist(outdir,'dir'); mkdir(outdir); end
        X = X2; X(:,5:50) = randn(1000,46);
        pic = sqrt(CP/(200*k2));
        fprintf('=== NN: dist=%d rho=%.3f L=%d k2=%d R=%d -> %s ===\n', dist, rho, Lval, k2, R, outdir);
        for j = 1:numel(pic)
            cfg = struct('eq',1,'replications',R,'bst',0,'dList',dist,'constT',1, ...
                         'tList',200,'X',X,'eList',[90101 90105]);
            Cm = zeros(k1,2); Cm(1,1)=-1; Cm(2,1)=-0.6; Cm(3:k1,2)=pic(j);
            cfg.model(1).Beta=[1 -0.200; 0 1]'; cfg.model(1).Lambda=[0 0;0 0]';
            cfg.model(1).C=Cm; cfg.model(1).Sigma=[1 rho; rho 1];
            oldp = cd(outdir);
            try, doMonteCarlosSimulation_CI(cfg); catch ME, cd(oldp); rethrow(ME); end
            cd(oldp);
            dd = dir(fullfile(outdir,'output_CI_*.mat')); [~,kk] = max([dd.datenum]);
            movefile(fullfile(outdir,dd(kk).name), fullfile(outdir,sprintf('cell_CP%02d.mat',CP(j))));
            fprintf('  rho=%.3f CP=%d done (%.0f s)\n', rho, CP(j), toc(t0));
        end
    end
    fprintf('TASK COMPLETE dist=%d L=%d in %.0f s\n', dist, Lval, toc(t0));
end
