% Table_Est_HC replica at the PAPER'S design: delta ~ 0.78 (HET_B=0.4985), proper Fuller
% sandwich. Produces Bias/Var/ratio_HC/ratio_BC for both beta and gamma, matching the format
% of the main-text table (caption: "balanced design, delta~0.78").
clear;
global HET_TYPE HET_Z HET_B
rng(20260628);
load X2.mat;
HET_TYPE = 1;
HET_B    = 0.4985;                 % delta = exp(-HET_B^2) ~ 0.78
HET_Z    = randn(200,1);
config.eq           = 1;
config.replications = 20000;
config.bst          = 0;
config.dList        = 1;
config.constT       = 1;
config.tList        = 200;
config.X            = X2;
config.X(:,3:50)    = randn(1000,48);
config.eList        = [90101 90105];
beta1 = -0.2; rho = 0.279;
CP1   = [8 10 15 25 70] - 1;
for k2 = [1 3 5]
    pii = sqrt(CP1/(config.tList*k2));
    for i = 1:length(pii)
        config.model(i).Beta   = [ 1.000  beta1 ; 0.00  1.000 ]';
        config.model(i).Lambda = [ 0.000  0.000 ; 0.000 0.000 ]';
        config.model(i).C      = [ -1.000 -0.60  zeros(1,k2) ; 0 0 pii(i)*ones(1,k2) ]';
        config.model(i).Sigma  = [ 1.000  rho ; rho 1.000 ];
    end
    doMonteCarlosSimulation_YX_Hetro(config);
    d = dir('output (*).mat'); [~,ix] = max([d.datenum]);
    movefile(d(ix).name, sprintf('sizeci78_L%d.mat', k2-1));
    fprintf('done L=%d\n', k2-1);
end
v = exp(HET_B*HET_Z); v = v/mean(v);
fprintf('design delta = %.4f\n', (mean(v))^2/mean(v.^2));
disp('SIZECI78_RUN_DONE');
