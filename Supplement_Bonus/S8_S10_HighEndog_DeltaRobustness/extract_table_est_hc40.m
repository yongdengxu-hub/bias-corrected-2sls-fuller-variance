% Extract Table_Est_HC-format rows (Bias, Var, ratio_HC, ratio_BC for beta & gamma)
% from sizeci40_L{0,2,4}.mat (delta~0.40, more-severe-heteroskedasticity robustness design).
clear; Ls=[0 2 4]; CPd=[8 10 15 25 70];
coefIdx = [1 3];   % [beta, gamma] column indices into V_All/Var_AS_All/Var_BC_All
truth   = [0.2 0.6];

fprintf('\n\\hline\n\\multicolumn{9}{c}{\\textbf{Panel A: 2SLS}} \\\\\n\\hline\n');
for li=1:3
    S = load(sprintf('sizeci40_L%d.mat', Ls(li)));
    V = S.V_All; A = S.Var_AS_All; B = S.Var_BC_All;
    fprintf('\\multicolumn{9}{c}{$L=%d$} \\\\\n', Ls(li));
    for cp=1:5
        row = sprintf('%d', CPd(cp));
        for est=1
            for ci=1:2
                co = coefIdx(ci);
                vb = squeeze(V(cp,co,est,:)); m = mean(vb,'omitnan'); varMC = var(vb,'omitnan');
                bias = m - truth(ci);
                ehc0 = mean(squeeze(A(cp,co,est,:)),'omitnan');
                ebc  = mean(squeeze(B(cp,co,est,:)),'omitnan');
                row = [row sprintf(' & %.3f & %.3f & %s & %s', bias, varMC, fmt(ehc0/varMC), fmt(ebc/varMC))];
            end
        end
        fprintf('%s \\\\\n', row);
    end
end

fprintf('\\hline\n\\multicolumn{9}{c}{\\textbf{Panel B: Fuller}} \\\\\n\\hline\n');
for li=1:3
    S = load(sprintf('sizeci40_L%d.mat', Ls(li)));
    V = S.V_All; A = S.Var_AS_All; B = S.Var_BC_All;
    fprintf('\\multicolumn{9}{c}{$L=%d$} \\\\\n', Ls(li));
    for cp=1:5
        row = sprintf('%d', CPd(cp));
        for est=2
            for ci=1:2
                co = coefIdx(ci);
                vb = squeeze(V(cp,co,est,:)); m = mean(vb,'omitnan'); varMC = var(vb,'omitnan');
                bias = m - truth(ci);
                ehc0 = mean(squeeze(A(cp,co,est,:)),'omitnan');
                ebc  = mean(squeeze(B(cp,co,est,:)),'omitnan');
                row = [row sprintf(' & %.3f & %.3f & %s & %s', bias, varMC, fmt(ehc0/varMC), fmt(ebc/varMC))];
            end
        end
        fprintf('%s \\\\\n', row);
    end
end
disp('EXTRACT40_DONE');

function s = fmt(x)
    if x > 999
        s = sprintf('%.2E', x);
    else
        s = sprintf('%.3f', x);
    end
end
