% Extract Table 1 (Table_Est_b: Bias/Var/ratio) and Table 2 (tab:size_main: size/CI/length)
% from reverttest_L{0,2,4}.mat, using the CURRENT per-coordinate-revert Var_BC_All
% (matches the fixed doMonteCarlosSimulation_YX.m, decision 2026-07-06).
clear; z = 1.959964;
Ls = [0 2 4]; CPd = [8 10 15 25 70]; estlab = {'2SLS','Fuller'};
coefIdx = [1 3]; truth = [0.2 0.6];

fprintf('=== TABLE 1 (Bias/Var/ratio) ===\n');
for est = 1:2
    fprintf('\n-- %s --\n', estlab{est});
    for li = 1:3
        S = load(sprintf('table1_2_L%d.mat', Ls(li)));
        V = S.V_All; A = S.Var_AS_All; B = S.Var_BC_All;
        fprintf('L=%d\n', Ls(li));
        for cp = 1:5
            row = sprintf('%d', CPd(cp));
            for ci = 1:2
                co = coefIdx(ci);
                vb = squeeze(V(cp,co,est,:)); m = mean(vb,'omitnan'); varMC = var(vb,'omitnan');
                bias = m - truth(ci);
                ehc = mean(squeeze(A(cp,co,est,:)),'omitnan');
                ebc = mean(squeeze(B(cp,co,est,:)),'omitnan');
                row = [row sprintf(' & %.3f & %.3f & %s & %s', bias, varMC, fmt(ehc/varMC), fmt(ebc/varMC))];
            end
            fprintf('%s \\\n', row);
        end
    end
end

fprintf('\n=== TABLE 2 (Size/CI/Length) ===\n');
for est = 1:2
    fprintf('\n-- %s --\n', estlab{est});
    for li = 1:3
        S = load(sprintf('table1_2_L%d.mat', Ls(li)));
        V = S.V_All; A = S.Var_AS_All; B = S.Var_BC_All;
        fprintf('L=%d\n', Ls(li));
        for cp = 1:5
            row = sprintf('%d', CPd(cp));
            for ci = 1:2
                co = coefIdx(ci);
                vb = squeeze(V(cp,co,est,:)); ab = squeeze(A(cp,co,est,:)); bb = squeeze(B(cp,co,est,:));
                ob = isfinite(vb) & isfinite(ab) & ab>0 & isfinite(bb) & bb>0;
                sA = mean(abs((vb(ob)-truth(ci))./sqrt(ab(ob)))>z);
                sB = mean(abs((vb(ob)-truth(ci))./sqrt(bb(ob)))>z);
                cA = mean(abs(vb(ob)-truth(ci))<z*sqrt(ab(ob)));
                cB = mean(abs(vb(ob)-truth(ci))<z*sqrt(bb(ob)));
                len = mean(sqrt(bb(ob)))/mean(sqrt(ab(ob)));
                if Ls(li)==0 && est==1
                    lenstr = 'n.a.';
                else
                    lenstr = sprintf('%.3f', len);
                end
                row = [row sprintf(' & %.3f & %.3f & %.3f & %.3f & %s', sA, sB, cA, cB, lenstr)];
            end
            fprintf('%s \\\n', row);
        end
    end
end
disp('EXTRACT_TABLE1_TABLE2_DONE');

function s = fmt(x)
    if x > 999; s = sprintf('%.2E', x); else; s = sprintf('%.3f', x); end
end
