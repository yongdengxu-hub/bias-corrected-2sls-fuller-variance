% Build two table bodies for the heteroskedastic high-endogeneity (rho=0.8, delta~0.78) run:
%   het78rho08_var_rows.txt    : Bias, Var, HC0/Var, BC/Var   (tab:het_var08 format)
%   het78rho08_sizeci_rows.txt : size (AS,BC) and 95% CI coverage (AS,BC)  (tab:het_sizeci08 format)
% Adapted 2026-07-04 from extract_het_rho08.m (which read the delta~0.6 default design, not
% the delta~0.78 both table captions actually specify) -- also drops the stale "revert BC to
% HC0 if negative or >3x" cap, redundant now that Var_BC is always positive by construction
% (the gated two-branch estimator, doMonteCarlosSimulation_YX_Hetro.m).
clear; z = 1.959964;
Ls = [0 2 4]; CPd = [8 10 15 25 70]; tb = 0.2; tg = 0.6;
estlab = {'Panel A: 2SLS','Panel B: Fuller'};
fv = fopen('het78rho08_var_rows.txt','w');
fs = fopen('het78rho08_sizeci_rows.txt','w');
for est = 1:2
    fprintf(fv,'\\hline\n\\multicolumn{9}{c}{\\textbf{%s}} \\\\\n\\hline\n', estlab{est});
    fprintf(fs,'\\hline\n\\multicolumn{9}{c}{\\textbf{%s}} \\\\\n\\hline\n', estlab{est});
    for li = 1:3
        S = load(sprintf('sizeci78rho08_L%d.mat', Ls(li)));
        V = S.V_All; A = S.Var_AS_All; B = S.Var_BC_All;
        cm = mean(mean(mean(V,4),1),3); cm = cm(:);
        [~,bi] = min(abs(cm-tb)); [~,gi] = min(abs(cm-tg));
        fprintf(fv,'\\multicolumn{9}{c}{$L=%d$} \\\\\n', Ls(li));
        fprintf(fs,'\\multicolumn{9}{c}{$L=%d$} \\\\\n', Ls(li));
        for cp = 1:5
            vb=squeeze(V(cp,bi,est,:)); ab=squeeze(A(cp,bi,est,:)); bb=squeeze(B(cp,bi,est,:));
            vg=squeeze(V(cp,gi,est,:)); ag=squeeze(A(cp,gi,est,:)); bg=squeeze(B(cp,gi,est,:));
            ob=isfinite(vb)&isfinite(ab)&ab>0; og=isfinite(vg)&isfinite(ag)&ag>0;
            % --- variance / bias table ---
            Bb=mean(vb(ob))-tb; Vb=var(vb(ob)); rH=mean(ab(ob))/Vb; rB=mean(bb(ob))/Vb;
            Bg=mean(vg(og))-tg; Vg=var(vg(og)); rHg=mean(ag(og))/Vg; rBg=mean(bg(og))/Vg;
            fprintf(fv,'%d & %.3f & %s & %s & %s & %.3f & %s & %s & %s \\\\\n', CPd(cp), ...
                Bb, f3(Vb), fr(rH), fr(rB), Bg, f3(Vg), fr(rHg), fr(rBg));
            % --- size / CI table ---
            sAb=mean(abs((vb(ob)-tb)./sqrt(ab(ob)))>z); sBb=mean(abs((vb(ob)-tb)./sqrt(bb(ob)))>z);
            cAb=mean(abs(vb(ob)-tb)<z*sqrt(ab(ob)));    cBb=mean(abs(vb(ob)-tb)<z*sqrt(bb(ob)));
            sAg=mean(abs((vg(og)-tg)./sqrt(ag(og)))>z); sBg=mean(abs((vg(og)-tg)./sqrt(bg(og)))>z);
            cAg=mean(abs(vg(og)-tg)<z*sqrt(ag(og)));    cBg=mean(abs(vg(og)-tg)<z*sqrt(bg(og)));
            fprintf(fs,'%d & %.3f & %.3f & %.3f & %.3f & %.3f & %.3f & %.3f & %.3f \\\\\n', ...
                CPd(cp), sAb,sBb,cAb,cBb, sAg,sBg,cAg,cBg);
        end
    end
end
fclose(fv); fclose(fs);
disp('EXTRACT78RHO08_DONE');
function s=f3(x), if x>=1000, s=sprintf('%.2E',x); else, s=sprintf('%.3f',x); end, end
function s=fr(x), if ~isfinite(x), s='n.a.'; elseif x>=100, s=sprintf('%.2E',x); else, s=sprintf('%.3f',x); end, end
