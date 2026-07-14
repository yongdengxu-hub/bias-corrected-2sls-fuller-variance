% Build the heteroskedastic size / CI-coverage table body (delta~0.78, rho=0.279) from the
% sizeci78_L*.mat files produced by main_HC_sizeci_delta78.m (proper Fuller sandwich +
% two-branch BC, always positive by construction -- no reversion cap needed).
% Adapted 2026-07-04 from extract_het_sizeci.m (which read the delta~0.6 default design,
% not the delta~0.78 that tab:het_sizeci's caption actually specifies).
clear; z = 1.959964;
Ls = [0 2 4]; CPd = [8 10 15 25 70]; tb = 0.2; tg = 0.6;
estlab = {'Panel A: 2SLS','Panel B: Fuller'};
fid = fopen('het78_sizeci_rows.txt','w');
for est = 1:2
    fprintf(fid,'\\hline\n\\multicolumn{9}{c}{\\textbf{%s}} \\\\\n\\hline\n', estlab{est});
    for li = 1:3
        S = load(sprintf('sizeci78_L%d.mat', Ls(li)));
        V = S.V_All; A = S.Var_AS_All; B = S.Var_BC_All;   % (CP,coeff,est,rep)
        cm = mean(mean(mean(V,4),1),3); cm = cm(:);         % mean per coeff
        [~,bi] = min(abs(cm-tb)); [~,gi] = min(abs(cm-tg));
        fprintf(fid,'\\multicolumn{9}{c}{$L=%d$} \\\\\n', Ls(li));
        for cp = 1:5
            vb=squeeze(V(cp,bi,est,:)); ab=squeeze(A(cp,bi,est,:)); bb=squeeze(B(cp,bi,est,:));
            vg=squeeze(V(cp,gi,est,:)); ag=squeeze(A(cp,gi,est,:)); bg=squeeze(B(cp,gi,est,:));
            ob=isfinite(vb)&isfinite(ab)&ab>0;
            og=isfinite(vg)&isfinite(ag)&ag>0;
            sAb=mean(abs((vb(ob)-tb)./sqrt(ab(ob)))>z);  sBb=mean(abs((vb(ob)-tb)./sqrt(bb(ob)))>z);
            cAb=mean(abs(vb(ob)-tb)<z*sqrt(ab(ob)));     cBb=mean(abs(vb(ob)-tb)<z*sqrt(bb(ob)));
            sAg=mean(abs((vg(og)-tg)./sqrt(ag(og)))>z);  sBg=mean(abs((vg(og)-tg)./sqrt(bg(og)))>z);
            cAg=mean(abs(vg(og)-tg)<z*sqrt(ag(og)));     cBg=mean(abs(vg(og)-tg)<z*sqrt(bg(og)));
            fprintf(fid,'%d & %.3f & %.3f & %.3f & %.3f & %.3f & %.3f & %.3f & %.3f \\\\\n', ...
                CPd(cp), sAb,sBb,cAb,cBb, sAg,sBg,cAg,cBg);
        end
    end
end
fclose(fid);
disp('EXTRACT78_SIZECI_DONE'); type('het78_sizeci_rows.txt');
