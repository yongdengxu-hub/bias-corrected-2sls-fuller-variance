function make_tables_nn(root)
% Emit CI-coverage LaTeX bodies for the non-normal-error robustness runs
% (run_CI_nn): t_5 (dist 3) and chi^2_4 (dist 4), rho=0.279, L=0,2,4.
% Same 9-column format as the S.4 CI tables; n.a. for L=0 2SLS length (no finite variance).
    if nargin<1||isempty(root), root='.'; end
    z=1.959964; cp=[8 10 15 25 70]; est={'2SLS','Fuller'};
    dists=[3 4]; dname={'t5','chi2'};
    for di=1:2
        d=dists(di);
        fprintf('\n@@@ CI_%s @@@\n', dname{di});
        for L=[0 2 4]
            fprintf('\\multicolumn{9}{c}{$L=%d$} \\\\\n', L);
            for p=1:2
                for c=1:5
                    f=sprintf('%s/raw_nn/dist%d_L%d/cell_CP%02d.mat',root,d,L,cp(c));
                    if ~isfile(f); continue; end
                    S=load(f); rb=cs(S,1,0.2,p,z); rg=cs(S,3,0.6,p,z); na=(L==0&&p==1);
                    fprintf('CP=%d & %s & %.3f & %.3f & %s && %.3f & %.3f & %s \\\\\n',...
                        cp(c),est{p},rb(1),rb(2),tn(na,'n.a.',sprintf('%.3f',rb(3))),...
                        rg(1),rg(2),tn(na,'n.a.',sprintf('%.3f',rg(3))));
                end
            end
        end
    end
    fprintf('\nDONE\n');
end
function s=tn(c,a,b), if c, s=a; else, s=b; end, end
function r=cs(S,col,tv,p,z)
    b=squeeze(S.V_All(1,col,p,:)); va=squeeze(S.Var_AS_All(1,col,p,:)); vb=squeeze(S.Var_BC_All(1,col,p,:));
    ok=isfinite(b)&isfinite(va)&va>0; b=b(ok);va=va(ok);vb=vb(ok);
    g=~isfinite(vb)|vb<=0|vb>3*va; vb(g)=va(g);
    covAS=mean(abs(b-tv)<z*sqrt(va)); covBC=mean(abs(b-tv)<z*sqrt(vb));
    lenr=mean(sqrt(vb))/mean(sqrt(va)); bias=mean(b)-tv; vr=var(b);
    r=[covAS covBC lenr bias vr mean(va)/vr mean(vb)/vr];
end
