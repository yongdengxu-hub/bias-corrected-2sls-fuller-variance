%% AJR Table 7: recompute BC variances with PER-COORDINATE REVERT (matching the fixed
%  homoskedastic Monte Carlo mechanism in doMonteCarlosSimulation_YX.m), replacing the
%  gated form used in the previous round. AJR is a homoskedastic application, so it should
%  follow the same convention as Table 1/2 (revert-to-uncorrected if the correction would be
%  negative), not the gate (which is specific to the Fuller heteroskedastic Monte Carlo).
clear; clc;
D = load_named_columns('Table7_Data.xlsx');
geo = {'logem4','lt100km','meantemp'};
specs = {
  'C1', {'avexpr'},          {'logem4'}, {'malfal94'},            @(D) isnan(D.malfal94)
  'C2', {'avexpr'},          {'logem4'}, {'malfal94','lat_abst'}, @(D) any(isnan([D.malfal94 D.lat_abst]),2)
  'C3', {'avexpr'},          {'logem4'}, {'leb95'},               @(D) isnan(D.leb95)
  'C4', {'avexpr'},          {'logem4'}, {'leb95','lat_abst'},    @(D) any(isnan([D.leb95 D.lat_abst]),2)
  'C5', {'avexpr'},          {'logem4'}, {'imr95'},               @(D) isnan(D.imr95)
  'C6', {'avexpr'},          {'logem4'}, {'imr95','lat_abst'},    @(D) any(isnan([D.imr95 D.lat_abst]),2)
  'C7', {'avexpr','malfal94'}, geo, {}, @(D) any(isnan([D.logem4 D.lt100km D.meantemp]),2)
  'C8', {'avexpr','leb95'},    geo, {}, @(D) any(isnan([D.logem4 D.lt100km D.meantemp D.avexpr D.leb95]),2)
  'C9', {'avexpr','imr95'},    geo, {}, @(D) any(isnan([D.logem4 D.lt100km D.meantemp D.avexpr D.imr95]),2)
};
n = size(specs,1);
R = struct();
for i = 1:n
    yv = D.logpgp95;
    Y2 = assemble(D, specs{i,2});
    X2 = assemble(D, specs{i,3});
    ctrl = assemble(D, specs{i,4});
    if isempty(ctrl), ctrl = zeros(numel(yv),0); end
    X1 = [ctrl ones(numel(yv),1)];
    X  = [X2 X1];
    drop = specs{i,5}(D);
    yv(drop)=[]; Y2(drop,:)=[]; X2(drop,:)=[]; ctrl(drop,:)=[]; X1(drop,:)=[]; X(drop,:)=[];
    [T,K] = size(X); k = size(X1,2); g = size(Y2,2); L = K-k-g; G = g+1; kz = size(X2,2);
    Feff = gweakivtest1(yv, Y2, [], X2);
    R.Feff(i) = Feff; R.L(i) = L; R.T(i) = T;
    for est = 1:2
        if est==1, [b, V, bias] = tsls_parts(yv, Y2, X, X1, k, T, g, L, G, K);
        else,      [b, V, bias] = fuller_parts(yv, Y2, X, X1, k, T, g, L, G, K); end
        Vr = V - bias;
        nonpos = (Vr <= 0);
        Vr(nonpos) = V(nonpos);
        for cc = 1:2
            R.b(i,est,cc)=b(cc); R.se(i,est,cc)=sqrt(V(cc)); R.seBC(i,est,cc)=sqrt(Vr(cc));
            R.t(i,est,cc)=b(cc)/sqrt(V(cc)); R.tBC(i,est,cc)=b(cc)/sqrt(Vr(cc));
            R.reverted(i,est,cc) = nonpos(cc);
        end
    end
end
fid = fopen('AJR_reverttable_out.txt','w');
blocks = {[1 1],'2SLS beta'; [2 1],'Fuller beta'; [1 2],'2SLS gamma'; [2 2],'Fuller gamma'};
for h = [1 fid]
  for bl = 1:4
    est=blocks{bl,1}(1); cc=blocks{bl,1}(2);
    fprintf(h,'\n%s\n', blocks{bl,2});
    fprintf(h,'coef '); fprintf(h,'& %.3f ', squeeze(R.b(:,est,cc))); fprintf(h,'\\\n');
    fprintf(h,'se   '); fprintf(h,'& (%.3f) ', squeeze(R.se(:,est,cc))); fprintf(h,'\\\n');
    fprintf(h,'seBC '); fprintf(h,'& {%.3f} ', squeeze(R.seBC(:,est,cc))); fprintf(h,'\\\n');
    fprintf(h,'t    '); fprintf(h,'& [%.3f] ', squeeze(R.t(:,est,cc))); fprintf(h,'\\\n');
    fprintf(h,'tBC  '); fprintf(h,'& <%.3f> ', squeeze(R.tBC(:,est,cc))); fprintf(h,'\\\n');
    fprintf(h,'revert '); fprintf(h,'& %d ', squeeze(R.reverted(:,est,cc))); fprintf(h,'\\\n');
  end
  fprintf(h,'\nF_eff '); fprintf(h,'& %.3f ', R.Feff); fprintf(h,'\\\n');
  fprintf(h,'(L,T) '); for i=1:n, fprintf(h,'& (%d,%d) ', R.L(i), R.T(i)); end; fprintf(h,'\\\n');
end
fclose(fid);
fprintf('\nSaved AJR_reverttable_out.txt\n');

function M = assemble(D, names)
    M = [];
    for c = 1:numel(names), M = [M D.(names{c})]; end %#ok<AGROW>
end
function [b, Var_AS, bias] = tsls_parts(y, Y2, X, X1, k, T, g, L, G, K)
    [b, ~] = kClass(y, Y2, X, X1, 1, T);
    Yb = X*(X\Y2);
    Q = [Yb'*Yb Yb'*X1; X1'*Yb X1'*X1];
    Z1 = [Y2 X1];
    e = y - Z1*b;
    I = eye(G-1+k);
    S11 = e'*e/(T-g-k);
    Var_AS = diag(S11*(I/Q));
    V2 = Y2 - Yb;
    v1 = y - X*(X\y);
    Tau = V2'*(v1-V2*b(1:g))/(T-K);
    C = zeros(G-1+k); C(1:G-1,1:G-1) = V2'*V2/(T-K);
    C1 = zeros(G-1+k); C1(1:G-1,1:G-1) = (Tau*Tau')/S11;
    bias = diag(S11*((L+1)*(Q\C1/Q) + trace(Q\C)*(I/Q)));
end
function [b, Var_AS, bias] = fuller_parts(y, Y2, X, X1, k, T, g, L, G, K)
    Yd = [y Y2];
    YdtYd = Yd'*Yd;
    Wsdd = YdtYd - (Yd'*X1*(X1\Yd));
    Wdd  = YdtYd - (Yd'*X*(X\Yd));
    lambda = min(eig(Wdd\Wsdd)) - 1/(T-K);
    [b, ~] = kClass(y, Y2, X, X1, lambda, T);
    V2 = Y2 - X*(X\Y2);
    UL = Y2'*Y2 - lambda*(V2'*V2);
    Q1 = [UL Y2'*X1; X1'*Y2 X1'*X1];
    Z1 = [Y2 X1];
    e = y - Z1*b;
    S11 = e'*e/(T-g-k);
    I = eye(size(Q1,1));
    Var_AS = diag(S11*(I/Q1));
    Yb = X*(X\Y2);
    Q = [Yb'*Yb Yb'*X1; X1'*Yb X1'*X1];
    v1 = y - X*(X\y);
    Tau = V2'*(v1-V2*b(1:g))/T;
    C = zeros(G-1+k); C(1:G-1,1:G-1) = V2'*V2/T;
    C1 = zeros(G-1+k); C1(1:G-1,1:G-1) = (Tau*Tau')/S11;
    C2 = C - C1;
    bias = diag(S11*((Q\C/Q) + 5*(Q\C1/Q) + 2*trace(Q\C)*(I/Q) - trace(Q\C1)*(I/Q) - trace(Q\C2)*(I/Q) - L*(Q\C2/Q)));
end
