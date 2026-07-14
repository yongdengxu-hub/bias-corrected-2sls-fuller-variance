%% Replication of Section 6.2 (Angrist & Krueger 1991) of the Unbiased-Variance paper
%  Three specifications following Staiger & Stock (1997) / Andrews & Armstrong
%  (2017), men born 1930-1939, 1980 Census (N = 329,509).
%    I  : 3 quarter-of-birth instruments,  controls = year-of-birth dummies
%    II : 30 instruments (QOB + QOBxYOB),  controls = year-of-birth dummies
%    III: 30 instruments,  controls = YOB dummies + age, age^2 + demographics
%  Reports the schooling coefficient beta (Educ) and the YOB-1930 control
%  coefficient gamma, by 2SLS and Fuller, with asymptotic and bias-corrected
%  standard errors and t-statistics, the first-stage effective F, and (L,T).

clear; clc;
rng('default');

D = readtable('AK1991_3039.xlsx');
y  = D.LWKLYWGE;  Y2 = D.EDUC;  n = height(D);

yobD = dummyvar(categorical(D.YOB)); yobD = yobD(:,1:end-1);   % 9 YOB dummies (YOB30..YOB38)
qobD = dummyvar(categorical(D.QOB)); qobD = qobD(:,1:end-1);   % 3 QOB dummies
yqob = zeros(n, size(qobD,2)*size(yobD,2));                    % 27 QOBxYOB interactions
for i = 1:n, yqob(i,:) = reshape(yobD(i,:)'*qobD(i,:),1,[]); end

% Spec III adds age and age^2 (age in quarters), which are nearly collinear
% with the year-of-birth dummies and so sharply weaken the first stage.
agec = [D.AGEQ D.AGEQ.^2];

specs = {
  'I',   qobD,           [yobD ones(n,1)]
  'II',  [qobD yqob],    [yobD ones(n,1)]
  'III', [qobD yqob],    [yobD agec ones(n,1)]
};

fid = fopen('replicate_AK1991_out.txt','w');
for h = [1 fid]
    fprintf(h,'Angrist & Krueger (1991) replication  -  paper Section 6.2\n');
    fprintf(h,'%s\n', repmat('=',1,64));
    fprintf(h,'%-10s %10s %10s %10s\n','','Spec I','Spec II','Spec III');
end
B = cell(3,1);
for s = 1:3
    X2 = specs{s,2};  X1 = specs{s,3};  X = [X2 X1];
    [T,K]=size(X); k=size(X1,2); g=1; G=2; L=K-k-g;
    [bT,~,~,VAS_T,VBC_T] = TSLS_variance_BC(y, Y2, X, X1, k, T, g, L, G, K);
    [bF,~,~,VAS_F,VBC_F] = fuller_variance_BC(y, Y2, X, X1, k, T, g, L, G, K);
    Feff = gweakivtest1(y, Y2, X1, X2);
    B{s} = struct('b2',bT(1:2),'seA2',sqrt(VAS_T(1:2)),'seB2',sqrt(VBC_T(1:2)), ...
                  'bf',bF(1:2),'seAf',sqrt(VAS_F(1:2)),'seBf',sqrt(VBC_F(1:2)), ...
                  'F',Feff,'L',L,'T',T);
end

for h = [1 fid]
    blk(h,'2SLS  beta (Educ)',   B,'b2','seA2','seB2',1);
    blk(h,'Fuller beta (Educ)',  B,'bf','seAf','seBf',1);
    blk(h,'2SLS  gamma (YOB30)',  B,'b2','seA2','seB2',2);
    blk(h,'Fuller gamma (YOB30)', B,'bf','seAf','seBf',2);
    fprintf(h,'\n%-10s %10.3f %10.3f %10.3f\n','F_hat', B{1}.F, B{2}.F, B{3}.F);
    fprintf(h,'%-10s %10s %10s %10s\n','(L,T)', ...
        sprintf('(%d,%d)',B{1}.L,B{1}.T), sprintf('(%d,%d)',B{2}.L,B{2}.T), sprintf('(%d,%d)',B{3}.L,B{3}.T));
end
fclose(fid);
fprintf('\nSaved replicate_AK1991_out.txt\n');

function blk(h, ttl, B, bf, af, bcf, row)
    fprintf(h,'\n%s\n', ttl);
    fprintf(h,'%-10s %10.3f %10.3f %10.3f\n','b',    B{1}.(bf)(row), B{2}.(bf)(row), B{3}.(bf)(row));
    fprintf(h,'%-10s %10.3f %10.3f %10.3f\n','se',   B{1}.(af)(row), B{2}.(af)(row), B{3}.(af)(row));
    fprintf(h,'%-10s %10.3f %10.3f %10.3f\n','se_BC',B{1}.(bcf)(row),B{2}.(bcf)(row),B{3}.(bcf)(row));
    fprintf(h,'%-10s %10.3f %10.3f %10.3f\n','t',    B{1}.(bf)(row)/B{1}.(af)(row), B{2}.(bf)(row)/B{2}.(af)(row), B{3}.(bf)(row)/B{3}.(af)(row));
    fprintf(h,'%-10s %10.3f %10.3f %10.3f\n','t_BC', B{1}.(bf)(row)/B{1}.(bcf)(row),B{2}.(bf)(row)/B{2}.(bcf)(row),B{3}.(bf)(row)/B{3}.(bcf)(row));
end
