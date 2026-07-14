function [Y LY Vb] = genY_Hetro(y0, X, rho, Gamma, Pi, T, G, dist)
    % Heteroskedastic DGP generator with selectable conditional-variance form.
    % HET_TYPE (global) selects the form; HET_Z (global) is a fixed auxiliary
    % vector (~N(0,1), drawn once, independent of X) used by types 1 and 6.
    % Default (HET_TYPE unset) = type 2 = the original design-aligned form,
    % so existing drivers (main_HC_v8) reproduce Table 3 unchanged.
    global HET_TYPE HET_Z HET_B
    switch dist
        case 1; e = randn(T, G);
        case 2; e = -sqrt(12)/2 + sqrt(12).*rand(T, G);
        otherwise;
    end

    Vb = e*rho';
    x1 = X(:,2);                                  % included exogenous regressor
    if isempty(HET_TYPE), ht = 2; else, ht = HET_TYPE; end

    switch ht
        case 0                                    % homoskedastic
            sdmult = ones(T,1);
        case 1                                    % idiosyncratic: sigma^2 indep of X (balance holds)
            if isempty(HET_B), bcoef = 0.7; else, bcoef = HET_B; end
            v = exp(bcoef*HET_Z); v = v/mean(v);  sdmult = sqrt(v);
        case 2                                    % design-aligned baseline (= original Table 3)
            scale = sqrt(abs(x1));
            scale = scale / sqrt(mean(scale.^2));
            sdmult = sqrt(scale);
        case 3                                    % exponential / Harvey on x1
            x1s = (x1-mean(x1))/std(x1);
            v = exp(0.7*x1s);     v = v/mean(v);  sdmult = sqrt(v);
        case 4                                    % quadratic on x1
            x1s = (x1-mean(x1))/std(x1);
            v = 0.1 + x1s.^2;     v = v/mean(v);  sdmult = sqrt(v);
        case 5                                    % instrument / first-stage aligned (balance violated)
            fs  = X*Pi(:,2);                      % first-stage fit of endog regressor (~ sum of IVs)
            fss = (fs-mean(fs))/std(fs);
            v = 0.1 + fss.^2;     v = v/mean(v);  sdmult = sqrt(v);
        case 6                                    % two-regime (grouped) het
            v = ones(T,1); v(HET_Z>0) = 4; v = v/mean(v);  sdmult = sqrt(v);
        case 7                                    % design-aligned (= type 2) on the STRUCTURAL
            scale = sqrt(abs(x1));                % equation only; first stage left homoskedastic
            scale = scale / sqrt(mean(scale.^2)); % so that CP* (robust) = CP (homoskedastic/nominal)
            sdmult = sqrt(scale);
        otherwise
            sdmult = ones(T,1);
    end
    if ht == 7
        Vb(:,1) = Vb(:,1) .* sdmult;              % het only in the y1 (structural) reduced-form error
    else
        Vb = Vb .* sdmult;                        % het in all reduced-form errors
    end

    Y  = zeros(T, G);
    if ~any(any(Gamma))
        Y = X*Pi + Vb;
        LY = [];
    else
        XV = X*Pi + Vb;
        Y(1,:) = y0*Gamma + XV(1,:);
        for t = 2:T
            Y(t,:) = Y(t-1,:)*Gamma + XV(t,:);
        end
        LY = [y0; Y(1:T-1,:)];
    end
end
