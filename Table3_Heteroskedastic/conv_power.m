function [PowerT,ConvP,t] = conv_power(V1,VAR,B1,T)          % V1: estimated coefficients; VAR: estimated variance; B1: true coefficients;
   
    [m,n] = size(V1);
    CV = tinv(0.975,T-m);                   % m: number of regressors in the first equation,  95% confidence intervals  
    CV_uu = tinv(0.975,T-m);                   % m: number of regressors in the first equation,  95% confidence intervals  
    CV_ll = tinv(0.025,T-m);
    count1= NaN(m,n);
    count2= NaN(m,n);
    t= NaN(m,n);
    if isempty(VAR)
        VAR=var(squeeze(V1)');
        VAR=VAR'*ones(1,n);
    end  
      for rptr = 1:n
          conf_intv_u = V1(:,rptr) + CV*sqrt(VAR(:,rptr));
          conf_intv_l = V1(:,rptr) - CV*sqrt(VAR(:,rptr));
          count1(:,rptr)= B1' >= conf_intv_l & B1' <=conf_intv_u;
            
          t(:,rptr) = (V1(:,rptr)-B1')./ sqrt(VAR(:,rptr));
          count2(:,rptr) =t(:,rptr) > CV_uu;
          count3(:,rptr) =t(:,rptr) < CV_ll;
      end
      ConvP = mean(count1,2);
      PowerT = mean(count3,2);   % left or right tail test
      %CV_l=quantile(t',0.05);
      %CV_u=quantile(t',0.95);
      PowerT = mean(t < CV_ll'*ones(1,n),2)+mean(t > CV_uu'*ones(1,n),2);       
end