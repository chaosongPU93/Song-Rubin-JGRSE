function [testrst,prob,stats] = ttest2_wt(x1, wt1, x2, wt2, alpha, WtTypeFlag, EqualVarFlag)
% 
% This function is essentially trying to do the same thing as
% as the built-in function 'ttest2' does, which is Two-sample
% student-t test. Tests if two sample sets have the same mean,
% while they might or might NOT have an equal variance 
%
% The difference is that, this function deals with the data that
% has unequal weights, i.e. not 1s. The implementation follows
% procedure in Numerical Recipes Ch. 14.2.1, so that it might
% not be the same as matlab. 
%
% Given the arrays data1[0..n1-1] and data2[0..n2-1], this 
% routine returns Student??s t as t, and its p-value as prob.
%
%   Small values of prob indicating that the arrays have 
%   significantly different means.
%
% The data arrays are allowed to be drawn from populations
% with unequal variances.
%
%
%
%   OUTPUT:
%       testrst:    test decision, 1 reject, 0 failure to reject
%       prob:       p-value probability p of the hypothesis testm
%       stats:      the structure stats containing information like
%                   t statistic 
%   
%   INPUT:
%       x1:       data set 1
%       wt1:      weight of set 1
%       x2:       data set 2
%       wt2:      weight of set 2
%       alpha:      significance level, to be compared with prob
%       WtTypeFlag:     flag to indicate which weighting scheme is 
%                       used to compute unbiased weighted variance
%                       and degree of freedom, reliability weights (1)
%                       frequency weights (2)
%       EqualVarFlag:   flag to indicate if assuming x1 and x2 have 
%                       equal (1) or unequal (0) variances     
%       
%
%
% By Chao Song, chaosong@princeton.edu
% First created date:   2020/10/12
% Last modified date:   2020/10/29

% reshape to one column
x1 = reshape(x1, [],1);
wt1 = reshape(wt1, [],1);
x2 = reshape(x2, [],1);
wt2 = reshape(wt2, [],1);

% compute the weighted mean
mean_wt1 = wt_mean(x1, wt1);
mean_wt2 = wt_mean(x2, wt2);

% compute the weighted variance
%%% NOTE:
%%% there is in fact a lot of debate here on 'reliability weights', check wikipedia article:
%%% https://en.wikipedia.org/wiki/Weighted_arithmetic_mean#Variance_weights, AND this:
%%% https://stats.stackexchange.com/questions/61225/correct-equation-for-weighted-unbiased-sample-covariance/61298#61298
%%% https://stats.stackexchange.com/questions/51442/weighted-variance-one-more-time
%%% on this topic
sqsum1 = sum(wt1 .* (x1 - mean_wt1).^2);
sqsum2 = sum(wt2 .* (x2 - mean_wt2).^2);

if EqualVarFlag == 1    % assume x1 and x2 have unknown but equal variance
    if WtTypeFlag == 1
        %%% From the article, it seems that our case is more close to 'reliability
        %%% weights', thus the unbiased weighted variance is as follows:
        n1 = sum(wt1);
        v1 = sum(wt1.^2);
        n2 = sum(wt2);
        v2 = sum(wt2.^2);
        n1_eff = n1^2 / v1;
        n2_eff = n2^2 / v2;
        var_wt1 = sqsum1 / (n1 - (v1/n1));
        var_wt2 = sqsum2 / (n2 - (v2/n2));
        
    elseif WtTypeFlag == 2
        %%% otherwise for 'frequency weights', e.g. observe value A for 5 times, then weight of 
        %%% A is 5, the corresponding unbiased weighted variance would be:
        n1 = sum(wt1);
        n2 = sum(wt2);
        n1_eff = n1;
        n2_eff = n2;
        var_wt1 = sqsum1 / (n1_eff - 1);
        var_wt2 = sqsum2 / (n2_eff - 1);
        
    else
        error('WRONG WtTypeFlag value, assign 1 or 2!');
        %%% The implementation in Matlab function 'var' if you add weight option is in fact 
        %%% computing the biased weighted variance:
        % var_wt1 = sqsum1 / sum(wt1);
        % var_wt2 = sqsum2 / sum(wt2);
    end

    df = n1_eff + n2_eff -2;       % degree of freedom

    %%% the following formula comes from NR
    sqsum = sqsum1 + sqsum2;
    SE = sqrt(sqsum / df * (1/n1_eff + 1/n2_eff));
    meandiff = mean_wt1 - mean_wt2;
    tstat = meandiff / SE;
    
    % in NR, the prob is incomplete beta function A(t|v) = 1-I_x(1/2*v, 1/2), where v is df, x is
    % v/(v+((t-mu)/sigma)^2)
    % in matlab, incomplete beta function is implemented as 'betainc(x,z,w)' for the formula I_x(z,w)
    prob = betainc(df/(df+tstat^2), 0.5*df, 0.5, 'lower');
    
%     % Alternatively, matlab uses 'tcdf' for computing prob on student's t distribution CDF
%     prob = 2 * tcdf(-abs(tstat),df);

    % compute the confidence interval based on alpha as well
    spread = tinv(1 - alpha / 2, df) * SE;
    ci = [ meandiff - spread, meandiff + spread ];
    
    
elseif EqualVarFlag == 0     % assume x1 and x2 have unknown and unequal variance   
    
    if WtTypeFlag == 1
        %%% From the article, it seems that our case is more close to 'reliability
        %%% weights', thus the unbiased weighted variance is as follows:
        n1 = sum(wt1);
        v1 = sum(wt1.^2);
        n2 = sum(wt2);
        v2 = sum(wt2.^2);
        n1_eff = n1^2 / v1;
        n2_eff = n2^2 / v2;
        var_wt1 = sqsum1 / (n1 - (v1/n1));
        var_wt2 = sqsum2 / (n2 - (v2/n2));
        
    elseif WtTypeFlag == 2
        %%% otherwise for 'frequency weights', e.g. observe value A for 5 times, then weight of 
        %%% A is 5, the corresponding unbiased weighted variance would be:
        n1 = sum(wt1);
        n2 = sum(wt2);
        n1_eff = n1;
        n2_eff = n2;
        var_wt1 = sqsum1 / (n1_eff - 1);
        var_wt2 = sqsum2 / (n2_eff - 1);
        
    else
        error('WRONG WtTypeFlag value, assign 1 or 2!');
        %%% The implementation in Matlab function 'var' if you add weight option is in fact 
        %%% computing the biased weighted variance:
        % var_wt1 = sqsum1 / sum(wt1);
        % var_wt2 = sqsum2 / sum(wt2);
    end

    %%% the following formula comes from NR
    df = (var_wt1/n1_eff + var_wt2/n2_eff)^2 / ...
         ((var_wt1/n1_eff)^2/(n1_eff-1) + (var_wt2/n2_eff)^2/(n2_eff-1));      % degree of freedom
    SE = sqrt(var_wt1/n1_eff + var_wt2/n2_eff);
    meandiff = mean_wt1 - mean_wt2;
    tstat = meandiff / SE;
    
    % in NR, the prob p-value is incomplete beta function I_x(1/2*v, 1/2), where v is df, x is
    % v/(v+((t-mu)/sigma)^2)
    % in matlab, incomplete beta function is implemented as 'betainc(x,z,w)' for the formula I_x(z,w)
    prob = betainc(df/(df+tstat^2), 0.5*df, 0.5, 'lower');
    
%     % Alternatively, matlab uses 'tcdf' for computing prob on student's t distribution CDF
%     prob = 2 * tcdf(-abs(tstat),df);

    % compute the confidence interval based on alpha as well
    spread = tinv(1 - alpha / 2, df) * SE;
    ci = [ meandiff - spread, meandiff + spread ];

end

% 'prob' is the probability of observing a test statistic as extreme as, or more extreme than, the
% observed value under the null hypothesis. Small values of p cast doubt on the validity of the null
% hypothesis.
% the test decision is made by comparing 'prob' with significance level 'alpha'
if prob <= alpha
    testrst = 1;    % reject the null hypothesis
else
    testrst = 0;    % do not reject the null hypothesis
end
    
stats.t = tstat;
stats.df = df;
stats.mean_wt1 = mean_wt1;
stats.var_wt1 = var_wt1;
stats.mean_wt2 = mean_wt2;
stats.var_wt2 = var_wt2;
stats.n1_eff = n1_eff;
stats.n2_eff = n2_eff;
stats.meandiff = meandiff;
stats.ci = ci;


















