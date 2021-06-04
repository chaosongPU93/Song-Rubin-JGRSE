function [var_wt] = wt_var(x, wt, WtTypeFlag)
% [var_wt] = wt_var(x, wt, WtTypeFlag)
% This function is to calculate the weighted sample variance of
% data with weights
%
%
%
%   OUTPUT:
%       var_wt:    
%   
%   INPUT:
%       x:        data set 
%       wt1:      weights
%       WtTypeFlag:     flag to indicate which weighting scheme is 
%                       used to compute unbiased weighted variance
%                       and degree of freedom, reliability weights (1)
%                       frequency weights (2), OR biased type (0)     
%       
%
%
% By Chao Song, chaosong@princeton.edu
% First created date:   2020/10/19
% Last modified date:   2020/10/19

mean_wt = wt_mean(x, wt);
sqsum = sum(wt .* (x - mean_wt).^2, 1);

% compute the weighted variance
%%% NOTE:
%%% there is in fact a lot of debate here on 'reliability weights', check wikipedia article:
%%% https://en.wikipedia.org/wiki/Weighted_arithmetic_mean#Variance_weights, AND this:
%%% https://stats.stackexchange.com/questions/61225/correct-equation-for-weighted-unbiased-sample-covariance/61298#61298
%%% https://stats.stackexchange.com/questions/51442/weighted-variance-one-more-time
%%% on this topic
if WtTypeFlag == 1
    %%% From the article, it seems that our case is more close to 'reliability
    %%% weights', thus the unbiased weighted variance is as follows:
    n1 = sum(wt, 1);
    v1 = sum(wt.^2, 1);
    var_wt = sqsum ./ (n1 - (v1./n1));
    
elseif WtTypeFlag == 2
    %%% otherwise for 'frequency weights', e.g. observe value A for 5 times, then weight of
    %%% A is 5, the corresponding unbiased weighted variance would be:
    n1 = sum(wt, 1);
    n1_eff = n1;
    var_wt = sqsum ./ (n1_eff - 1);
    
elseif WtTypeFlag == 0
    %%% The implementation in Matlab function 'var' if you add weight option is in fact
    %%% computing the biased weighted variance:
    n1 = sum(wt, 1);
    n1_eff = n1;
    var_wt = sqsum ./ n1_eff;
else
    error('WRONG WtTypeFlag value, assign 0, 1 or 2!');

end
