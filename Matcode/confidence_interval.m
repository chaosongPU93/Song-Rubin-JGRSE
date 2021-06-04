function muCI = confidence_interval(mu,sigma,neff,conf)
% 
% This function is to compute the upper and lower confidence interval of mean of
% a sample set.
%
% INPUT:
%   mu: sample mean
%   sigma:  sample standard deviation
%   neff:   effective sample size, for data with weights of 1, it is the size of
%           samples; otherwise it is the sum of weights
%   conf:   percentage of confidence you want, e.g. 95
%
% OUTPUT:
%   muCI:   the desired confidence interval
%
%
% By Chao Song, chaosong@princeton.edu
% First created date:   2020/10/21
% Last modified date:   2020/10/21

se = sigma / sqrt(neff);    % standard error
nu = neff - 1;  % degree of freedom
pval = (100-conf)/100/2;    % desired prob. value 
crit = tinv([pval, 1-pval], nu);    % critical value of x based on pval on t pdf
CI = se .* crit;    % the confidence interval
muCI = CI + mu;

% keyboard