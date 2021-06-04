function estCI = confidence_interval_general(est,se,nu,conf)
% 
% This function is to compute the upper and lower confidence interval of mean of
% a sample set.
%
% INPUT:
%   est:    estimated value
%   se:     standard error of the estimate
%   nu:     degree of freedom of the t-distribution 
%   conf:   percentage of confidence you want, e.g. 95
%
% OUTPUT:
%   muCI:   the desired confidence interval
%
%
% By Chao Song, chaosong@princeton.edu
% First created date:   2020/11/25
% Last modified date:   2020/11/25

pval = (100-conf)/100/2;    % desired prob. value 
crit = tinv([pval, 1-pval], nu);    % critical value of x based on pval on t pdf
CI = se .* crit;    % the confidence interval
estCI = CI + est;

% keyboard