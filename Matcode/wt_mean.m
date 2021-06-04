function mean_wt = wt_mean(x, wt)
% mean_wt = wt_mean(x, wt)
% This function is to calculate the weighted sample mean of the data
%
%
%
%   OUTPUT:
%       mean_wt:    weighted sample mean
%   
%   INPUT:
%       x:        data set 
%       wt:      corresponding weight
%       
%
%
% By Chao Song, chaosong@princeton.edu
% First created date:   2020/10/19
% Last modified date:   2020/10/19

mean_wt = sum(wt .* x, 1) ./ sum(wt, 1);