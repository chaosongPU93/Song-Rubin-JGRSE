function [testrst,prob,stats] = kstest2_wt(x1, wt1, x2, wt2, alpha)
% 
% This function is essentially trying to do the same thing as
% as the built-in function 'kstest2' does, which is Two-sample
% Kolmogorov-Smirnov test. Tests if two samples come from the
% same continuous distribution, against the alternative that 
% they do not come from the same distribution.
%
% The difference is that, this function deals with the data that
% has unequal weights, i.e. not 1s. The implementation follows
% procedure in Numerical Recipes Ch. 14.3.3, so that it might
% not be the same as matlab. 
%
% Given an array data1[0..n1-1], and an array data2[0..n2-1], 
% this routine returns the K??S statistic d and the p-value 
% prob for the null hypothesis that the data sets are drawn
% from the same distribution.
%
%   Small values of prob show that the cumulative distribution
%   function of data1 is significantly different from that of data2.
%
% The arrays data1 and data2 are modified by being sorted into ascending order.
%
%   OUTPUT:
%       testrst:    test decision, 1 reject, 0 failure to reject
%       prob:       p-value probability p of the hypothesis testm
%       stats:      the structure stats containing information like
%                   chi-square statistic 
%   
%   INPUT:
%       x1:       binned data set 1
%       wt1:      expected data set whose distribution is known
%       x2:       bin width, for additional check for the number 
%                 of constraints
%       wt2:    
%       alpha:      significance level, to be compared with prob
%       
%
%
% By Chao Song, chaosong@princeton.edu
% First created date:   2020/10/07
% Last modified date:   2020/10/07

% reshape to one column
x1 = reshape(x1, [],1);
wt1 = reshape(wt1, [],1);
x2 = reshape(x2, [],1);
wt2 = reshape(wt2, [],1);

% define bin edges
binEdges    =  [-inf ; sort([x1;x2]) ; inf];

% compute the weights (counts) in each bin
binCounts1 = zeros(length(binEdges)-1, 1);
binCounts2 = zeros(length(binEdges)-1, 1);
for i = 1: length(binEdges)-1
    wt1obj = wt1(x1>=binEdges(i) & x1<binEdges(i+1));
    binCounts1(i) = sum(wt1obj);
    wt2obj = wt2(x2>=binEdges(i) & x2<binEdges(i+1));
    binCounts2(i) = sum(wt2obj);
end

empcdf1  =  cumsum(binCounts1)./sum(binCounts1);
empcdf2  =  cumsum(binCounts2)./sum(binCounts2);

% K-S statistic
ks2D = max(abs(empcdf1 - empcdf2));

% effective number of data points
n1 = sum(wt1);    % sum of the weights;
n2 = sum(wt2);
en = n1 * n2 /(n1 + n2);
ensqrt = sqrt(en);

% the probability is the Kolmogorov-Smirnov complementary cumulative distribution function (cdf)
% Q_KS(z). Call function 'kscdf_q' written by Chao Song
z = (ensqrt+0.12+0.11/ensqrt)*ks2D;
prob = kscdf_q(z);

% 'prob' is the probability of observing a test statistic as extreme as, or more extreme than, the
% observed value under the null hypothesis. Small values of p cast doubt on the validity of the null
% hypothesis.
% the test decision is made by comparing 'prob' with significance level 'alpha'
if prob <= alpha
    testrst = 1;    % reject the null hypothesis
else
    testrst = 0;    % do not reject the null hypothesis
end
    

stats = struct('effectN',en,'ks2stat',ks2D,'empcdf1',empcdf1,'empcdf2',empcdf2,'binEdges',binEdges);











