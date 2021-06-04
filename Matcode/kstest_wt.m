function [testrst,prob,stats] = kstest_wt(x, wt, nullOBJ, alpha)
% 
% This function is essentially trying to do the same thing as
% as the built-in function 'kstest' does, which is One-sample
% Kolmogorov-Smirnov test. Tests if a sample comes from a 
% continuous distribution with specified parameters, against
% the alternative that it does not come from that distribution.
%
% The difference is that, this function deals with the data that
% has unequal weights, i.e. not 1s. The implementation follows
% procedure in Numerical Recipes Ch. 14.3.3, so that it might
% not be the same as matlab. 
%
% Given an array data[0..n-1], and given a user-supplied function
% of a single variable func that is a cumulative distribution 
% function ranging from 0 (for smallest values of its argument) 
% to 1 (for largest values of its argument), this routine returns
% the KS statistic d and the p-value prob.
%
%   Small values of prob show that the cumulative distribution 
%   function of data is significantly different from func.
%
% The array data is modified by being sorted into ascending order.
%
%   OUTPUT:
%       testrst:    test decision, 1 reject, 0 failure to reject
%       prob:       p-value probability p of the hypothesis testm
%       stats:      the structure stats containing information like
%                   chi-square statistic 
%   
%   INPUT:
%       x:       binned data set 1
%       wt:       expected data set whose distribution is known
%       refobj:       bin width, for additional check for the number 
%                   of constraints
%       alpha:      significance level, to be compared with prob
%       
% By Chao Song, chaosong@princeton.edu
% First created date:   2020/10/08
% Last modified date:   2020/10/08

% reshape to one column
x = reshape(x, [],1);
wt = reshape(wt, [],1);
data = [x wt];

% sort it 
data = sortrows(data,1);

% % empirical cumulative density function, CDF
empcdf = cumsum(data(:,2)/sum(data(:,2)));
xloc = data(:,1);

% or use 'ecdf', basically when number of samples are large, this two not distinguishable
% [empcdf,xloc] = ecdf(data(:,1),'Frequency',data(:,2));

% the CDF from a known reference distribution
% xi = min(xloc(1)):0.01:max(xloc(end));
% nullcdfi = cdf(nullOBJ, xi);
% nullcdf = interp1(xi,nullcdfi,xloc,'linear');
nullcdf = cdf(nullOBJ, xloc);

% K-S statistic
ksD = max(abs(empcdf - nullcdf));

% effective number of data points
en = sum(wt);    % sum of the weights;
ensqrt = sqrt(en);

% the probability is the Kolmogorov-Smirnov complementary cumulative distribution function (cdf)
% Q_KS(z). Call function 'kscdf_q' written by Chao Song
z = (ensqrt+0.12+0.11/ensqrt)*ksD;
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


stats = struct('effectN',en,'ksstat',ksD,'empcdf',empcdf);









