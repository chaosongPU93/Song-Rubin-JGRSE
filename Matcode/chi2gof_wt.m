function [testrst,prob,stats] = chi2gof_wt(bin1, ebin, binw, nparam, alpha)
% 'CHI2GOF_WT' One-sample chi-square goodness-of-fitness test with weights
%
% This function is essentially trying to do the same thing as
% as the built-in function 'chi2gof' does, which is Chi-square
% goodness-of-fit test. Tests if a sample comes from a specified
% distribution, against the alternative that it does not come 
% from that distribution.
%
% The difference is that, this function deals with the data that
% has unequal weights, i.e. not 1s. The implementation follows
% procedure in Numerical Recipes Ch. 14.3.1, so that it might
% not be the same as matlab. 
%
% Given the array bins[0..nbins-1] containing the observed numbers
% of events, and an array ebins[0..nbins-1] containing the expected
% numbers of events, and given the number of constraints knstrn 
% (normally one), this routine returns (trivially) the number of 
% degrees of freedom df, and (nontrivially) the chi-square chsq and
% the p-value prob.
%
%   A small value of prob indicates a significant 
%   difference between the distributions bins and ebins.
%
% Note that bins and ebins are both double arrays, although bins 
% will normally contain integer values.
%
%   OUTPUT:
%       testrst:    test decision, 1 reject, 0 failure to reject
%       prob:       p-value probability p of the hypothesis testm
%       stats:      the structure stats containing information like
%                   chi-square statistic 
%   
%   INPUT:
%       bin1:       binned data set 1
%       ebin:       expected data set whose distribution is known
%       binw:       bin width, for additional check for the number 
%                   of constraints
%     * nparam:     number of parameters that are used if the expected
%                   counts 'ebin' depends on estimated parameters. If
%                   not estimated by fitting a distribution to data, it
%                   would be 0 otherwise.
%       alpha:      significance level, to be compared with prob
%       
% By Chao Song, chaosong@princeton.edu
% First created date:   2020/10/08
% Last modified date:   2020/10/08

% vectorization 
bin1 = reshape(bin1,[],1);
ebin = reshape(ebin,[],1);
nbins = length(bin1);   % number of bins

% the number of constraints
knstrn = 0 + nparam;     % default
% additional check for the number of constraints
if sum(bin1) == sum(ebin) && sum(bin1) == 1     % i.e. the value of each bin is probability
    knstrn = 1 + nparam;     % according to NR, the sum of counts for two sets are same 
elseif sum(bin1)*binw == sum(ebin)*binw && sum(bin1)*binw == 1  % i.e. the value of each bin is pdf
    knstrn = 1 + nparam;     % according to NR, the sum of counts for two sets are same 
elseif sum(bin1) == sum(ebin) && sum(bin1)~=1 && sum(ebin)~=1  % pseudo counts, i.e., sum of weights in that bin
    knstrn = 1 + nparam;     % according to NR, the sum of counts for two sets are same  
elseif sum(bin1) ~= sum(ebin) && sum(bin1)~=1 && sum(ebin)~=1  % pseudo counts, but not same, due to different number of detections
    knstrn = 0 + nparam;
end

%  degree of freedom    
df = nbins - knstrn;    %  degree of freedom

% chi-square statistic formulation from NR
chisq = 0;
for i = 1: nbins
    if ebin(i) < 0 || (ebin(i) == 0 && bin1(i) > 0)
        warning('Bag expected number in some bins');
    end
    if bin1(i) == 0 && ebin(i) == 0     % no data means one less degree of freedom
        df = df-1;
    else
        tmp = bin1(i) - ebin(i);
        chisq = chisq + tmp.^2 ./ ebin(i); 
    
    end
end

% the probability is the complementary incomplete gamma function Q(a,x), where a = df/2, x =  chi2/2
% in matlab, Y = gammainc(X,A,'upper') returns the upper incomplete gamma function Q(x,a)
prob = gammainc(chisq/2, df/2, 'upper'); 

%%% alternatively, you could use 'chi2cdf' to calculate the probability which gives same result
% prob = chi2cdf(chisq, df, 'upper');
%%% Also, another matlab function 'chi2pval' uses EXACTLY the same syntax as we use here.
% prob = chi2pval(chisq, df);

% 'prob' is the probability of observing a test statistic as extreme as, or more extreme than, the
% observed value under the null hypothesis. Small values of p cast doubt on the validity of the null
% hypothesis.
% the test decision is made by comparing 'prob' with significance level 'alpha'
if prob <= alpha
    testrst = 1;    % reject the null hypothesis
else
    testrst = 0;    % do not reject the null hypothesis
end
    
stats = struct('chi2stat',chisq,'knstrn',knstrn,'df',df);









