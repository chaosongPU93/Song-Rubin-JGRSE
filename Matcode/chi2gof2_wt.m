function [testrst,prob,stats] = chi2gof2_wt(bin1, bin2, binw, alpha)
% 'CHI2GOF2_WT' Two-sample chi-square goodness-of-fitness test with weights
%
% This function is essentially trying to do the same thing as
% as the built-in function 'chi2gof' does, which is Chi-square
% goodness-of-fit test. Tests if 2 samples come from the same
% distribution, against the alternative that they are significant
% different
%
% The difference is that, this function deals with the data that
% has unequal weights, i.e. not 1s. The implementation follows
% procedure in Numerical Recipes Ch. 14.3.1, so that it might
% not be the same as matlab. 
%
% Given the arrays bins1[0..nbins-1] and bins2[0..nbins-1], 
% containing two sets of binned data, and given the number of
% constraints knstrn (normally 1 or 0), this routine returns the
% number of degrees of freedom df, the chi-square chsq, and the 
% p-value prob. 
%
%   A small value of prob indicates a significant 
%   difference between the distributions bins1 and bins2.
%
% Note that bins1 and bins2 are both double arrays, although 
% they will normally contain integer values.
%
% NOTE:
%   Different from 'chi2gof_wt.m', this function doesn't need 'nparam'
%   since both data sets are NOT estimated from some known distribution,
%   such that 'nparam' is always 0.
%
%
%   OUTPUT:
%       testrst:    test decision, 1 reject, 0 failure to reject
%       prob:       p-value probability p of the hypothesis testm
%       stats:      the structure stats containing information like
%                   chi-square statistic 
%   
%   INPUT:
%       bin1:       binned data set 1
%       bin2:       binned data set 2
%       binw:       bin width, for additional check for the number 
%                   of constraints
%       alpha:      significance level, to be compared with prob
%       
% By Chao Song, chaosong@princeton.edu
% First created date:   2020/10/07
% Last modified date:   2020/10/07

% vectorization 
bin1 = reshape(bin1,[],1);
bin2 = reshape(bin2,[],1);
nbins = length(bin1);   % number of bins

% the number of constraints
knstrn = 0;     % default
% additional check for the number of constraints
if sum(bin1) == sum(bin2) && sum(bin1) == 1     % i.e. the value of each bin is probability
    knstrn = 1;     % according to NR, the sum of counts for two sets are same 
elseif sum(bin1)*binw == sum(bin2)*binw && sum(bin1)*binw == 1  % i.e. the value of each bin is pdf
    knstrn = 1;     % according to NR, the sum of counts for two sets are same 
elseif sum(bin1) == sum(bin2) && sum(bin1)~=1 && sum(bin2)~=1   % pseudo counts, i.e., sum of weights in that bin
    knstrn = 1;     % according to NR, the sum of counts for two sets are same  
elseif sum(bin1) ~= sum(bin2) && sum(bin1)~=1 && sum(bin2)~=1   % pseudo counts, but not same, due to different number of detections
    knstrn = 0;
end

%  degree of freedom
df = nbins - knstrn;

% chi-square statistic formulation from NR
chisq = 0;
for i = 1: nbins
    if bin1(i) == 0 && bin2(i) == 0     % no data means one less degree of freedom
        df = df-1;
    else
        tmp = bin1(i) - bin2(i);
        chisq = chisq + tmp^2 / (bin1(i) + bin2(i));    
    
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





