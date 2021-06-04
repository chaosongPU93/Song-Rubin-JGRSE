function pks = kscdf_p(z)
% KSCDF_P Kolmogorov-Smirnov cumulative distribution function (cdf) P_KS(z).
% 
% This function is to compute the Kolmogorov-Smirnov cumulative distribution
% function. Its probability density function does not directly enter into the
% test and is virtually never even written.
% 
% mainly just a direct transformation from Numerical Recipes Ch. 6.14.12
%
%       
% By Chao Song, chaosong@princeton.edu
% First created date:   2020/10/07
% Last modified date:   2020/10/07

if z < 0
    warning('bad input, has to be positive');
end

if z == 0
    pks = 0;
end

if z < 1.18
    y = exp(-1.23370055013616983/(z^2));
    pks = 2.25675833419102515*sqrt(-log(y))*(y + y^9 + y^25 + y^49);
    
else
    x = exp(-2.0*z^2);
    pks = 1.0 - 2.0*(x-x^4+x^9);
end