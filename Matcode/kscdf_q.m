function pks = kscdf_q(z)
% KSCDF_Q Kolmogorov-Smirnov complementary cumulative distribution function (cdf) Q_KS(z).
% 
% This function is to compute the complementary Kolmogorov-Smirnov cumulative
% distribution function. Its probability density function does not directly
% enter into the test and is virtually never even written. Q_KS(z)=1-P_KS(z)
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
    pks = 1.0;
end

if z < 1.18
    pks = 1.0- kscdf_p(z);
    
else
    x = exp(-2.0*z^2);
    pks = 2.0*(x-x^4+x^9);
end