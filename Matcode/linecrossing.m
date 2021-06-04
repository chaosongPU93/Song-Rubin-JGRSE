function [x0,y0] = linecrossing(a1,b1,a2,b2)
% y = linefcn(x,a,b)
%
% This function outputs the crossing point (x0,y0) based on
% the input slope a and intercept b of 2 lines
%
% y = a1.* x + b1
% y = a2.* x + b2
%
% combine this 2 equations to solve for the crossing point x0,y0
%
% By Chao Song, chaosong@princeton.edu
% First created date:   2020/11/18
% Last modified date:   2020/11/18

if a1 == a2 && b1 == b2
    disp('Two lines overlap.');
    x0 = [];
    y0 = [];
elseif a1 == a2 && b1 ~= b2
    disp('Two lines are parallel.');
    x0 = [];
    y0 = [];
else
    x0 = (b2-b1)/(a1-a2);
    % y0 = (b2*a1-b1*a2)/(a1-a2);
    y0 = linefcn(x0,a1,b1);
end