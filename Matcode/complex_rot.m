function [newx,newy] = complex_rot(x0,y0,angle)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function is used to rotate the original point by an angle, and output the coordinates in the
% original system  
%
% INPUT:  
%   x0:     original x
%   y0:     original y
%   angle:  calculated counter-clockwise. in degree, add '-' to incate clockwise rotation
% 
% OUTPUT:
%   newx:   rotated x in the same coordinate system
%   newy:   rotated y in the same coordinate system
%
% Chao Song, chaosong@princeton.edu
% First created date:   2019/10/02
% Last modified date:   2019/10/02
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

c = x0+1i.*y0;
rotc = c.*exp(1i * deg2rad(angle));   % rotate counter-clockwise by angle
newx = real(rotc);
newy = imag(rotc);