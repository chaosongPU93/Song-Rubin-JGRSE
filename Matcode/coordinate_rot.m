function [newx,newy] = coordinate_rot(x0,y0,angle,sftx,sfty)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function is used to rotate the coodinate system by angle to get the new expression of points
% in the new coord system
%
% INPUT:  
%   x0:     original x, one point or vector
%   y0:     original y, one point
%   angle:  calculated counter-clockwise. in degree, add '-' to incate clockwise rotation
%   sftx:   shift of x in new coord
%   stfy:   shift of y in new coord
% 
% OUTPUT:
%   newx:   rotated x in the new coordinate system
%   newy:   rotated y in the new coordinate system
%
% Chao Song, chaosong@princeton.edu
% First created date:   2019/10/02
% Last modified date:   2019/10/02
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

radang = deg2rad(angle);
R = [cos(radang) sin(radang); -sin(radang) cos(radang)];    % rotation matrix
x0 = reshape(x0,1,[]);
y0 = reshape(y0,1,[]);
newcoord = R * [x0; y0];    % matrix multiplication
newx = newcoord(1,:)'+sftx;
newy = newcoord(2,:)'+sfty;