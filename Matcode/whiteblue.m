function c = whiteblue(m)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%WHITEBLUE    Shades of blue color map
%   WHITEBLUE(M), is an M-by-3 matrix that defines a colormap.
%   The colors begin with white, and then through shades of 
%   blue to bright blue.
%   WHITEBLUE, by itself, is the same length as the current figure's
%   colormap. If no figure exists, MATLAB creates one.
%
%   For example, to reset the colormap of the current figure:
%
%             colormap(whitebue)
%
%   See also HSV, GRAY, HOT, BONE, COPPER, PINK, FLAG, 
%   COLORMAP, RGBPLOT.
%
% Chao Song, chaosong@princeton.edu
% Inspired by REDBLUE.m 
% First created date:   2021/06/02
% Last modified date:   2021/06/02
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 1, m = size(get(gcf,'colormap'),1); end

% From [1 1 1] to [0 0 1] 
r = (m-1:-1:0)'/max(m-1,1);
g = r;
b = ones(m,1);

c = [r g b];