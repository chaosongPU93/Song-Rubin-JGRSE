function c = whitered(m)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%WHITERED    Shades of red color map
%   WHITERED(M), is an M-by-3 matrix that defines a colormap.
%   The colors begin with white, and then through shades of 
%   red to bright red.
%   WHITERED, by itself, is the same length as the current figure's
%   colormap. If no figure exists, MATLAB creates one.
%
%   For example, to reset the colormap of the current figure:
%
%             colormap(whitered)
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
 
% From [1 1 1] to [1 0 0] 
r = ones(m,1);
g = (m-1:-1:0)'/max(m-1,1);
b = g;

c = [r g b]; 
