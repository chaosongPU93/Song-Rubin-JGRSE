function [scrsz, res] = pixelperinch(imo)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This function is to get how many pixels per inch of your monitor so that
% when you plot figures, you can easily convert between to make sure the
% saved pictures, especially pdfs won't exceed the real paper margin 
%
%
% By Chao Song, chaosong@princeton.edu
% First created date:   2020/02/10
% Last modified date:   2020/02/10
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% get the monitor size in inches
set(0,'units','inches');
tmp1 = get(0,'MonitorPositions');
scrszin = tmp1(imo,:);

% get the monitor size in pixels
set(0,'units','pixels'); 
tmp2 = get(0,'MonitorPositions');
scrsz = tmp2(imo,:);

res = scrsz(4)/scrszin(4);












