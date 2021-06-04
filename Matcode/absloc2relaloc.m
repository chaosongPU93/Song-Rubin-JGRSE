function [dx, dy] = absloc2relaloc(lon,lat,lon0,lat0) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function is used to convert the longitude and latitude of detections
% inverted by hypoinverse to locations in km relative to offset (0,0)
% 
%
% Modified by Chao Song, chaosong@princeton.edu
% First created date:   2019/08/20
% Last modified date:   2019/08/20
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rads=pi/180.;
erad=6372.028;
srad=erad*cos(lat0*rads);
dy=rads*(lat-lat0)*erad;    % N--S
dx=rads*(lon-lon0)*srad;   % W--E