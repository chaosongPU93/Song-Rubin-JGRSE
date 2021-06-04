%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Read the one day of Polaris station data, rotate according to split
% correction angle, polarization angle, to the optimal & orthogonal
% direction. version 2 for data removed from station response
% 
% USE data removed from station response 
%
% Modified from readpols.m by Allan Rubin.
%
%
% By Chao Song, chaosong@princeton.edu
% Last modified date:   2019/06/29
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [STAEd, STANd, timepola] = readpols_norotsv2(datafnm, POLSTA, idx, sps, fact)

STAEdat = [datafnm,'.',POLSTA(idx,1:4),'..HHE.D.SAC']; %HHE for polstas.
STANdat = [datafnm,'.',POLSTA(idx,1:4),'..HHN.D.SAC'];
[STAE, ~, ~, ~, timepola] = readsac(STAEdat, 0, 'l');   %%% Function 'readsac'
[STAN, ~, ~, ~, ~] = readsac(STANdat,0, 'l');    %%% 'l' means linux

%%% SCALE
STAE = fact* detrend(STAE(:)); 
STAN = fact* detrend(STAN(:));

% resample if needed, by Chao
[num, denom] = rat(sps/100);
STAEd = resample(STAE,num,denom);   % times = num/denom
STANd = resample(STAN,num,denom);
timepola = linspace(round(timepola(1)), round(timepola(end)),...
        round(length(timepola)*sps/40));


