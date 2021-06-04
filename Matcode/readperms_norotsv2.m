%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Read the one day of Permanent station data without station response,
% read n and e component, ready to rotate to get rots parameter
% 
% USE data removed from station response 
%
% By Chao Song, chaosong@princeton.edu
% Last modified date:   2019/07/08
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [STAEd, STANd, timeperm] = readperms_norotsv2(datafnm, PERMSTA, idx, sps, fact)

STAEdat = [datafnm,'.',PERMSTA(idx,1:3),'..BHE.D.SAC']; %HHE for polstas.
STANdat = [datafnm,'.',PERMSTA(idx,1:3),'..BHN.D.SAC'];
[STAE, HdrDataSTA, ~, ~, timeperm] = readsac(STAEdat, 0, 'l');   %%% Function 'readsac'
[STAN, ~, ~, ~, ~] = readsac(STANdat,0, 'l');    %%% 'l' means linux

% Catch the bug that data at LZB on 2003/mar has sps of 100 instead of 40,
% but actually fix all the potential inconsistency in sps of data
spstrue = round(1./HdrDataSTA.DELTA);
if spstrue ~= 40    % permannet station should have sps of 40;
    % resample to 40 sps
    [num, denom] = rat(40/spstrue);
    STAE = resample(STAE,num,denom);   % times = num/denom
    STAN = resample(STAN,num,denom);
    timeperm = linspace(round(timeperm(1)), round(timeperm(end)),...
        round(length(timeperm)*40/spstrue));
end

%%% SCALE
STAE = fact* detrend(STAE(:)); 
if strcmp(PERMSTA(idx,1:3),'PGC')
%     STAN = 3*fact*detrend(STAN(:)); 
    STAN = fact*detrend(STAN(:));

else
    STAN = fact*detrend(STAN(:)); 
end

% resample if needed, by Chao
[num, denom] = rat(sps/40);
STAEd = resample(STAE,num,denom);   % times = num/denom
STANd = resample(STAN,num,denom);
timeperm = linspace(round(timeperm(1)), round(timeperm(end)),...
        round(length(timeperm)*sps/40));
   
