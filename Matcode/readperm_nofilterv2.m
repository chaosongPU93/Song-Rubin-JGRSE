%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Read the one day of Permanent station data, rotate according to split
% correction angle, polarization angle, to the optimal & orthogonal
% direction. version 2 for data removed from station response
% 
% Modified from readperms.m by Allan Rubin.
%
%
% By Chao Song, chaosong@princeton.edu
% Last modified date:   2019/03/02
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [opt, ort, timeperm] = readperm_nofilterv2(datafnm, PERMSTA, ...
    PERMROTS, idx, sps, fact)

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
%     STAN = 3*fact*STAN(:);
    STAN = fact*STAN(:);
else
    STAN = fact*STAN(:);
%     STAN = fact*detrend(STAN(:));
end

% resample if needed, by Chao
[num, denom] = rat(sps/40);
STAEd = resample(STAE,num,denom);   % times = num/denom
STANd = resample(STAN,num,denom);
timeperm = linspace(round(timeperm(1)), round(timeperm(end)),...
        round(length(timeperm)*sps/40));
    
tracelen = length(STAEd); %after decimation/upsampling

% convert the offset in 40 sps unit to the sps specified.
%%% 1st column of PERMROTS == offset in samples between fast and slow components
fastslowoff = round(PERMROTS(idx,1)*sps/40);  % this offset is based on 40sps data, even for polaris stations
%%% 4th column of PERMROTS == offset in samples between the arrival times at different stations
STAsoff = round(PERMROTS(idx,4)*sps/40);  %input offsets given at 40 sps.

%Rotate & split-correct
STA = STAEd+ 1i* STANd;
STAfastslow = STA* exp(-1i* PERMROTS(idx,2));   %%% 2th column of PERMROTS, == rotation angle to get fast/slow direction
STAslow = real(STAfastslow);
STAfast = imag(STAfastslow);
len = length(STA);
off = round(10*sps/40);
%%% 1st column of PERMROTS == offset in samples between fast/slow direction
STAslow(off: len-off) = STAslow(off+ fastslowoff: len- off+ fastslowoff); 
STAsplitcorrected = (STAslow+ 1i* STAfast)* exp(1i* PERMROTS(idx,2));
STAscrot = STAsplitcorrected* exp(-1i* PERMROTS(idx,3));

%Timeshift
if STAsoff > -1
    STAscrot(1: tracelen- STAsoff) = STAscrot(STAsoff+ 1: tracelen);
    STAscrot(tracelen- STAsoff+ 1: tracelen) = 0;
else
    STAscrot(-STAsoff+ 1: tracelen) = STAscrot(1: tracelen+ STAsoff);
    STAscrot(1: -STAsoff) = 0;
end

opt = real(STAscrot);
ort = imag(STAscrot);
