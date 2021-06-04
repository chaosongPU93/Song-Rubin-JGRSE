function [opt,ort,timsSTAperm]=readpermsv3(prename,PERMSTA,PERMROTS,idx,sps,lo,hi,npo,...
                                        npa,fact)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function is an updated version of reading data from permanent
% stations. Specially for data that has already been rmean, rtrend, taper
% and removed station response. Thus the similar operations in the previous
% readperms can be deleted.
% 
% version 3, delete the find glitches part
% 
% Feartures:
%   1. read new data in data-no-resp
%   2. no taper needed, no need to multiply a factor to match the amplitude
%   between PGC and LZB in the main script.
%    
% 
%
% Modified by Chao Song, chaosong@princeton.edu
% First created date:   2020/08/13
% Last modified date:   2020/08/13
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

STAEdat=[prename,'.',PERMSTA(idx,1:3),'..BHE.D.SAC']; %HHE for polstas.
STANdat=[prename,'.',PERMSTA(idx,1:3),'..BHN.D.SAC'];
[STAE,HdrDataSTA,~,~,timsSTAperm]=readsac(STAEdat,0,'l');   %%% Function 'readsac'
[STAN,~,~,~,~]=readsac(STANdat,0,'l');

% Catch the bug that data at LZB on 2003/mar has sps of 100 instead of 40,
% but actually fix all the potential inconsistency in sps of data
spstrue = round(1./HdrDataSTA.DELTA);
if spstrue ~= 40    % permannet station should have sps of 40;
    % resample to 40 sps
    [num, denom] = rat(40/spstrue);
    STAE = resample(STAE,num,denom);   % times = num/denom
    STAN = resample(STAN,num,denom);
    timsSTAperm = linspace(round(timsSTAperm(1)), round(timsSTAperm(end)),...
        round(length(timsSTAperm)*40/spstrue));

end
    
%Filter data:   ALWAYS filter first, then resample
[STAEf]=fact*Bandpass(STAE,40,lo,hi,npo,npa,'butter');   %%% Function 'bandpass'
% It is learned from Michael Bostock that station PGC has something wrong with
% the scale of the N component. To correct the problem you should multiply the N component by 3.
if strcmp(PERMSTA(idx,1:3),'PGC')  
%     [STANf]=3*fact*Bandpass(STAN,40,lo,hi,npo,npa,'butter');
    [STANf]=fact*Bandpass(STAN,40,lo,hi,npo,npa,'butter');
else
    [STANf]=fact*Bandpass(STAN,40,lo,hi,npo,npa,'butter');
end    

% resample if needed, by Chao
[num, denom] = rat(sps/40);     % 40 is the true sampling rate of data
STAEfd = resample(STAEf,num,denom);   % times = num/denom
STANfd = resample(STANf,num,denom);
timsSTAperm = linspace(round(timsSTAperm(1)), round(timsSTAperm(end)),...
        round(length(timsSTAperm)*sps/40));

tracelen=length(STAEfd); %after decimation/upsampling

% convert the offset in 40 sps unit to the sps specified.
%%% 1st column of PERMROTS == offset in samples between fast and slow components
fastslowoff = round(PERMROTS(idx,1)*sps/40);  % this offset is based on 40sps data, even for polaris stations
%%% 4th column of PERMROTS == offset in samples between the arrival times at different stations
STAsoff = round(PERMROTS(idx,4)*sps/40);  %input offsets given at 40 sps.

%Rotate & split-correct
STA=STAEfd+1i*STANfd;
STAfastslow=STA*exp(-1i*PERMROTS(idx,2));   %%% 2th column of PERMROTS, == rotation angle to get fast/slow direction
STAslow=real(STAfastslow);
STAfast=imag(STAfastslow);
len=length(STA);
off=round(10*sps/40);
%%% 1st column of PERMROTS == offset in samples between fast/slow direction
% STAslow(off:len-off)=STAslow(off+PERMROTS(idx,1):len-off+PERMROTS(idx,1));
STAslow(off:len-off)=STAslow(off+fastslowoff: len-off+fastslowoff);
STAsplitcorrected=(STAslow+1i*STAfast)*exp(1i*PERMROTS(idx,2));
STAscrot=STAsplitcorrected*exp(-1i*PERMROTS(idx,3));

%Timeshift 
if STAsoff > -1
    STAscrot(1:tracelen-STAsoff)=STAscrot(STAsoff+1:tracelen);
    STAscrot(tracelen-STAsoff+1:tracelen)=0;
else
    STAscrot(-STAsoff+1:tracelen)=STAscrot(1:tracelen+STAsoff);
    STAscrot(1:-STAsoff)=0;
end

opt=real(STAscrot);
ort=imag(STAscrot);
