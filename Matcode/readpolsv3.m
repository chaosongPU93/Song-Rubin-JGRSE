function [opt,ort]=readpolsv3(prename,POLSTA,POLROTS,idx,sps,lo,hi,npo,npa,fact,nwin,winlen,winoff,igstart)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function is an updated version of reading data from polaris
% stations. Specially for data that has already been rmean, rtrend, taper
% and removed station response. Thus the similar operations in the previous
% readpols can be deleted.
% 
% version 3, delete the find glitches part
%
% Feartures:
%   1. read new data in data-no-resp
%   2. no taper needed, a factor is still needed, since in year 2003 and 
%      after, amplitude in the polaris station still has large difference
%      which is not shown in permanent stations, even after removing
%      station response.
%    
% 
%
% Modified by Chao Song, chaosong@princeton.edu
% First created date:   2020/08/13
% Last modified date:   2020/08/13
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

STAsoffspsorig=round(POLROTS(idx,4)*100/40);  %input offsets given at 40 sps, POLARIS stations record at 100.
if(POLSTA(idx,5:5)==' ')   %%% The name has 4 charcters
    STAEdat=[prename,'.',POLSTA(idx,1:4),'..HHE.D.SAC']; %BHE for permstas.
    STANdat=[prename,'.',POLSTA(idx,1:4),'..HHN.D.SAC'];
    [STAE,~,~,~,~]=readsac(STAEdat,0,'l');
    [STAN,~,~,~,~]=readsac(STANdat,0,'l');

    %Filter data:
    [STAEf]=fact*Bandpass(STAE,100,lo,hi,npo,npa,'butter');
    [STANf]=fact*Bandpass(STAN,100,lo,hi,npo,npa,'butter');
    
    % resample if needed, by Chao
    [num, denom] = rat(sps/100);    % 100 is the true sampling rate of data
    STAEfd = resample(STAEf,num,denom);   % times = num/denom
    STANfd = resample(STANf,num,denom);
    
    tracelen=length(STAEfd); %after decimation
    
    % convert the offset in 40 sps unit to the sps specified.
    %%% 1st column of PERMROTS == offset in samples between fast and slow components
    fastslowoff = round(POLROTS(idx,1)*sps/40);  % this offset is based on 40sps data, even for polaris stations
    %%% 4th column of PERMROTS == offset in samples between the arrival times at different stations
    STAsoff = round(POLROTS(idx,4)*sps/40);  %input offsets given at 40 sps.
    
    %Rotate & split-correct
    STA=STAEfd+1i*STANfd;
    STAfastslow=STA*exp(-1i*POLROTS(idx,2));
    STAslow=real(STAfastslow);
    STAfast=imag(STAfastslow);
    len=length(STA);
    off=round(10*sps/40);
%     STAslow(off:len-off)=STAslow(off+POLROTS(idx,1):len-off+POLROTS(idx,1));
    STAslow(off:len-off)=STAslow(off+fastslowoff: len-off+fastslowoff);
    STAsplitcorrected=(STAslow+1i*STAfast)*exp(1i*POLROTS(idx,2));
    STAscrot=STAsplitcorrected*exp(-1i*POLROTS(idx,3));

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
    
    
else    %%% The name has 5 characters??? ONLY read Z component
    STAZdat=[prename,'.',POLSTA(idx,1:4),'..HHZ.D.SAC']; %BHE for permstas.
    [STAZ,~,~,~,~]=readsac(STAZdat,0,'l');
    
    %find glitches
    nzeros=glitchesPOL(STAZ,nwin,winlen,winoff,igstart,sps,0); %Last variable is offset; make it zero for vertical component.
    
    %Filter data:
    [STAZf]=fact*bandpass(STAZ,100,lo,hi,npo,npa,'butter');

    % resample if needed, by Chao
    [num, denom] = rat(sps/100);
    STAZfd=resample(STAZf,num,denom);   % times = num/denom
    tracelen=length(STAZfd); %after decimation
    
    %Timeshift
    STAsoff=round(POLROTS(idx,4)*sps/40);  %input offsets given at 40 sps
    if STAsoff > -1
        STAZfd(1:tracelen-STAsoff)=STAZfd(STAsoff+1:tracelen);
        STAZfd(tracelen-STAsoff+1:tracelen)=0;
    else
        STAZfd(-STAsoff+1:tracelen)=STAZfd(1:tracelen+STAsoff);
        STAZfd(1:-STAsoff)=0;
    end
    opt=STAZfd;
    ort=0*opt;
end
