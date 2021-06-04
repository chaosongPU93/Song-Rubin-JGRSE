function [opt,ort,nzeros]=readpols(prename,POLSTA,POLROTS,idx,sps,lo,hi,npo,npa,fact,nwin,winlen,winoff,igstart)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
% Altered 07/03/18 to shift the stations before counting consecutive zeros.
% (For traveltime but not splitting.)  Z comp. not shifted.

STAsoffspsorig=round(POLROTS(idx,4)*100/40)  %input offsets given at 40 sps, POLARIS stations record at 100.
if(POLSTA(idx,5:5)==' ')   %%% The name has 4 charcters
    STAEdat=[prename,'.',POLSTA(idx,1:4),'..HHE.D.SAC']; %BHE for permstas.
    STANdat=[prename,'.',POLSTA(idx,1:4),'..HHN.D.SAC'];
    [STAE,HdrDataSTA,tnuSTA,pobjSTA,timsSTApol]=readsac(STAEdat,0,'l');
    [STAN,~,~,~,~]=readsac(STANdat,0,'l');
    tracelen=length(STAE); %temporary; before decimation

    %find glitches
    nzerosE=glitchesPOL(STAE,nwin,winlen,winoff,igstart,sps,STAsoffspsorig);
    nzerosN=glitchesPOL(STAN,nwin,winlen,winoff,igstart,sps,STAsoffspsorig);
    nzeros=max(nzerosE,nzerosN);

    %cosine taper over 1 second (at 100 sps) before filtering:
    x=(0:pi/200:pi/2-pi/200)';
    STAE(1:100)=sin(x).*STAE(1:100);
    STAN(1:100)=sin(x).*STAN(1:100);
    x=flipud(x);
    STAE(tracelen-99:tracelen)=sin(x).*STAE(tracelen-99:tracelen);
    STAN(tracelen-99:tracelen)=sin(x).*STAN(tracelen-99:tracelen);

    %Filter data:
    [STAEf]=fact*Bandpass(STAE,100,lo,hi,npo,npa,'butter');
    [STANf]=fact*Bandpass(STAN,100,lo,hi,npo,npa,'butter');
    
    %%%%%%% the following was Allan's code, commented out by Chao for a better
    %%%%%%% expression
%     if sps == 100
%         %no decimation
%         STAEfd=STAEf;
%         STANfd=STANf;
%     else
%         %decimate the 100 sps data to 40 (assume one or the other)
%         STAEfd=resample(STAEf,2,5);
%         STANfd=resample(STANf,2,5);
%     end
    %%%%%%%%
    
    % resample if needed, by Chao
    [num, denom] = rat(sps/100);
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
    [STAZ,HdrDataSTA,tnuSTA,pobjSTA,timsSTApol]=readsac(STAZdat,0,'l');
    tracelen=length(STAZ); %temporary; before decimation
    %find glitches
    nzeros=glitchesPOL(STAZ,nwin,winlen,winoff,igstart,sps,0); %Last variable is offset; make it zero for vertical component.
    %cosine taper over 1 second (at 100 sps) before filtering:
    x=(0:pi/200:pi/2-pi/200)';
    STAZ(1:100)=sin(x).*STAZ(1:100);
    x=flipud(x);
    STAZ(tracelen-99:tracelen)=sin(x).*STAZ(tracelen-99:tracelen);
    %Filter data:
    [STAZf]=fact*bandpass(STAZ,100,lo,hi,npo,npa,'butter');
    %%%%%%% the following was Allan's code, commented out by Chao for a better
    %%%%%%% expression
%     if sps == 100
%         %no decimation
%         STAZfd=STAZf;
%     else
%         %decimate the 100 sps data to 40 (assume one or the other)
%         STAZfd=resample(STAZf,2,5);
%     end
    %%%%%%%%
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
