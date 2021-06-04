% function plotaddcheck_hf_LZB_exp
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% THIS is used to plot the relationship of the offset 1-2, 1-3 & 1-4/5, ...
% station 4/5, ... is to check the high freq (1.25-6.5hz) detections by the 
% 3-sta trio LZB
%
% This is created for the migration example #13, since we want to lower the 
% the HF detection CC threshold to redo the detection for day 2004 198 (Jul 14),
% in order to see if during the earlier time we can see HF. This is to answer
% why we missed HF in the current detection regime (2021/01/26)
%
% DIFFerence with low freq is mainly on the detection parameters!!
%
% INPUT:
%   ifam:       fam number
%   hitallow:   maximum hit counts at one location
%   fcflag:     flag to do the filtering effect correction
%   convflag:   flag to convert all offsets to same ref fam
%   iup:        number of times that offsets are upsampled
%
% Features:
%   1. x-axis: 1-2 offset
%      y-axis: 1-3 offset
%      z-axis: color coded by 1-4 offset
%   2. scatter only the offsets that meet the criteria in station adding
%       check script 'hf_LZBtrio_detection_check.m'
%   3. save results
%   4. add option flag to correct filtering effect (1/0)
%   5. add option flag to convert all fams to the same frame (1/0)
%
% NOTES:
%   1. NO set 'countmin' to the original detections (set to 1), 2019/08/23 
%   2. fcflag is usually set to be 1, while being 0 in merge_fam
%   3. if convflag is set to be 1, then will convert all fam to the same frame,
%       so that it is possible to concatenate all rsts to generate a summary phase
%       input files for all fam in hypoinverse; 
%       if set to be 0, then individual operations need to be done. Now (2019/11/20)
%       the conversion algorithm is checked to be true, so it is safe to set this
%       flag to 1.
%
% Conversion algorithm:
%       (off12')_fam - (off12')_ref = (off12)_fam - (off12)_ref
%   where ' means the detection in that family, no ' means the original
%   reference point in each family, i.e. the location of that family 
%   (considered as 0,0 in each fam). Since we already have the travel time
%   difference between sta 12 & 13 in each family (4th col in rots para), 
%   what we need to do here is to convert them to a uniform frame.
%   Thus, relative to a selected fam, all detections in other fam can be
%   written in the new frame as:
%       (off12')_ref = (off12')_fam - [(off12)_fam - (off12)_ref]
%


% Chao Song, chaosong@princeton.edu
% First created date:   2021/01/26
% Last modified date:   2021/01/26
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% %% default value for easy debugging
% defval('ifam', 11);
% defval('hitallow', 400);
% defval('fcflag',1);
% defval('convflag',0);
% defval('iup',4);


%% Initialization
%%% SAME if focusing on the same region (i.e. same PERMROTS and POLROTS)
%%% AND if using the same family, same station trio

format short e   % Set the format to 5-digit floating point
clear
close all

% set(0,'DefaultFigureVisible','on');
set(0,'DefaultFigureVisible','off');   % switch to show the plots or not

scrsz=get(0,'ScreenSize');

% WHEN CHANGING FAMILIES CHANGE: (1)dates (2)Bostnames (3)hilo,frequency band 
% (4)mshift (5)bostsec (6)stas (7)PERMROTS and POLROTS (8)tempoffs

workpath = getenv('ALLAN');
%%% choose the data path here 
% datapath = workpath;
datapath = strcat(workpath,'/data-no-resp');
temppath = strcat(datapath, '/templates/LZBtrio/');
rstpath = strcat(datapath, '/LZBtrio');

nfampool = ['002';
            '043';
            '141';
            '047';
            '010';
            '144';
            '099';
            '068';
            '125';
            '147';
            '017'];

% combine all fams into an pool, 12 fams in total, final version        
nfampool = ['002';
            '043';
            '141';
            '047';
            '010';
            '144';
            '099';
            '068';
            '125';
            '147';
            '017';
            '006'];
        
nfampool = ['002';
            '043';
            '141';
            '047';
            '010';
            '144';
            '099';
            '068';
            '125';
            '147';
            '017';
            '006';
            '001';
%             '158';      % 158, 20200916,testing purpose
            ];     % 2020/09/07, i am adding a new family 006 to the pool, final version
        
freqflag='hf';  % flag to indicate whether to do hf or lf;

% get days, 3-sta trio name, bostock's template name, zero crossing index
% of templates
timoffrot= [
            2004 198;
            2005 256;
            2005 259;
           ];
        
nday = size(timoffrot, 1);

%%% corresponds to PERMROTS
PERMSTA=['PGC'        % permanent station names
         'LZB'];
POLSTA=['SSIB '           % polaris station names
        'SILB '
        'KLNB '
        'MGCB '
        'TWKB '];

stas=['TWKB '
      'LZB  '
      'MGCB '];     % determine the trio and order

for ifam = 1: size(nfampool,1) 
    
    fam = nfampool(ifam, :);
    hitallow = 400;  % hitallow to 300 is efficient to all fam,
    fcflag = 1;  % set fcflag=1
    convflag = 0;   % convflag=0
    iup = 4;    % upsample number,times
    
%% Important parameters same as that during detections    
nsta=size(stas,1);         %  number of stations
sps=40;     % samples per second

%%% IMPORTANT, NEED to change sometimes
%Basics of the cross-correlation:  Window length, number of windows, filter parameters, etc.
winlensec=4;
winoffsec=1;        % window offset in sec, which is the step of a moving window
winlen=winlensec*sps;      % length in smaples
hi=6.5;
lo=1.25;
npo=2;     % poles, passes of filters
npa=2;
cyclskip = 0;
mshift=14+cyclskip; %19; %maximum shift for the x-correlations. 19 for 002 Stanford,    % in sps, 0.5s*40sps=20
loopoffmax=2.5; %1.5 for standard 1.5-6Hz; 4 for 0.5-1.5Hz.  2 for non-interpolated.   % what is loopoffmax, the circuit of time offsets
xcmaxAVEnmin=0.01; %0.44; %0.44 for 002 Stanford %0.45; %0.36 for 4s 1-12 Hz; %0.4 for 4s 1.5-6 Hz and 6s 0.5-1.5Hz; 0.36 for 4s2-8 Hz ; 0.38 for 4s0.75-6 Hz; 0.094 for 128s 2-8Hz;  0.1 for 128s 1.5-6Hz; 0.44 for 3-s window?

upmshift = mshift*iup; 
matrange = 2*upmshift+1; 

%     hitallow = 230;     % max number of hit allowed to iniate the matrix
offmat14 = NaN(matrange, matrange, hitallow);     % +mshiftnew+1 is just for easy visuallization
countmat14 = zeros(matrange, matrange);     % count of hit at the each spot
meanmat14 = NaN(matrange, matrange);  % mean of offset14 at each spot,
stdmat14 = NaN(matrange, matrange);   % standard deviation of offset14 at each spot
tbgmat14 = zeros(matrange, matrange, hitallow);
tctmat14 = zeros(matrange, matrange, hitallow);
tdatemat14 = zeros(matrange, matrange, hitallow);
ccmat14 = zeros(matrange, matrange, hitallow);


offmat15 = NaN(matrange, matrange, hitallow);
countmat15 = zeros(matrange, matrange);
meanmat15 = NaN(matrange, matrange, hitallow);
stdmat15 = NaN(matrange, matrange, hitallow);
tbgmat15 = zeros(matrange, matrange, hitallow);
tctmat15 = zeros(matrange, matrange, hitallow);
tdatemat15 = zeros(matrange, matrange, hitallow);
ccmat15 = zeros(matrange, matrange, hitallow);


offmat16 = NaN(matrange, matrange, hitallow);
countmat16 = zeros(matrange, matrange);
meanmat16 = NaN(matrange, matrange, hitallow);
stdmat16 = NaN(matrange, matrange, hitallow);
tbgmat16 = zeros(matrange, matrange, hitallow);
tctmat16 = zeros(matrange, matrange, hitallow);
tdatemat16 = zeros(matrange, matrange, hitallow);
ccmat16 = zeros(matrange, matrange, hitallow);


offmat17 = NaN(matrange, matrange, hitallow);
countmat17 = zeros(matrange, matrange);
meanmat17 = NaN(matrange, matrange, hitallow);
stdmat17 = NaN(matrange, matrange, hitallow);
tbgmat17 = zeros(matrange, matrange, hitallow);
tctmat17 = zeros(matrange, matrange, hitallow);
tdatemat17 = zeros(matrange, matrange, hitallow);
ccmat17 = zeros(matrange, matrange, hitallow);


countmat23 = zeros(matrange, matrange);
tbgmat23 = zeros(matrange, matrange, hitallow);
tctmat23 = zeros(matrange, matrange, hitallow);
tdatemat23 = zeros(matrange, matrange, hitallow);
ccmat23 = zeros(matrange, matrange, hitallow);

if fcflag
    %%% load the template filtering effect result
    lolf = 0.5;
    hilf = 1.25;
    lohf = 1.25;
    hihf = 6.5;
    winsechf = winlensec;
    winseclf = 16;
    lofflf = 4;
    ccminlf = 0.35;
    PREFIX = strcat(fam,'.lf',num2str(lolf),'-',num2str(hilf),'.hf',num2str(lohf),'-',num2str(hihf),...
        '.hw',num2str(winsechf),'.lw',num2str(winseclf),'.lo',num2str(lofflf),'.ccm', ...
        num2str(ccminlf),'.','80sps');
    fname = strcat(rstpath, '/MAPS/tempfeff_LZB_enc','_',PREFIX);
    filtcor = load(fname);
    spsratio = sps/80;
    filthf = filtcor(:,1)*spsratio;   % lf template shift due to filterig effect, sign is the same
else
    filthf = zeros(2,1);
end


%% START TO LOOP FOR every day
%cycle over each day:
for nd=1:length(timoffrot(:,1))      % num of rows, also num of days

    %Which days of data to read?
    year=timoffrot(nd,1);
    YEAR=int2str(year);
    jday=timoffrot(nd,2);
    if jday <= 9
        JDAY=['00',int2str(jday)];
    elseif jday<= 99
        JDAY=['0',int2str(jday)];
    else
        JDAY=int2str(jday);
    end
    
    IDENTIF=[YEAR,'.',JDAY,'.',fam,'.lo',num2str(loopoffmax),'.cc',num2str(xcmaxAVEnmin),'.',...
             int2str(npo),int2str(npa),'.ms', int2str(mshift)]    

%%%%%% Uncomment below WHEN using merge_fam_befhypo.m, means removing duplicates BEFORE hypo %%%%%%%        
%     if iup == 1     
%         fname1 = strcat(rstpath,'/MAPS/mapallrstNoDou',IDENTIF,'_',num2str(lo),'-',num2str(hi),'_',...
%                         num2str(winlen/sps),'s',num2str(sps),'sps', '4nsta');   
%         fname2 = strcat(rstpath,'/MAPS/mapallrstNoDou',IDENTIF,'_',num2str(lo),'-',num2str(hi),'_',...
%                         num2str(winlen/sps),'s',num2str(sps),'sps', '3nsta');
%     else
%         fname1 = strcat(rstpath,'/MAPS/mapallrstUpNoDou',IDENTIF,'_',num2str(lo),'-',num2str(hi),'_',...
%                         num2str(winlen/sps),'s',num2str(sps),'sps', '4nsta');   
%         fname2 = strcat(rstpath,'/MAPS/mapallrstUpNoDou',IDENTIF,'_',num2str(lo),'-',num2str(hi),'_',...
%                         num2str(winlen/sps),'s',num2str(sps),'sps', '3nsta');
%     end
%%%%%% Uncomment above WHEN using merge_fam_befhypo.m, means removing duplicates BEFORE hypo %%%%%%%

%%%%%% Uncomment below WHEN using merge_fam_afthypo.m, means removing duplicates AFTER hypo %%%%%%%        
    if iup == 1     
        fname1 = strcat(rstpath,'/MAPS/mapallrst',IDENTIF,'_',num2str(lo),'-',num2str(hi),'_',...
                        num2str(winlen/sps),'s',num2str(sps),'sps', '4nsta');   
        fname2 = strcat(rstpath,'/MAPS/mapallrst',IDENTIF,'_',num2str(lo),'-',num2str(hi),'_',...
                        num2str(winlen/sps),'s',num2str(sps),'sps', '3nsta');
    else
        fname1 = strcat(rstpath,'/MAPS/mapallrst_up',IDENTIF,'_',num2str(lo),'-',num2str(hi),'_',...
                        num2str(winlen/sps),'s',num2str(sps),'sps', '4nsta');   
        fname2 = strcat(rstpath,'/MAPS/mapallrst_up',IDENTIF,'_',num2str(lo),'-',num2str(hi),'_',...
                        num2str(winlen/sps),'s',num2str(sps),'sps', '3nsta');
    end
%%%%%% Uncomment above WHEN using merge_fam_befhypo.m, means removing duplicates AFTER hypo %%%%%%%
    
    if isfile(fname1)
        allrst_new = load(fname1);
        %%% additional stations for checking detections
        stasnew=['PGC  '
             'SSIB '
             'SILB '
             'KLNB '];
        if year==2003   % in 2003, station KLNB is named KELB, so use KELB to replace KLNB in 2003
            POLSTA(3,:)='KELB ';
            stasnew(4,:)='KELB ';
        else
            POLSTA(3,:)='KLNB ';  % remember to change it back
            stasnew(4,:)='KLNB ';
        end
        nstanew=size(stasnew,1);
    elseif isfile(fname2)
        allrst_new = load(fname2);
        stasnew=['PGC  '
             'SSIB '
             'SILB '];
        nstanew=size(stasnew,1);
    else
        fprintf('Day %s %s of fam %s is skipped because of no detections.\n',YEAR, JDAY, fam);
    end 

    %%%%% this part is only for original 3 station trio before adding check
    % load those offsets
    ndSTA12off = allrst_new(:, 3);    % 3rd col of mapfile
    ndSTA13off = allrst_new(:, 2);    % 2nd col of mapfile
    ndtct23 = allrst_new(:, 1);    % timswin == time at the center of each win
    ndtbg23 = allrst_new(:, 7);  % begin time in sec of the strongest arrival diapole
    ndcc23 = allrst_new(:, 4);  % ave cc coef of original trio
    
    % upsample the offset to make sure it is integer
    ndSTA12off = ndSTA12off*iup;
    ndSTA13off = ndSTA13off*iup;

    date = floor(str2double(strcat(YEAR,JDAY)));
    
    for i = 1: length(ndSTA12off)
        
        countmat23(ndSTA13off(i)+upmshift+1, ndSTA12off(i)+upmshift+1) = ...
            countmat23(ndSTA13off(i)+upmshift+1, ndSTA12off(i)+upmshift+1)+ 1;
        % temp is used to indicate the times one spot has been hit
        temp = countmat23(ndSTA13off(i)+upmshift+1, ndSTA12off(i)+upmshift+1);
        % temp is also the index of the 3rd dimension
        tctmat23(ndSTA13off(i)+upmshift+1, ndSTA12off(i)+upmshift+1, temp) = ...
            ndtct23(i);
        
        tbgmat23(ndSTA13off(i)+upmshift+1, ndSTA12off(i)+upmshift+1, temp) = ...
            ndtbg23(i);
        
        tdatemat23(ndSTA13off(i)+upmshift+1, ndSTA12off(i)+upmshift+1, temp) = ...
            date;
        
        ccmat23(ndSTA13off(i)+upmshift+1, ndSTA12off(i)+upmshift+1, temp) = ...
            ndcc23(i);
        
    end
    %%%%%

    % get the 3D offset and count matrix
    [temp_count14,temp_off14,temp_tct14,temp_tbg14,temp_tdate14,temp_cc14] = ...
        Addcheck_offcount(1,nstanew,date,upmshift,allrst_new,countmat14,offmat14,tctmat14,tbgmat14,...
        tdatemat14,ccmat14,iup);
    countmat14 = temp_count14;
    offmat14 = temp_off14;
    tctmat14 = temp_tct14;
    tbgmat14 = temp_tbg14;
    tdatemat14 = temp_tdate14;
    ccmat14 = temp_cc14;
    
    [temp_count15,temp_off15,temp_tct15,temp_tbg15,temp_tdate15,temp_cc15] = ...
        Addcheck_offcount(2,nstanew,date,upmshift,allrst_new,countmat15,offmat15,tctmat15,tbgmat15,...
        tdatemat15,ccmat15,iup);
    countmat15 = temp_count15;
    offmat15 = temp_off15;
    tctmat15 = temp_tct15;
    tbgmat15 = temp_tbg15;
    tdatemat15 = temp_tdate15;
    ccmat15 = temp_cc15;
    
    if nstanew > 2
        
        [temp_count16,temp_off16,temp_tct16,temp_tbg16,temp_tdate16,temp_cc16] = ...
            Addcheck_offcount(3,nstanew,date,upmshift,allrst_new,countmat16,offmat16,tctmat16,...
            tbgmat16,tdatemat16,ccmat16,iup);
        countmat16 = temp_count16;
        offmat16 = temp_off16;
        tctmat16 = temp_tct16;
        tbgmat16 = temp_tbg16;
        tdatemat16 = temp_tdate16;
        ccmat16 = temp_cc16;
        
        if nstanew > 3
            
            [temp_count17,temp_off17,temp_tct17,temp_tbg17,temp_tdate17,temp_cc17] = ...
                Addcheck_offcount(4,nstanew,date,upmshift,allrst_new,countmat17,offmat17,tctmat17,...
                tbgmat17,tdatemat17,ccmat17,iup);
            countmat17 = temp_count17;
            offmat17 = temp_off17;
            tctmat17 = temp_tct17;
            tbgmat17 = temp_tbg17;
            tdatemat17 = temp_tdate17;
            ccmat17 = temp_cc17;
            
        end
        
    end
    
end

% re-define stasnew and nstanew
stasnew=['PGC  '
         'SSIB '
         'SILB '
         'KLNB '];  % in 2003, KLNB was named KELB
nstanew=size(stasnew,1);

%%
maxvect = [max(max(countmat14)),max(max(countmat15)),max(max(countmat16)),max(max(countmat17)),...
           max(max(countmat23))];
hitmax = max(maxvect);
if hitmax > hitallow
    fprintf(strcat("ERROR, need to increase the 'hitallow' at least to ", num2str(hitmax), '\n'))
    return;
end

%%%%% this part is only for original 3 station trio before adding check
% get the index that has offsets, which is also the spot that is activated
% for at least once
% ind23 = find(countmat23 ~= 0);

% get the subscripts of those activated spots
[ioff13, ioff12] = ind2sub(size(countmat23), find(countmat23 ~= 0));

%%% ADD 1 more constraint to check here.
%%% 1. total hit count in the vicinity (3*3-sample rectangle) is reasonable
countmin = 10;   % min counts allowed in the vicinity
countmin = 1;
newflag = zeros(matrange, matrange);
for i = 1: length(ioff13)
    
    totcount = sum(reshape(countmat23(max(ioff13(i)-1,1): min(ioff13(i)+1,matrange), ...
                                      max(ioff12(i)-1,1): min(ioff12(i)+1,matrange)), [], 1));
    if totcount >= countmin
       
        newflag(ioff13(i), ioff12(i)) = 1;

    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% still use ioff13, 12, in order to comment out the following line to
% compare the results before and after applying the new constraints easily

[ioff13, ioff12] = ind2sub(size(newflag), find(newflag ~= 0));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% get the mean & std matrix and vector
STA23count = zeros(length(ioff13), 1);
STA23tbg = [];
STA23tct = [];
STA23tdate = [];
STA23cc = [];

for i = 1: length(ioff13)
    %     i=13;
    count = countmat23(ioff13(i), ioff12(i));
    STA23count(i) = count;
    
    tmp1(1:count,1) = ioff12(i);
    tmp1(1:count,2) = ioff13(i);
    
    tmp2 = tbgmat23(ioff13(i), ioff12(i), 1:count);
    tmp2 = reshape(tmp2, count,1);
    STA23tbg = [STA23tbg; tmp2];
    clear tmp2 ;
    
    tmp2 = tctmat23(ioff13(i), ioff12(i), 1:count);
    tmp2 = reshape(tmp2, count,1);
    STA23tct = [STA23tct; tmp2];
    clear tmp2 ;
    
    tmp2 = ccmat23(ioff13(i), ioff12(i), 1:count);
    tmp2 = reshape(tmp2, count,1);
    STA23cc = [STA23cc; tmp2];
    clear tmp2 ;
    
    tmp2 = tdatemat23(ioff13(i), ioff12(i), 1:count);
    tmp2 = reshape(tmp2, count,1);
    tmp3 = [tmp1 tmp2];
    STA23tdate = [STA23tdate; tmp3];
    clear tmp1 tmp2 tmp3;
    
end

% get the offset 12 & 13 vector
STA13off = (ioff13 -(upmshift+1))/iup;  % divided by iup to recover the original offsets
STA12off = (ioff12 -(upmshift+1))/iup;

%%% filtering correction
if size(filthf,1) > 2   % which means containing the effects of additional stations
    STA12off = STA12off - filthf(2);  % the original correction contains all seven with diff algorithm
    STA13off = STA13off - filthf(3);  % ranges with order
else
    STA12off = STA12off - filthf(1);  % since the total offset 12 contains filtering effect, so substract to correct
    STA13off = STA13off - filthf(2);
end

%%% whether to convert offsets from different fams to a uniform frame
%%% ALGORITHM is confirmed to be true!
[PERMROTS, POLROTS] = GetRotsChao(fam,datapath,0,0);
reftfam = POLROTS(5,4);     % reftime is the reference time at the 1st station,depends on the choice of 1st station
PERMROTS(:,4) = PERMROTS(:,4)-reftfam;    % to make sure that 1st station is 0
POLROTS(:,4) = POLROTS(:,4)-reftfam;

reffam = '043';     %043
[refperm, refpol] = GetRotsChao(reffam,datapath,0,0);
reftref = refpol(5,4);     % reftime is the reference time at the 1st station,depends on the choice of 1st station
refperm(:,4) = refperm(:,4)-reftref;    % to make sure that 1st station is 0
refpol(:,4) = refpol(:,4)-reftref;

% conversion time constants, NOTE that the slow/fast time offset and offset between stas
% is all calculated based on 40 sps
conv(1) = PERMROTS(2,4) - refperm(2,4);   % station 2, LZB
conv(2) = POLROTS(4,4) - refpol(4,4);    % station 3, MGCB
if convflag
    disp('Offsets from different fams are converted to a uniform frame ...');
    spsrat = sps/40;
    STA12off = STA12off-conv(1)*spsrat;   % get the new off based on the ref fam frame
    STA13off = STA13off-conv(2)*spsrat;
end

STA13offsec = STA13off/sps;
STA12offsec = STA12off/sps;


STA23tdate(:, 5) = STA23tdate(:, 3);
STA23tdate(:,1) = (STA23tdate(:,1)-(upmshift+1))/iup; % divided by iup to recover the original offsets
STA23tdate(:,2) = (STA23tdate(:,2)-(upmshift+1))/iup;
if size(filthf,1) > 2
    STA23tdate(:,1) = STA23tdate(:,1)-filthf(2);     % add filtering correction as well
    STA23tdate(:,2) = STA23tdate(:,2)-filthf(3);
else
    STA23tdate(:,1) = STA23tdate(:,1)-filthf(1);     
    STA23tdate(:,2) = STA23tdate(:,2)-filthf(2);
end

if convflag
    STA23tdate(:,1) = STA23tdate(:,1)-conv(1)*spsrat;   % get the new off based on the ref fam frame
    STA23tdate(:,2) = STA23tdate(:,2)-conv(2)*spsrat;
end

STA23tdate(:,3:4)= STA23tdate(:,1:2)/sps;

%%% STA23time structure 9 cols
%%% off12 | off13 | off12sec | off13sec | year-day timing-of-strongest-arrival |
%%% timing-of-center-of-detection-window | cc | famnum
famarr = ceil(str2double(fam))*ones(size(STA23tdate(:,1))); % label of fam number array
STA23time = [STA23tdate STA23tbg STA23tct STA23cc famarr];
%%%%%


% get the mean & std matrix and vector of the additional stations
stdtol = 0.1;   % tolerence of std
countmin = 5;   % min counts allowed in the vicinity
[STA14meansec,STA14stdsec,STA134offsec,STA124offsec,STA134off,STA124off,STA14count,STA14tbg,...
    STA14tct,STA14tdate,STA14cc] = Addcheck_meanstd(1,filthf,countmat14,offmat14,tbgmat14,...
    tctmat14,tdatemat14,ccmat14,upmshift,matrange,sps,stdtol,countmin,convflag,conv,iup);

%%% STA14time has the same structure as STA23time
famarr = ceil(str2double(fam))*ones(size(STA14tct));
STA14time = [STA14tdate STA14tbg STA14tct STA14cc famarr];


[STA15meansec,STA15stdsec,STA135offsec,STA125offsec,STA135off,STA125off,STA15count,STA15tbg,...
    STA15tct,STA15tdate,STA15cc] = Addcheck_meanstd(2,filthf,countmat15,offmat15,tbgmat15,...
    tctmat15,tdatemat15,ccmat15,upmshift,matrange,sps,stdtol,countmin,convflag,conv,iup);
famarr = ceil(str2double(fam))*ones(size(STA15tct));
STA15time = [STA15tdate STA15tbg STA15tct STA15cc famarr];

       
[STA16meansec,STA16stdsec,STA136offsec,STA126offsec,STA136off,STA126off,STA16count,STA16tbg,...
    STA16tct,STA16tdate,STA16cc] = Addcheck_meanstd(3,filthf,countmat16,offmat16,tbgmat16,...
    tctmat16,tdatemat16,ccmat16,upmshift,matrange,sps,stdtol,countmin,convflag,conv,iup);
famarr = ceil(str2double(fam))*ones(size(STA16tct)); 
STA16time = [STA16tdate STA16tbg STA16tct STA16cc famarr];


[STA17meansec,STA17stdsec,STA137offsec,STA127offsec,STA137off,STA127off,STA17count,STA17tbg,...
    STA17tct,STA17tdate,STA17cc] = Addcheck_meanstd(4,filthf,countmat17,offmat17,tbgmat17,...
    tctmat17,tdatemat17,ccmat17,upmshift,matrange,sps,stdtol,countmin,convflag,conv,iup);
famarr = ceil(str2double(fam))*ones(size(STA17tct));
STA17time = [STA17tdate STA17tbg STA17tct STA17cc famarr];


%% PLOT
%%% the final integrated map mean & standard deviation of offset for all days
if convflag == 0    % plot only for individual fam, not the merged
    figure('Position',[scrsz(3)/10 scrsz(4)/10 3*scrsz(3)/5 4.5*scrsz(4)/5]);
    
    nrow = nstanew+1;
    ncol = 3;
    msize = 1.5;
    
    meancran = [-1 1];
    stdcran = [-0.5 0.5];
    hitcran = [0 5];
    tothitcran = [0 10];
    axran = [-1 1 -1 1];
    
    % additional sta 1, mean offset
    ax1=subplot(nrow, ncol, 1,'align');
    scatter(STA124offsec, STA134offsec, msize, STA14meansec, 'filled','s');    %, 'MarkerEdgeColor', 'w')
    colormap(ax1,jet)
    colorbar
    caxis(meancran);
    box on
    axis equal
    xlabel('STA12 offset (s)')
    ylabel('STA13 offset (s)')
    axis(axran);
    title(strcat(num2str(lo),'-',num2str(hi),'\_mean\_', strtrim(stasnew(1, :))));
    
    % additional sta 1, offset std
    ax2=subplot(nrow, ncol, 2,'align');
    scatter(STA124offsec, STA134offsec, msize, STA14stdsec, 'filled','s');    %, 'MarkerEdgeColor', 'w')
    colormap(ax2,jet)
    colorbar
    caxis(stdcran);
    box on
    axis equal
    xlabel('STA12 offset (s)')
    ylabel('STA13 offset (s)')
    axis(axran);
    title(strcat(num2str(lo),'-',num2str(hi),'\_std\_', strtrim(stasnew(1, :))));
    
    % additional sta 1, hit count
    ax3=subplot(nrow, ncol, 3,'align');
    scatter(STA124offsec, STA134offsec, msize, STA14count, 'filled','s');    %, 'MarkerEdgeColor', 'w')
    oldcmap = colormap(ax3, hot);
    colormap(ax3, flipud(oldcmap) );
    colorbar
    caxis(hitcran);
    box on
    axis equal
    xlabel('STA12 offset (s)')
    ylabel('STA13 offset (s)')
    axis(axran);
    title(strcat(num2str(lo),'-',num2str(hi),'\_counts\_', strtrim(stasnew(1, :))));
    
    % additional sta 2, mean offset
    ax4=subplot(nrow, ncol, 4,'align');
    scatter(STA125offsec, STA135offsec, msize, STA15meansec, 'filled','s');    %, 'MarkerEdgeColor', 'w')
    colormap(ax4,jet)
    colorbar
    caxis(meancran);
    box on
    axis equal
    xlabel('STA12 offset (s)')
    ylabel('STA13 offset (s)')
    axis(axran);
    title(strcat(num2str(lo),'-',num2str(hi),'\_mean\_', strtrim(stasnew(2, :))));
    
    % additional sta 2, offset std
    ax5=subplot(nrow, ncol, 5,'align');
    scatter(STA125offsec, STA135offsec, msize, STA15stdsec, 'filled','s');    %, 'MarkerEdgeColor', 'w')
    colormap(ax5,jet)
    colorbar
    caxis(stdcran);
    box on
    axis equal
    xlabel('STA12 offset (s)')
    ylabel('STA13 offset (s)')
    axis(axran);
    title(strcat(num2str(lo),'-',num2str(hi),'\_std\_', strtrim(stasnew(2, :))));
    
    % additional sta 2, hit count
    ax6=subplot(nrow, ncol, 6,'align');
    scatter(STA125offsec, STA135offsec, msize, STA15count, 'filled','s');    %, 'MarkerEdgeColor', 'w')
    oldcmap = colormap(ax6, hot);
    colormap(ax6, flipud(oldcmap) );
    colorbar
    caxis(hitcran);
    box on
    axis equal
    xlabel('STA12 offset (s)')
    ylabel('STA13 offset (s)')
    axis(axran);
    title(strcat(num2str(lo),'-',num2str(hi),'\_counts\_', strtrim(stasnew(2, :))));
    
    % additional sta 3, mean offset
    ax7=subplot(nrow, ncol, 7,'align');
    scatter(STA126offsec, STA136offsec, msize, STA16meansec, 'filled','s');    %, 'MarkerEdgeColor', 'w')
    colormap(ax7,jet)
    colorbar
    caxis(meancran);
    box on
    axis equal
    xlabel('STA12 offset (s)')
    ylabel('STA13 offset (s)')
    axis(axran);
    title(strcat(num2str(lo),'-',num2str(hi),'\_mean\_', strtrim(stasnew(3, :))));
    
    % additional sta 3, offset std
    ax8=subplot(nrow, ncol, 8,'align');
    scatter(STA126offsec, STA136offsec, msize, STA16stdsec, 'filled','s');    %, 'MarkerEdgeColor', 'w')
    colormap(ax8,jet)
    colorbar
    caxis(stdcran);
    box on
    axis equal
    xlabel('STA12 offset (s)')
    ylabel('STA13 offset (s)')
    axis(axran);
    title(strcat(num2str(lo),'-',num2str(hi),'\_std\_', strtrim(stasnew(3, :))));
    
    % additional sta 3, hit count
    ax9=subplot(nrow, ncol, 9,'align');
    scatter(STA126offsec, STA136offsec, msize, STA16count, 'filled','s');    %, 'MarkerEdgeColor', 'w')
    oldcmap = colormap(ax9, hot);
    colormap(ax9, flipud(oldcmap) );
    colorbar
    caxis(hitcran);
    box on
    axis equal
    xlabel('STA12 offset (s)')
    ylabel('STA13 offset (s)')
    axis(axran);
    title(strcat(num2str(lo),'-',num2str(hi),'\_counts\_', strtrim(stasnew(3, :))));
    
    % additional sta 4, mean offset
    ax10=subplot(nrow, ncol, 10,'align');
    scatter(STA127offsec, STA137offsec, msize, STA17meansec, 'filled','s');    %, 'MarkerEdgeColor', 'w')
    colormap(ax10,jet)
    colorbar
    caxis(meancran);
    box on
    axis equal
    xlabel('STA12 offset (s)')
    ylabel('STA13 offset (s)')
    axis(axran);
    title(strcat(num2str(lo),'-',num2str(hi),'\_mean\_', strtrim(stasnew(4, :))));
    
    % additional sta 4, offset std
    ax11=subplot(nrow, ncol, 11,'align');
    scatter(STA127offsec, STA137offsec, msize, STA17stdsec, 'filled','s');    %, 'MarkerEdgeColor', 'w')
    colormap(ax11,jet)
    colorbar
    caxis(stdcran);
    box on
    axis equal
    xlabel('STA12 offset (s)')
    ylabel('STA13 offset (s)')
    axis(axran);
    title(strcat(num2str(lo),'-',num2str(hi),'\_std\_', strtrim(stasnew(4, :))));
    
    % additional sta 4, hit count
    ax12=subplot(nrow, ncol, 12,'align');
    scatter(STA127offsec, STA137offsec, msize, STA17count, 'filled','s');    %, 'MarkerEdgeColor', 'w')
    oldcmap = colormap(ax12, hot);
    colormap(ax12, flipud(oldcmap) );
    colorbar
    caxis(hitcran);
    box on
    axis equal
    xlabel('STA12 offset (s)')
    ylabel('STA13 offset (s)')
    axis(axran);
    title(strcat(num2str(lo),'-',num2str(hi),'\_counts\_', strtrim(stasnew(4, :))));
    
    % original trio hit count
    ax13=subplot(nrow, ncol, 13,'align');
    scatter(STA12offsec, STA13offsec, msize, STA23count, 'filled','s');    %, 'MarkerEdgeColor', 'w')
    oldcmap = colormap(ax13, hot);
    colormap(ax13, flipud(oldcmap) );
    colorbar
    box on
    axis equal
    caxis(tothitcran);
    xlabel('STA12 offset (s)')
    ylabel('STA13 offset (s)')
    axis(axran);
    title(strcat(num2str(lo),'-',num2str(hi),'\_trio-counts\_'));
end

%% save the image & the matfile
if convflag == 0 
    if iup == 1
        PREFIX = strcat(fam,'.lo',num2str(loopoffmax),'.ccm',num2str(xcmaxAVEnmin),'.',int2str(npo),...
                        int2str(npa),'.ms', int2str(mshift));
    else
        PREFIX = strcat(fam,'.up.lo',num2str(loopoffmax),'.ccm',num2str(xcmaxAVEnmin),'.',int2str(npo),...
                        int2str(npa),'.ms', int2str(mshift));
    end                
    print('-depsc',strcat(rstpath,'/FIGS/',PREFIX,'_','AddOffCheck', '_',num2str(lo),'-',...
          num2str(hi),'_',strtrim(stasnew(1, :)),'&',strtrim(stasnew(2, :)),'&',...
          strtrim(stasnew(3, :)),'&',strtrim(stasnew(4, :)),'.eps'));
else
    if iup == 1
        PREFIX = strcat(fam,'.conv.lo',num2str(loopoffmax),'.ccm',num2str(xcmaxAVEnmin),'.',...
                        int2str(npo),int2str(npa),'.ms', int2str(mshift));
    else
        PREFIX = strcat(fam,'.upconv.lo',num2str(loopoffmax),'.ccm',num2str(xcmaxAVEnmin),'.',...
                        int2str(npo),int2str(npa),'.ms', int2str(mshift));
    end
end

famarr = ceil(str2double(fam))*ones(size(STA12off));   % label of fam number array
if ~isempty(famarr)
    savefile23 = [STA12off STA13off STA12offsec STA13offsec STA23count famarr]; % 6 cols
    fid = fopen(strcat(rstpath, '/MAPS/countori_',PREFIX,'_',num2str(lo),'-',num2str(hi),'_',...
        num2str(winlen/sps),'s',num2str(sps),'sps', num2str(nstanew), 'nsta'),'w+');
    fprintf(fid,'%d %d %.4f %.4f %d %d \n',savefile23');
    fclose(fid);
    fid = fopen(strcat(rstpath, '/MAPS/timeori_',PREFIX,'_',num2str(lo),'-',num2str(hi),'_',...
        num2str(winlen/sps),'s',num2str(sps),'sps', num2str(nstanew), 'nsta'),'w+');
    fprintf(fid,'%d %d %.4f %.4f %d %.4f %.1f %.4f %d \n',STA23time');  % 9 cols
    fclose(fid);
end


famarr = ceil(str2double(fam))*ones(size(STA124off));
if ~isempty(famarr)
    savefile14 = [STA124off STA134off STA124offsec STA134offsec STA14meansec STA14stdsec STA14count ...
        famarr];
    fid = fopen(strcat(rstpath, '/MAPS/countadd_', stasnew(1,:), '_',PREFIX,'_',num2str(lo),'-',...
        num2str(hi),'_',num2str(winlen/sps),'s',num2str(sps),'sps', num2str(nstanew), 'nsta'),...
        'w+');
    fprintf(fid,'%d %d %.4f %.4f %.4f %.4f %d %d \n',savefile14');  % 8 cols
    fclose(fid);
    fid = fopen(strcat(rstpath, '/MAPS/timeadd_', stasnew(1,:), '_',PREFIX,'_',num2str(lo),'-',...
        num2str(hi),'_',num2str(winlen/sps),'s',num2str(sps),'sps', num2str(nstanew), 'nsta'),...
        'w+');
    fprintf(fid,'%d %d %.4f %.4f %d %.4f %.1f %.4f %d \n',STA14time');  % 9 cols
    fclose(fid);
end

famarr = ceil(str2double(fam))*ones(size(STA125off));
if ~isempty(famarr)
    savefile15 = [STA125off STA135off STA125offsec STA135offsec STA15meansec STA15stdsec STA15count ...
        famarr];
    fid = fopen(strcat(rstpath, '/MAPS/countadd_', stasnew(2,:), '_',PREFIX,'_',num2str(lo),'-',...
        num2str(hi),'_',num2str(winlen/sps),'s',num2str(sps),'sps', num2str(nstanew), 'nsta'),...
        'w+');
    fprintf(fid,'%d %d %.4f %.4f %.4f %.4f %d %d \n',savefile15');
    fclose(fid);
    fid = fopen(strcat(rstpath, '/MAPS/timeadd_', stasnew(2,:), '_',PREFIX,'_',num2str(lo),'-',...
        num2str(hi),'_',num2str(winlen/sps),'s',num2str(sps),'sps', num2str(nstanew), 'nsta'),...
        'w+');
    fprintf(fid,'%d %d %.4f %.4f %d %.4f %.1f %.4f %d \n',STA15time');
    fclose(fid);
end

famarr = ceil(str2double(fam))*ones(size(STA126off));
if ~isempty(famarr)
    savefile16 = [STA126off STA136off STA126offsec STA136offsec STA16meansec STA16stdsec STA16count ...
        famarr];
    fid = fopen(strcat(rstpath, '/MAPS/countadd_', stasnew(3,:), '_',PREFIX,'_',num2str(lo),'-',...
        num2str(hi),'_',num2str(winlen/sps),'s',num2str(sps),'sps', num2str(nstanew), 'nsta'),...
        'w+');
    fprintf(fid,'%d %d %.4f %.4f %.4f %.4f %d %d \n',savefile16');
    fclose(fid);
    fid = fopen(strcat(rstpath, '/MAPS/timeadd_', stasnew(3,:), '_',PREFIX,'_',num2str(lo),'-',...
        num2str(hi),'_',num2str(winlen/sps),'s',num2str(sps),'sps', num2str(nstanew), 'nsta'),...
        'w+');
    fprintf(fid,'%d %d %.4f %.4f %d %.4f %.1f %.4f %d \n',STA16time');
    fclose(fid);
end

famarr = ceil(str2double(fam))*ones(size(STA127off));
if ~isempty(famarr)
    savefile17 = [STA127off STA137off STA127offsec STA137offsec STA17meansec STA17stdsec STA17count ...
        famarr];
    fid = fopen(strcat(rstpath, '/MAPS/countadd_', stasnew(4,:), '_',PREFIX,'_',num2str(lo),'-',...
        num2str(hi),'_',num2str(winlen/sps),'s',num2str(sps),'sps', num2str(nstanew), 'nsta'),...
        'w+');
    fprintf(fid,'%d %d %.4f %.4f %.4f %.4f %d %d \n',savefile17');
    fclose(fid);
    fid = fopen(strcat(rstpath, '/MAPS/timeadd_', stasnew(4,:), '_',PREFIX,'_',num2str(lo),'-',...
        num2str(hi),'_',num2str(winlen/sps),'s',num2str(sps),'sps', num2str(nstanew), 'nsta'),...
        'w+');
    fprintf(fid,'%d %d %.4f %.4f %d %.4f %.1f %.4f %d \n',STA17time');
    fclose(fid);
end

end   % loop end for fam

% keyboard









