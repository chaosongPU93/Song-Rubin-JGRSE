%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% THIS is the low-frequency detection with LZB-TWKB-MGCB trio in Cascadia
% Using data at 40 sps, but down-sampled to 20 sps, then add extra stations
% to check the validity of detections, version 2 (Currently used)
%
% USE the same maximum offset in samples between stations in order to get
% the same detected area, considering the different sampling rate ofc
%
% Features:
%   1. Low freq. naturally has low resolution, so 40 sps is redundant,
%       downsample to 20 sps.
%   2. Down-sample has a similar to increase the mshift allowed to resolve
%       a large area, but significantly save the computation cost
%   3. Allow an extra mshift, but if final location is inside the extra
%        range, reject it in case cycle interference
%   4. To exmine the lf detections, add new stations to do xcorr
%   5. Avoid double-counting detections
%   6. Modify the code to save or plot all results after implementing step 5
%   7. Adding the option to use data without station response
%   8. modify the codes to use less memory by avoiding saving cumsum
%   matrix and do more computation instead.
%   9. STICK on PGC N*1
%   10. Use the days the same as in Peng et al. 2015 instead of different days
%       that bond to different fams;
%
% Difference between hf_LZBtrio_detection_check:
%   1. Previously, i am using different days to stack for lfe templates as well as
%       detection that are active inside that fam (set a min threshold of lfe
%       detections in that day based on bostock's catalog); It is reasonable in 
%       obtaining stacks or rotation parameters, however could cause problems in 
%       detection. For example, there might be detections in day 1 fam 1 and fam 2,
%       but if do detection in day 1 only in fam 1, the merging results would have
%       the bias that the complete pattern would not be resolved. So now use the same
%       days as that in Peng et al. 2015
%   2. Comment lines that seems inrelevant to current analysis, no wiggles plot as well
%
% NOTES:
%   16/8/4/0.4/19; 16/8/4/0.35/14;
%   fam 010, amp of TWKB is very small
%
% Similar to lf_PGCtrio_detection_resp_check_v2.m
%        
% 
%
% Modified by Chao Song, chaosong@princeton.edu
% First created date:   2021/02/24
% Last modified date:   2021/02/24
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Initialization
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
% add some families at the southwestern part of our study region
nfampool = ['001';
            '019';
            '045'];
        
% the fams in the SW part of the currect study region
nfampool = ['158';
            '234'];
        
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
            '006';
            '001';
%             '158';      % adding 158 mainly to test
            ];
        
        
for ifam = 1: size(nfampool,1)
fam = nfampool(ifam, :)
% fam='043';     % family number

freqflag='lf';  % flag to indicate whether to do hf or lf;

% get days, 3-sta trio name, bostock's template name, zero crossing index
% of templates
timoffrot= [2003 060;
            2003 061;
            2003 062; % dates of the year, two columns: year Julian_day
            2003 063;
            2003 064;
%             2003 065;   % LZB only has data from 2003 060-064
%             2003 066;
%             2003 067;
            2004 194;
            2004 195;
            2004 196;
            2004 197;
            2004 198;
            2004 199;
            2004 200;
            2004 201;
            2004 202;
            2004 203;
            2005 254;
            2005 255;
            2005 256;
            2005 257;
            2005 258;
            2005 259;
            2005 260;
            2005 261];

nday = size(timoffrot, 1);

% the following is gained by looking at the templates with eyes 
if isequal(fam,'002')
    tempoffs = [612 612 612];
elseif isequal(fam,'043')
    tempoffs = [599 599 599];
elseif isequal(fam,'141')
    tempoffs = [592 592 592];
elseif isequal(fam,'047')
    tempoffs = [610 610 610];
elseif isequal(fam,'010')
    tempoffs = [589 589 589];
elseif isequal(fam,'144')
    tempoffs = [597 597 597];
elseif isequal(fam,'099')
    tempoffs = [591 591 591];
elseif isequal(fam,'068')
    tempoffs = [598 598 598];
elseif isequal(fam,'125')
    tempoffs = [601 601 601];
elseif isequal(fam,'147')
    tempoffs = [597 597 597];
elseif isequal(fam,'017')
    tempoffs = [595 595 595];
    elseif isequal(fam,'001')
    tempoffs = [586 586 586];
elseif isequal(fam,'019')
    tempoffs = [597 597 597];    
elseif isequal(fam,'045')
    tempoffs = [590 590 590];
else
    tempoffs = [597 597 597];   
end

% get permanent and polaris station rotation parameters
sft2=0;     % centroid shift of station 2
sft3=0;     % centroid shift of station 3
[PERMROTS, POLROTS] = GetRotsChao(fam,datapath,sft2,sft3);
reftime = POLROTS(5,4);     % reftime is the reference time at the 1st station,depends on the choice of 1st station
PERMROTS(:,4) = PERMROTS(:,4)-reftime;    % to make sure that 1st station is 0
POLROTS(:,4) = POLROTS(:,4)-reftime;

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

nsta=size(stas,1);         %  number of stations
PERMROTS(:,2:3)=pi*PERMROTS(:,2:3)/180.;     % convert 2-3 columns to rad 
POLROTS(:,2:3)=pi*POLROTS(:,2:3)/180.;
sps=20;     % samples per second

tempwinsec = 60;
tempwinlen = tempwinsec*sps;    % template window length, 60s, 60*sps samples, 1 min

%%% IMPORTANT, NEED to change sometimes
%Basics of the cross-correlation:  Window length, number of windows, filter parameters, etc.
% winlensec=20;     % 20 was used in first-year report and so was 10 winoffsec 
winlensec=16;     % cc window length in sec
winoffsec=8;        % window offset in sec, which is the step of a moving window
winlen=winlensec*sps;      % length in smaples
winoff=winoffsec*sps;      % offset in samples
tracelen=86400*sps; %one day of data at 40 sps, overall trace length, 24*3600
%%% winbig could be problematic if mshift+1 >= 2*sps
cutsec = 2;
winbig=2*(tracelen/2-(cutsec*sps)); %ignore 2 seconds at each end of day, bigger window contains positive and negative, why tracelen/2?, -4s
timbig=winbig/(2*sps); %half that time, in seconds, half day - 2s
igstart=floor(tracelen/2-winbig/2)+1; %start counting seis data from here, 2*sps+1, floor(1.8)==1
nwin=floor((winbig-winlen)/winoff)+1;    % number of windows, first one not included, ADD +1 by Chao, 2019/02/17

%UPGRADING SINCE MODIFYING READPOLS & READPERMS STOPED HERE
%hi=6.5;  %002 Stanford
%lo=1.25; %002 Stanford
%hi=6;
%lo=1.5;
hi=1.25;
lo=0.5;
npo=2;     % poles, passes of filters
npa=2;
cyclskip = 20;
mshift=14+cyclskip; %19; %maximum shift for the x-correlations. 19 for 002 Stanford,    % in sps, 0.5s*40sps=20
mshift=7;
loopoffmax=4; %1.5 for standard 1.5-6Hz; 4 for 0.5-1.5Hz.  2 for non-interpolated.   % what is loopoffmax, the circuit of time offsets
% xcmaxAVEnmin=0.5 was used in first-year report
xcmaxAVEnmin=0.35; %0.44; %0.44 for 002 Stanford %0.45; %0.36 for 4s 1-12 Hz; %0.4 for 4s 1.5-6 Hz and 6s 0.5-1.5Hz; 0.36 for 4s2-8 Hz ; 0.38 for 4s0.75-6 Hz; 0.094 for 128s 2-8Hz;  0.1 for 128s 1.5-6Hz; 0.44 for 3-s window?
%%% xcmaxAVEnmin = x-correlation max average min ?

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Get the templates for the stations. Won't be always required.  PGCopt_002_1.25-6.5Hz_2pass_100sps.le14shift.resamp
% %%% THIS loads the templates from Allan
% for ista=1:nsta
%     if stas(ista,4) == ' '          % which means it has a 3-letter code, permanent stations
%         temptemps(ista,:)=load([stas(ista,1:3),'opt_',fam,'_',num2str(lo),'-',num2str(hi), ...
%                          'Hz_',int2str(npa),'pass_',int2str(sps),'sps.input'],'w');         % write permission to data file
%     else          % which means it has a 4 (maybe more)-letter code, Polaris stations
%         temptemps(ista,:)=load([stas(ista,1:4),'opt_',fam,'_',num2str(lo),'-',num2str(hi), ...
%                          'Hz_',int2str(npa),'pass_',int2str(sps),'sps.input'],'w');
%     end
% end
% %%%

% %%% THIS loads the  templates from Chao, this is the new template generated by Make_templatesv2.m
% %%% using the same rotation para from calc_rots_noresp.m, different from the PGC trio. But pay
% %%% attention here, templates DO NOT affect the detection itself!!!
% for ista=1:nsta
%     fname = strcat(temppath, fam, '_', strtrim(stas(ista, :)), '_', num2str(40), 'sps_', ...
%                    num2str(tempwinsec), 's_', 'BBCCS_', 'opt_Nof_Non_Chao');
%     data = load(fname);
%     [num, denom] = rat(sps/40);
% %     temp1(ista,:) = Bandpass(data, sps, lo, hi, npo, npa, 'butter');    % bandpass with wrong sps
% %     temp2(ista,:) = resample(temp1(ista,:), num, denom);    % continue to downsample
%     %%% ATTENTION: the sample rate in bandpass must be the same as the data 
%     data = Bandpass(data, 40, lo, hi, npo, npa, 'butter');
%     temptemps(ista,:) = resample(data, num, denom); 
%     
% end
% % %%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% %%% a simple plot to find the temoffs
% figure(111)
% for ista=1:nsta
%     plot(temptemps(ista,:)); hold on
% end
% %%%
% 
% figure
% for ista = 1: nsta
%     %%% normalization
%     ampmax = max(temptemps(ista, :));
%     ampmin = min(temptemps(ista, :));
%     norm = max(ampmax, -ampmin);
%     tempno(ista, :) = temptemps(ista, :)/ norm;
%     plot(tempno(ista,:)); hold on
% end
% vertical_cursors
% 
% 
% % tempbef=70;
% % tempaft=89;
% tempbef=39;        % template before time in samples
% tempaft=40;         % template after time in samples
% templen=tempbef+tempaft+1;        % template length in sec
% for ista=1:nsta
%     STAtemps(ista,:)=temptemps(ista,tempoffs(ista)-tempbef:tempoffs(ista)+tempaft);   % 120 samples, 3s
%     snips(templen*(ista-1)+1: templen*ista)=STAtemps(ista,:);       % concatenate all data to one vector
% end

%%% 'scalefact' scales templates; 'scaleseisms' scales seisms.  Strategy changes with family.
if isequal(fam,'002')
    tempoffs=tempoffs-1; %Center of strongest window is 1 or 2 samples ahead of zero crossing (002); make it 1.
    scaleseisms=[1.0 0.76 0.95];       % scaleseisms scales seismograms
elseif isequal(fam,'043')
    tempoffs=tempoffs-1;
    scaleseisms=[1.0 1.0 1.0];
else
    tempoffs=tempoffs-1;
    scaleseisms=[1.0 1.0 1.0];
end
% minses=-min(STAtemps,[],2); %STAtemps is (currently) 3 by 120     % give the min/max of the each trace(row) of the STAtemps
% maxses= max(STAtemps,[],2);
% plustominus=maxses./minses;        % ratio of +/-
% scalefact=minses*max(plustominus); %This is used to scale templates, just for plotting purposes
% for ista=1:nsta
%     STAtemps(ista,:)=STAtemps(ista,:)/scalefact(ista); %This plots the templates with the largest positive value (of any) at +1
% end
%%%%%%%%%%% EXAMPLE of scaling:
% STAtemps = [5 -4 1; 10 -5 2];
% minses=-min(STAtemps,[],2);
% maxses= max(STAtemps,[],2);
% plustominus=maxses./minses;
% scalefact=minses*max(plustominus);
% for ista=1:2
%     STAtemps(ista,:)=STAtemps(ista,:)/scalefact(ista);
% end
% STAtemps =
% 
%     0.6250   -0.5000    0.1250
%     1.0000   -0.5000    0.2000
%%%%%%%%%%%
% figure
% plot(STAtemps(1,:),'r', 'linewidth', 1.5)
% hold on
% plot(STAtemps(2,:),'b', 'linewidth', 1.5)
% plot(STAtemps(3,:),'k', 'linewidth', 1.5)
% drawnow     % update figures immediately

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% START TO LOOP FOR every day
% get the all catalog LFEs in that family
bostname = ('/BOSTOCK/total_mag_detect_0000_cull_NEW.txt');
catalog = load(strcat(workpath, bostname));
famnum = str2double(fam);
dateall = catalog(famnum == catalog(:, 1), :);

% cycle over each day:
for nd=1:nday      % num of rows, also num of days
%     nd=1;
tic
    close all
    %Bostock's Detections:
    year = timoffrot(nd,1);
    YEAR = int2str(year);
    jday = timoffrot(nd,2);
    if year == 2003
        date = jday-62+30303;
    elseif year == 2004
        date = jday-196+40714;
    elseif year == 2005
        date = jday-254+50911;
    end
    bostocks = dateall(date == dateall(:, 2), :);
    bostsec = 3600*(bostocks(:,3)-1)+bostocks(:,4);
    bostsamp = round(bostsec*sps);
    
    if jday <= 9
        JDAY=['00',int2str(jday)];
    elseif jday<= 99
        JDAY=['0',int2str(jday)];
    else
        JDAY=int2str(jday);
    end
%     MO = 'SEP';
    MO=day2month(jday,year);     % EXTERNAL function, day2month, get the month of one particular date
    IDENTIF=[YEAR,'.',JDAY,'.',fam,'.lo',num2str(loopoffmax),'.cc',num2str(xcmaxAVEnmin),'.',...
             int2str(npo),int2str(npa),'.ms', int2str(mshift)]
    direc=[datapath, '/', YEAR,'/',MO,'/'];     % directory name
%     fprintf('%s \n',direc);
    prename=[direc,YEAR,'.',JDAY,'.00.00.00.0000.CN'];    %  path plus prefix of data file,
    
    if year==2003   % in 2003, station KLNB is named KELB, so use KELB to replace KLNB in 2003
        POLSTA(3,:)='KELB ';
    else
        POLSTA(3,:)='KLNB ';  % remember to change it back
    end
    
    %Read the data; find glitches (many consecutive zeros) and flag those windows in STAnzeros.
    %Get timsSTA from the permanent stations (last one over-writes):
    STAopt=zeros(nsta,tracelen);    % one day of samples
    STAort=zeros(nsta,tracelen);
    STAnzeros=zeros(nsta,nwin);    % used to flag windows with glitches
    
    fileflag = 1;   % 1 means all files exist, the file status is normal
    for ista=1:nsta
        found=0;
        
        [LIA,idx]=ismember(stas(ista,:),PERMSTA,'rows');     % ismember, to determine whether each row of stas is contained in PERMSTA, return logical value 1/0 and index
        if LIA
            found=found+LIA;
            if strcmp(PERMSTA(idx,1:3),'PGC')     % string compare
                %%% WAHT is the meanning of 'fact', similar to instrument
                %%% response
                fact=1.0e-3;
            elseif strcmp(PERMSTA(idx,1:3),'LZB')
                fact=1.8e-3;
            end
            %%% readperms is an EXTERNAL FUNCTION
            % opt: optimal seismogram after rotations
            % ort: orthogonal seismogram after rotations
            % nzeros: number of zeros in the trace
            % timsSTA: time sequence
            fname = strcat(prename,'.',PERMSTA(idx,1:3),'..BHE.D.SAC');
            if isfile(fname)    % if have the data file
%                 % 1. this is for data without removing station response
%                 [opt,ort,nzeros,timsSTA]=readperms(prename,PERMSTA,PERMROTS,idx,sps,lo,hi,npo,npa,fact,nwin,winlen,winoff,igstart);
                % 2. this is for data with no response
                [opt,ort,nzeros,timsSTA]=readpermsv2(prename,PERMSTA,PERMROTS,idx,sps,lo,hi,npo,...
                                                     npa,fact,nwin,winlen,winoff,igstart);                
            else
                fileflag = 0;   % change the file flag to 0, meaning abnormal
                fprintf('No data for station %s in day %s %s, this day will be omitted. \n',...
                        PERMSTA(idx,1:3), YEAR, JDAY);
                break   % break the entire station loop
            end                
        end
        
        [LIA,idx]=ismember(stas(ista,:),POLSTA,'rows');    
        if LIA        % if are in POLSTA
            found=found+LIA; %better be 1
            if year==2003 && jday<213        % should result from some criteria
                fact=7.5e-3;
            else
                fact=1.5e-3; 
            end
            
            %%% readpols is an EXTERNAL FUNCTION
            % opt:
            % ort:
            % nzeros:
            fname = strcat(prename,'.',POLSTA(idx,1:4),'..HHE.D.SAC');
            if isfile(fname)    % if have the data file
%                 % 1. this is for data without removing station response
%                 [opt,ort,nzeros]=readpols(prename,POLSTA,POLROTS,idx,sps,lo,hi,npo,npa,fact,nwin,winlen,winoff,igstart);
                % 2. this is for data with no response
                [opt,ort,nzeros]=readpolsv2(prename,POLSTA,POLROTS,idx,sps,lo,hi,npo,npa,fact,...
                                            nwin,winlen,winoff,igstart);
            else
                fileflag = 0;   % change the file flag to 0, meaning abnormal
                fprintf('No data for station %s in day %s / %s. \n',POLSTA(idx,1:4), YEAR, JDAY);
                break   % break the entire station loop
            end                  
        end
%         found=found         % could be a benchmark
        %factr1(ista)=prctile(abs(opt),90); %Not normalized
        %factr2(ista)=factr1(ista)/factr1(1); %This is what is used; keeps 1st station unchanged but scales the others
        STAopt(ista,:)=opt/scaleseisms(ista);   % what is STAopt for??
        STAort(ista,:)=ort/scaleseisms(ista);
        STAnzeros(ista,:)=nzeros;    % STAnzeros is to count the number of zeros in each window of each station
    end
    
    if fileflag == 0    % means there are missing files
        fprintf('Day %s / %s will be omitted because of missing files. \n', YEAR, JDAY);
        continue    % continue to the next day
    end
        
    
    %NOW for broader band (bb)
    %%% the broader frequency band is the main difference compared to above
    lobb=1.25;
    hibb=6.5;
%    lobb=4;
%    hibb=10;
    STAoptbb=zeros(nsta,tracelen);
    for ista=1:nsta
        [LIA,idx]=ismember(stas(ista,:),PERMSTA,'rows');
        if LIA
            if strcmp(PERMSTA(idx,1:3),'PGC')
                fact=1.0e-3;
            elseif strcmp(PERMSTA(idx,1:3),'LZB')
                fact=1.8e-3;
            end
%             [opt,~,~,~]=readperms(prename,PERMSTA,PERMROTS,idx,sps,lobb,hibb,npo,npa,fact,nwin,winlen,winoff,igstart);
            [opt,~,~,~]=readpermsv2(prename,PERMSTA,PERMROTS,idx,sps,lobb,hibb,npo,npa,fact,...
                                    nwin,winlen,winoff,igstart);
        end
        [LIA,idx]=ismember(stas(ista,:),POLSTA,'rows');
        if LIA
            if year==2003 && jday<213
                fact=7.5e-3;
            else
                fact=1.5e-3; 
            end          
%             [opt,~,~]=readpols(prename,POLSTA,POLROTS,idx,sps,lobb,hibb,npo,npa,fact,nwin,winlen,winoff,igstart);
            [opt,~,~]=readpolsv2(prename,POLSTA,POLROTS,idx,sps,lobb,hibb,npo,npa,fact,nwin,...
                                 winlen,winoff,igstart);
        end
        STAoptbb(ista,:)=opt/scaleseisms(ista); 
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %   Autocorrelation of stations.  Those that end in "sq" are the running
    %   cumulative sum, to be used later by differencing the window edpoints.
    %   (Used to be PGCauto, PGC2, SSIBauto, SSIB2, etc.)
    %   Station to itself is in a 3 x tracelen array
    %   Cross-station measurements are in their own linear array
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    STAauto=STAopt.^2;
%     STAsq=cumsum(STAauto,2);    % self cumulative sum of each trace(row)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %  Cross-correlation between stations, with small offsets up to +/- mshift.
    %  First index is pointwise multiplication of traces; second is shifting offset.
    %  lenx is shorter than tracelen by mshift at each end (see notebook sketch)
    %  For STA12 and PGSI, SSI and SIL are shifted relative to PGC, by 1 each time through loop.
    %  For SISS, SSI is shifted relative to SILB.
    %  cumsumSTA12 etc. are the running cumulative sum of the x-correlation.
    %  PGSSx becomes STA12x, PGSI -> STA13, SISS -> STA32
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    lenx=tracelen-2*mshift;     % 86400*sps-19*2
    STA12x=zeros(lenx, 2*mshift+1);    % 19*2+1
    STA13x=zeros(lenx, 2*mshift+1);    % stas: 1->PGC, 2->SSIB, 3->SILB
    STA32x=zeros(lenx, 2*mshift+1);
    %%% SEE NOTES, Chao
    for n=-mshift:mshift
        % PGC corr SSIB, 1+mshift:tracelen-mshift == 1+mshift-n:tracelen-mshift-n == lenx
        STA12x(:,n+mshift+1)=STAopt(1,1+mshift:tracelen-mshift).* ...
            STAopt(2,1+mshift-n:tracelen-mshift-n);
        % PGC corr SILB
        STA13x(:,n+mshift+1)=STAopt(1,1+mshift:tracelen-mshift).* ...
            STAopt(3,1+mshift-n:tracelen-mshift-n);
        % SILB corr SSIB
        STA32x(:,n+mshift+1)=STAopt(3,1+mshift:tracelen-mshift).* ...
            STAopt(2,1+mshift-n:tracelen-mshift-n);
    end
%     cumsumSTA12=cumsum(STA12x);  %prev cumsumPGSS    % ==cumsum(STA12x,1),column
%     cumsumSTA13=cumsum(STA13x);  %prev cumsumPGSI
%     cumsumSTA32=cumsum(STA32x);  %prev cumsumSISS

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %  "winbig"             is now the whole day, minus 2 sec at each end (apparently).
    %  "timbig"             is the time of half that.
    %  'igstart'            is the index of the starting sample.
    %  'winlen and winoff'  refer to the small windows.
    %  'timswin'            refers to the central times of those small windows.
    %  'sumsPGSS' (etc.)    is the cross-correlation sum over the window.  The first
    %    index refers to the window number and the second the shift over +/-mshift.
    %
    %  Normalized x-correlation:
    %    For PGSS and PGSI, for a given window PGC does not shift but SSI and 
    %    SIL do.  So can compute sumsPGC2 (from the running cum. sum PGC2) just
    %    once for each window.  Same for sumsSILB2b.  But for the stations that
    %    shift, SSI and SIL (for PGC) and SSI (for SIL), must compute sumsSSIB2 
    %    and sumsSILB2 upon each shift (actually, this is is easy book-keeping
    %    but not efficient).  Again, the first index refers to the window
    %    number and the second the shift over +/-mshift.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    timswin=zeros(nwin,1);    %%% timswin refers to the central times of those small windows.
    sumsSTA12=zeros(nwin,2*mshift+1);   % PGC-SSIB
    sumsSTA13=zeros(nwin,2*mshift+1);   % PGC-SILB
    sumsSTA32=zeros(nwin,2*mshift+1);   % SILB-SSIB
    sumsSTA1sq=zeros(nwin,2*mshift+1);  % "sq" are the running cumulative sum, to be used later by differencing the window edpoints
    sumsSTA2sq=zeros(nwin,2*mshift+1);
    sumsSTA3sq=zeros(nwin,2*mshift+1);
    sumsSTA3Bsq=zeros(nwin,2*mshift+1); % refers the shifting SSI for SIL is not moving?
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %  sumsPGSS is shorter than sumsPGC2 by 2*mshift.  This is why sumsPGC2 etc
    %  is shifted by +mshift.  cumsumPGSS(1,:)=cumsum(PGSSx)(1,:) starts mshift
    %  to the right of the first data sample.  igstart is how many to the right
    %  of that.
    %  6/29/2018:  I'm pretty sure I want to subtract mshift from every iend or
    %  istart index to the right of the equal signs in the following FOR loop.
    %  The following looks like it (previously) had the proper relative shifts but the
    %  absolute time wasn't registered to the global igstart properly.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for n=1:nwin
        istart=igstart+(n-1)*winoff;     % 2*sps+1 + (n-1)*3*sps, igstart is the index of the starting sample; istart is index of each sample window
        iend=istart+winlen-1;              % + 12.5*sps, ADD -1 by Chao, 2019/02/17
        timswin(n)=timsSTA(istart+winlen/2);    % timsSTA, time serie data of STA; timswin, center of win, also == (istart+iend)/2
        sumsSTA12(n,:)=sum(STA12x(istart-mshift: iend-mshift, :));
        sumsSTA13(n,:)=sum(STA13x(istart-mshift: iend-mshift, :));
        sumsSTA32(n,:)=sum(STA32x(istart-mshift: iend-mshift, :));
        sumsSTA1sq(n,:)=sum(STAauto(1, istart: iend));
        sumsSTA3Bsq(n,:)=sum(STAauto(3, istart: iend));
        for m=-mshift:mshift
            sumsSTA2sq(n,m+mshift+1)=sum(STAauto(2, istart-m: iend-m)); %+m??? (yes).
            sumsSTA3sq(n,m+mshift+1)=sum(STAauto(3, istart-m: iend-m));
        end
        
%         sumsSTA12(n,:)=cumsumSTA12(iend-mshift,:)-cumsumSTA12(istart-mshift-1,:);  %Yes, -mshift (6/29/18), sum of square of points between istart-mshift and iend-shift
%         sumsSTA13(n,:)=cumsumSTA13(iend-mshift,:)-cumsumSTA13(istart-mshift-1,:);  %%% result in the sum from istart-mshift to iend-mshift
%         sumsSTA32(n,:)=cumsumSTA32(iend-mshift,:)-cumsumSTA32(istart-mshift-1,:);
%         sumsSTA1sq(n,:)=STAsq(1,iend)-STAsq(1,istart-1);  %PGC2 is cumsummed. Yes, +mshift.  No, no mshift (6/29/18)
%         sumsSTA3Bsq(n,:)=STAsq(3,iend)-STAsq(3,istart-1); %Similar, for the SILB-SSIB connection, here SILB does not shift
%         for m=-mshift:mshift
%             sumsSTA2sq(n,m+mshift+1)=STAsq(2,iend-m)-STAsq(2,istart-1-m); %+m??? (yes).
%             sumsSTA3sq(n,m+mshift+1)=STAsq(3,iend-m)-STAsq(3,istart-1-m);
%         end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %  Denominator for the normalization.  A 2D array, nwin by 2*mshift+1.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %An attempt to bypass glitches in data.  Min value of good data typically ~10^{-2}
    glitches=1.e-7;
    sumsSTA1sq=max(sumsSTA1sq,glitches);    
    sumsSTA2sq=max(sumsSTA2sq,glitches);
    sumsSTA3sq=max(sumsSTA3sq,glitches);    % return maximum between A and B
    %
    denomSTA12n=realsqrt(sumsSTA1sq.*sumsSTA2sq);    % Real square root, An error is produced if X is negative
    denomSTA13n=realsqrt(sumsSTA1sq.*sumsSTA3sq);
    denomSTA32n=realsqrt(sumsSTA3Bsq.*sumsSTA2sq);
    %
    sumsSTA12n=sumsSTA12./denomSTA12n;   % suffix 'n' means normalized
    sumsSTA13n=sumsSTA13./denomSTA13n;
    sumsSTA32n=sumsSTA32./denomSTA32n;
    [xcmaxSTA12n,imaxSTA12]=max(sumsSTA12n,[],2);   %Integer-offset max cross-correlation
    [xcmaxSTA13n,imaxSTA13]=max(sumsSTA13n,[],2);   % along row, max cc val and index in each window
    [xcmaxSTA32n,imaxSTA32]=max(sumsSTA32n,[],2);
    %Parabolic fit:
    [xmaxSTA12n,ymaxSTA12n,aSTA12]=parabol(nwin,mshift,sumsSTA12n,imaxSTA12); %Interpolated max cross-correlation
    [xmaxSTA13n,ymaxSTA13n,aSTA13]=parabol(nwin,mshift,sumsSTA13n,imaxSTA13);
    [xmaxSTA32n,ymaxSTA32n,aSTA32]=parabol(nwin,mshift,sumsSTA32n,imaxSTA32);

    %h=figure('Position',[0.1*wid 1 2.5*wid hite]); %center

    ix=sub2ind(size(denomSTA12n),(1:nwin)',imaxSTA12); %Find the linear index of the largest denominator
    ampSTA12=sqrt(denomSTA12n(ix)); %This makes amplitude linear rather than quadratic with counts.   % JUST FOR EASIER USAGE
    ampSTA1sq=sumsSTA1sq(ix); %by construction PGC2 is the same for all shifts  % sumsPGC2 becomes sumsSTA1sq
    ampSTA2sq=sumsSTA2sq(ix); % sumsSSIB2 becomes sumsSTA2sq, NOTICE: here ampSTA1sq are still sum of sqaures, i.e., quadratic
    ix=sub2ind(size(denomSTA13n),(1:nwin)',imaxSTA13);
    ampSTA13=sqrt(denomSTA13n(ix));
    ampSTA3sq=sumsSTA3sq(ix);
    ix=sub2ind(size(denomSTA32n),(1:nwin)',imaxSTA32);
    ampSTA32=sqrt(denomSTA32n(ix));
    AmpComp(1:4)=0;       % amplitude compare
    %AmpComp seems to be amplitude squared in 4s window minus amp squared in prior window,
    %divided by sum of amp squared in the two windows.  And why?
    AmpComp(5:nwin)=((ampSTA1sq(5:nwin)+ampSTA2sq(5:nwin)+ampSTA3sq(5:nwin))- ...
                     (ampSTA1sq(1:nwin-4)+ampSTA2sq(1:nwin-4)+ampSTA3sq(1:nwin-4)))./ ...
                    ((ampSTA1sq(5:nwin)+ampSTA2sq(5:nwin)+ampSTA3sq(5:nwin))+ ...
                     (ampSTA1sq(1:nwin-4)+ampSTA2sq(1:nwin-4)+ampSTA3sq(1:nwin-4)));
    
    %%%%%%%%%%%%
    % add by Chao, to clear unnecessary variables to free memory in case of
    % crash
    clear denomSTA12n denomSTA13n denomSTA32n sumsSTA12 sumsSTA13 sumsSTA32
    clear sumsSTA1sq sumsSTA2sq sumsSTA3sq sumsSTA3Bsq
    %%%%%%%%%%%%
                    
    %Center them
    imaxSTA12cent=imaxSTA12-mshift-1;  % "cent" is "centered"; imaxSTA12 is original 1: 2*mshift+1, corresponds to -mshift: mshift
    imaxSTA13cent=imaxSTA13-mshift-1;  
    imaxSTA32cent=imaxSTA32-mshift-1;
    %%% NOTICE: the right order of a closed 3-sta pair is +13, -12, +32, where 13 means 1-->3 
    iloopoff=imaxSTA13cent-imaxSTA12cent+imaxSTA32cent; %How well does the integer loop close?
    %
    xmaxSTA12n=xmaxSTA12n-mshift-1;
    xmaxSTA13n=xmaxSTA13n-mshift-1;
    xmaxSTA32n=xmaxSTA32n-mshift-1;
    loopoff=xmaxSTA13n-xmaxSTA12n+xmaxSTA32n; %How well does the interpolated loop close?
    xcmaxAVEn=(xcmaxSTA12n+xcmaxSTA13n+xcmaxSTA32n)/3;     % arithmetic average, == x-corr max average normalized
    % xcnshifts=cputime-t
    % t=cputime;
    %%% ampmax == max amplitude among all windows of all 3 2-station pairs
    ampmax=max([ampSTA12; ampSTA13; ampSTA32]);  % ';' means another row, size of the concatenate is 3*nwin, 2*mshift+1
    medxcmaxAVEn=median(xcmaxAVEn);
    xmaxSTA12ntmp=xmaxSTA12n;    % tmp == temporary
    xmaxSTA13ntmp=xmaxSTA13n;
    xmaxSTA32ntmp=xmaxSTA32n;

    %% Detection
    % iup=4;
    nin=0;      % subscript flag
    zerosallowed=20*winlen/160;   % 20*12.5*sps/160, BUT what is the meaning of this eqn, 5 zeros /sec == 5*12.5 /win ???  
    concentration=1.1; %in seconds; how concentrated is the coherent energy within the window?
    cncntr=concentration*sps;   % in samples, 20
    offset=round(0.5*cncntr);   % +- 1/2*cncntr, 10 samples, 0.25s
    for n=1:nwin
%         n=1;

        %%%%%%%%%%%%%%% Detection Rejection Criteria %%%%%%%%%%%%%%%%%%%%%%
        % 1. if < min threshold == 0.3, or
        % 2. abs(loopoff) > loopoffmax == 2.1, or
        % 3. imaxSTA12/13/32 is 1 or 2*mshift+1, which is located at the edge of the range [1, 2*mshift+1], or
        % 4. too much zeros in the trace
        % NEW for low-freq detection:
        % 5. imaxSTA12/13/32cent > mshift-5, i.e. 19, 5 is in case of cycle skipping
        %   || abs(imaxSTA12cent(n))> mshift-5 || abs(imaxSTA13cent(n))> mshift-5 ...
        %   || abs(imaxSTA32cent(n))> mshift-5
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %%% xmaxSTA13n(n)-xmaxSTA12n(n)+xmaxSTA32n(n) == loopoff(n)        
        if xcmaxAVEn(n)<xcmaxAVEnmin || abs(loopoff(n))>loopoffmax ...
                || isequal(abs(imaxSTA12cent(n)),mshift) || isequal(abs(imaxSTA13cent(n)),mshift) ...
                || isequal(abs(imaxSTA32cent(n)),mshift) || max(STAnzeros(:,n))>zerosallowed
            xmaxSTA12ntmp(n)=mshift+1; xmaxSTA13ntmp(n)=mshift+1; xmaxSTA32ntmp(n)=mshift+1; %dummy them, if these criteria are met
        else
            xcmaxconprev=-99999.;  %used to be 0; not good with glitches
            %%% NOTE here are using offset imaxSTA12/13/32 without centering, i.e. [1, 2shift+1]  
            imaxSTA12n=imaxSTA12(n); %This "n" for nth window; other "n's" for "normalized".  Unfortunately.
            imaxSTA13n=imaxSTA13(n);
            imaxSTA32n=imaxSTA32(n);
            sumsSTA12nn=sumsSTA12n(n,:);   % xxx'n' for normalized,
            sumsSTA13nn=sumsSTA13n(n,:);
            sumsSTA32nn=sumsSTA32n(n,:);
            % Usually, the loop happens between imaxSTA12n +- floor(loopoffmax+1), 
            % width of which is 2*floor(loopoffmax+1)
            %%% floor(2.5)=2; floor(-2.6)=-3
            for iSTA12=max(1,imaxSTA12n-floor(loopoffmax+1)): min(imaxSTA12n+floor(loopoffmax+1),...
                                                                  2*mshift+1)
                for iSTA13=max(1,imaxSTA13n-floor(loopoffmax+1)): min(imaxSTA13n+floor(loopoffmax+1),...
                                                                  2*mshift+1)
                    ibangon = (mshift+1)-iSTA13+iSTA12;
                    %%% i.e., -mshift <= -iSTA13+iSTA12 <= mshift
                    if ibangon >= 1 && ibangon <= 2*mshift+1
                        xcmaxcon=sumsSTA12nn(iSTA12)+sumsSTA13nn(iSTA13)+sumsSTA32nn(ibangon);
                        if xcmaxcon > xcmaxconprev
                            xcmaxconprev=xcmaxcon;
                            iSTA12bang=iSTA12;
                            iSTA13bang=iSTA13;
                        end
                    end
                end
            end
            %%% will result in the max xcmaxcon and corresponding iSTA12,
            %%% iSTA13, and save them into xcmaxconprev, iSTA12bang and iSTA13bang
            
            iSTA32bang=(mshift+1)-iSTA13bang+iSTA12bang;
            tmp12 = iSTA12bang-(mshift+1);  % these 3 lines are added by Chao
            tmp13 = iSTA13bang-(mshift+1);
            tmp32 = iSTA32bang-(mshift+1);
            if abs(iSTA12bang-imaxSTA12n) <= loopoffmax && ...  %not sure if these 3 lines are satisfied automatically ...
               abs(iSTA13bang-imaxSTA13n) <= loopoffmax && ...  % SHOULD not be, i think, could be floor(loopoffmax+1) > loopoffmax 
               abs(iSTA32bang-imaxSTA32n) <= loopoffmax && ...           
               sumsSTA12n(n,iSTA12bang)+sumsSTA13n(n,iSTA13bang)+sumsSTA32n(n,iSTA32bang) >= ...
               3*xcmaxAVEnmin
               %%% ALSO, sumsSTA12n(n,iSTA12bang) == sumsSTA12nn(iSTA12bang)
               %%% sumsSTA12n(n,iSTA12bang)+sumsSTA13n(n,iSTA13bang)+sumsSTA32n(n,iSTA32bang)
               %%% == xcmaxconprev according to codes
               
                xmaxSTA12ntmp(n)=tmp12; %without interpolation this is just centering.
                xmaxSTA13ntmp(n)=tmp13;
                xmaxSTA32ntmp(n)=tmp32;
                
                %%% let us assume xmaxSTA12ntmp and so is the mean, mu of the distribution of the
                %%% real value, and now we need some estimates of standard deviation, sigma
                % 1st kind of estimate, the deviation between forced circuit and real circuit and
                % then take it as the same to both 12 and 13
                tmp12 = (imaxSTA12n-(mshift+1));
                tmp13 = (imaxSTA13n-(mshift+1));
                tmp32 = (imaxSTA32n-(mshift+1));
                sigma(nin+1) = abs(tmp13+tmp32-tmp12);
                
                % 2nd kind of estimate, the deviation between CC shift and forced shift in each
                % pair, 12 and 13, so sigma is different for each pair 
                sigma12(nin+1) = abs(iSTA12bang-imaxSTA12n);
                sigma13(nin+1) = abs(iSTA13bang-imaxSTA13n);
                
                
                %for plotting traces
                imaxSTA12wr=round(xmaxSTA12ntmp(n)); %without interpolation this is not needed.
                imaxSTA13wr=round(xmaxSTA13ntmp(n));
    % 
                istart=igstart+(n-1)*winoff; %+mshift; %a better way might exist?  %ADDED mshift 10/20/12; DELETED IT 1/19/17.
                    %ADDED IT BACK 10/4/2017 to fix bug.  PGC is offset from igstart by mshift before first x-correlation.
                    %Not sure why mshift was added.  It changes STA12tr, STA1file etc. relative to the window that was used
                    %in the original x-correlation.  This will affect the stated time of max energy (through idiff).
                    %GOT RID of the mshift, yet again, 6/29/18, but only after subtracing mshift from all those istarts and
                    %iends in lines 342-355.
                iend=istart+winlen-1;
                imid=round((istart+iend)/2);
                %Check power spectrum for reasonableness
                %%% pwelch is a built-in function, [Pxx F] = pwelch(X, WINDOW, NOVERLAP, NFFT, Fs)
                [STA1xx fp] = pwelch(STAopt(1,istart:iend),[],[],[],sps); %40 is sps   
                STA1xx=STA1xx/max(STA1xx);    % normalization
                [STA2xx fp] = pwelch(STAopt(2,istart-imaxSTA12wr:iend-imaxSTA12wr),[],[],[],sps);  
                STA2xx=STA2xx/max(STA2xx);
                [STA3xx fp] = pwelch(STAopt(3,istart-imaxSTA13wr:iend-imaxSTA13wr),[],[],[],sps);
                STA3xx=STA3xx/max(STA3xx);
                flo=find(fp > lo,1)-1;    % find(fp > lo,1) finds the first 1 indice that satisfies fp > lo
                fhi=find(fp > hi,1)+1;    %extra 1 for good measure
                belowcut=median([STA1xx(2:flo); STA2xx(2:flo); STA3xx(2:flo)]);   
                ppeaksSTA1=findpeaks(STA1xx(flo+1:fhi));   % PKS = findpeaks(Y) finds local peaks in the data vector Y
                if length(ppeaksSTA1)>=1
                    maxppeakSTA1=max(ppeaksSTA1);
                else
                    maxppeakSTA1=0.;
                end
                ppeaksSTA2=findpeaks(STA2xx(flo+1:fhi));   % for STA2, use exactly the same procedure as STA1 
                if length(ppeaksSTA2)>=1
                    maxppeakSTA2=max(ppeaksSTA2);
                else
                    maxppeakSTA2=0.;
                end
                ppeaksSTA3=findpeaks(STA3xx(flo+1:fhi));   % for STA3, still the same
                if length(ppeaksSTA3)>=1
                    maxppeakSTA3=max(ppeaksSTA3);
                else
                    maxppeakSTA3=0.;
                end
                abovecut=median([maxppeakSTA1 maxppeakSTA2 maxppeakSTA3]);   % relative to belowcut, remain [belowcut, abovecut]
                if abovecut > 0.9*belowcut %-1 %This checks for frequency range; make sure it's not too narrow?
                    STA12tr=STAopt(1,istart:iend).*STAopt(2,istart-imaxSTA12wr:iend-imaxSTA12wr);   % see line 554, imax is already centered
                    STA13tr=STAopt(1,istart:iend).*STAopt(3,istart-imaxSTA13wr:iend-imaxSTA13wr);   % see line 556
                    STA32tr=STAopt(3,istart-imaxSTA13wr:iend-imaxSTA13wr).*...
                    STAopt(2,istart-imaxSTA12wr:iend-imaxSTA12wr);
                    cumsumtr=cumsum(STA12tr)+cumsum(STA13tr)+cumsum(STA32tr);    % sum of the cumsum of all traces
                    %%% first get the squared sum of each 0.5s window, then get the maximum and start indice
                    %%% this gives the proxy of Energy of 0.5s window, i.e. the integral of \dot(E)
                    %%% of eqn 1 in Rubin et al. 2013
                    %%% idiff<==>win [idiff+1, idiff+20] in regional ind<==>win [istart+idiff, istart+idiff+19] in global index
                    [cumsumtrdiff, idiff]=max(cumsumtr(cncntr+1:winlen)-cumsumtr(1:winlen-cncntr));
                    

                    %%%% Added by Chao, implemention of R_o/o in Peng et al. 2015
                    imdiff = istart+idiff+offset;
                    isdiff = imdiff-cncntr; %Start of strongest 0.5s, DELETE -1 by Chao, 2019/02/17, see NOTES
                    iediff = imdiff+cncntr-1;  % 2*cncntr is 40sps ~= 1s in current hf

%                     isdiff = istart+idiff; %Start of strongest 0.5s, DELETE -1 by Chao, 2019/02/17, see NOTES
%                     iediff = istart+idiff-1+2*cncntr; 
                    
                    STA11opt = STAopt(1,isdiff:iediff).* STAopt(1,isdiff:iediff);
                    STA22opt = STAopt(2,isdiff-imaxSTA12wr:iediff-imaxSTA12wr).* ...
                               STAopt(2,isdiff-imaxSTA12wr:iediff-imaxSTA12wr);
                    STA33opt = STAopt(3,isdiff-imaxSTA13wr:iediff-imaxSTA13wr).* ...
                               STAopt(3,isdiff-imaxSTA13wr:iediff-imaxSTA13wr);
                    STAsumopt = STA11opt+STA22opt+STA33opt;      
                    STA11ort = STAort(1,isdiff:iediff).* STAort(1,isdiff:iediff);
                    STA22ort = STAort(2,isdiff-imaxSTA12wr:iediff-imaxSTA12wr).* ...
                               STAort(2,isdiff-imaxSTA12wr:iediff-imaxSTA12wr);
                    STA33ort = STAort(3,isdiff-imaxSTA13wr:iediff-imaxSTA13wr).* ...
                               STAort(3,isdiff-imaxSTA13wr:iediff-imaxSTA13wr);
                    STAsumort = STA11ort+STA22ort+STA33ort;        
                    Eopt = sum(STAsumopt);  % in hf, it is ~ 1s interval
                    Eort = sum(STAsumort);
                    Roo = Eopt/Eort;
                    
                    Roomin = 1.5;                    
                    %%%% Added by Chao, implemention of R_o/o in Peng et al. 2015
                    
                    if Roo >= Roomin
                        
                        %%% xcmaxAVEnbang is added by Chao, to distinguish from
                        %%% xcmaxAVEn, because it is max average CC coef
                        xcmaxAVEnbang(nin+1)=(sumsSTA12n(n,iSTA12bang)+...
                                              sumsSTA13n(n,iSTA13bang)+...
                                              sumsSTA32n(n,iSTA32bang))/3;
                                          
%                         % display the number of overlapping window, easy for testing
%                         disp(n)
%                         keyboard
                
%                     %%% NOTICE! HERE, istart=igstart+(n-1)*winoff >= igstart,
%                     %%% So, isdiff >= istart >= igstart 
%                     isdiff=istart+idiff; %Start of strongest 0.5s, DELETE -1 by Chao, 2019/02/17, see NOTES
%                     iediff=istart+idiff-1+cncntr;  % cncntr is 20sps == 0.5s
%                     
%                     %What is amp squared in strongest coherent 1/2 sec?
%                     %%% amp squared is a self to self operation
%                     dummy=STAopt(1,isdiff:iediff).^2+ ...   % point square 
%                           STAopt(2,isdiff-imaxSTA12wr:iediff-imaxSTA12wr).^2+ ...
%                           STAopt(3,isdiff-imaxSTA13wr:iediff-imaxSTA13wr).^2;
%                     dum2=cumsum(dummy);
%                     Ampsq(nin+1)=dum2(end);   % Ampsq(1) == amplitude square = cumsum of all dummy = squared sum
                        Ampsq(nin+1)=0;                    

%                     %%Energy in prior 2.5*1.25s == 1.25*sps == 2.5*cncntr, with offset (assuming 0.5cncntr)
%                     %%% energy is a self to self operation, and is proportional to amp squared
%                     if isdiff > round(2.5*cncntr)+mshift+offset
%                         dummy=STAopt(1,isdiff-round(2.5*cncntr)-offset:isdiff-offset).^2+ ...
%                               STAopt(2,isdiff-round(2.5*cncntr)-imaxSTA12wr-offset:isdiff-imaxSTA12wr-offset).^2+ ...
%                               STAopt(3,isdiff-round(2.5*cncntr)-imaxSTA13wr-offset:isdiff-imaxSTA13wr-offset).^2;
%                     else
%                         dummy=STAopt(1,mshift:isdiff-1-offset).^2+ ...
%                               STAopt(2,mshift-imaxSTA12wr:isdiff-imaxSTA12wr-1-offset).^2+ ...
%                               STAopt(3,mshift-imaxSTA13wr:isdiff-imaxSTA13wr-1-offset).^2;
%                     end
%                     dum2=cumsum(dummy);
% %                     Prev(nin+1)=dum2(end);    % Prev(1) == previous 1.25s window before the strongest window, length 2.5*cncntr
%                     Prev(nin+1)=0;
%                     clear dummy
                        Prev(nin+1)=0;
                    
%                     %CC in same window (test)
%                     %%% cc is a cross-station operation
%                     dummy(1,:)=STAopt(1,istart:iend);
%                     dummy(2,:)=STAopt(2,istart-imaxSTA12wr:iend-imaxSTA12wr);
%                     dummy(3,:)=STAopt(3,istart-imaxSTA13wr:iend-imaxSTA13wr);
%                     denoms=dot(dummy,dummy,2);    % dot(A,B,DIM) returns the summed scalar product of A and B in the dimension DIM, 2 is row 
%                     cc(nin+1)=(dot(dummy(1,:),dummy(2,:))/sqrt(denoms(1)*denoms(2))+dot(dummy(2,:),dummy(3,:))/sqrt(denoms(2)*denoms(3))+ ...
%                                dot(dummy(3,:),dummy(1,:))/sqrt(denoms(3)*denoms(1)))/3;
%                     clear dummy
                    cc(nin+1)=0;    
                    
%                     %CC in prior 2.5*1.25 seconds, with offset
%                     if isdiff > round(2.5*cncntr)+mshift+offset                        
%                         dummy(1,:)=STAopt(1,isdiff-round(2.5*cncntr)-offset:isdiff-1-offset);   % 2.5*cnctr == 1.25s *sps
%                         dummy(2,:)=STAopt(2,isdiff-imaxSTA12wr-round(2.5*cncntr)-offset:isdiff-imaxSTA12wr-1-offset);
%                         dummy(3,:)=STAopt(3,isdiff-imaxSTA13wr-round(2.5*cncntr)-offset:isdiff-imaxSTA13wr-1-offset);
%                     else
%                         dummy(1,:)=STAopt(1,mshift:isdiff-1-offset);   
%                         dummy(2,:)=STAopt(2,mshift-imaxSTA12wr:isdiff-imaxSTA12wr-1-offset);
%                         dummy(3,:)=STAopt(3,mshift-imaxSTA13wr:isdiff-imaxSTA13wr-1-offset);
%                     end                      
%                     denoms=dot(dummy,dummy,2);   % dot(A,B,DIM) returns the summed scalar product of A and B in the dimension DIM, 2 is row
%                     ccprior125(nin+1)=(dot(dummy(1,:),dummy(2,:))/sqrt(denoms(1)*denoms(2))+dot(dummy(2,:),dummy(3,:))/sqrt(denoms(2)*denoms(3))+ ...
%                                        dot(dummy(3,:),dummy(1,:))/sqrt(denoms(3)*denoms(1)))/3;  % ccprior125(1) means cc in prior 1.25s
%                     clear dummy
                        ccprior125(nin+1)=0;
                    
%                     %CC in prior 4 seconds, with offset
%                     if isdiff > 4*sps+mshift+offset   % >winsec*sps+19+10
%                         dummy(1,:)=STAopt(1,isdiff-4*sps-offset:isdiff-1-offset);   % this is 12.5s much more than 4s, maybe the previous winlen is 4s
%                         dummy(2,:)=STAopt(2,isdiff-imaxSTA12wr-4*sps-offset:isdiff-imaxSTA12wr-1-offset);
%                         dummy(3,:)=STAopt(3,isdiff-imaxSTA13wr-4*sps-offset:isdiff-imaxSTA13wr-1-offset);
%                     else
%                         dummy(1,:)=STAopt(1,mshift:isdiff-1-offset);   % 19: isdiff-1    % MIGHT this part be a mistake? no symmetry, with offset
%                         dummy(2,:)=STAopt(2,mshift-imaxSTA12wr:isdiff-imaxSTA12wr-1-offset);
%                         dummy(3,:)=STAopt(3,mshift-imaxSTA13wr:isdiff-imaxSTA13wr-1-offset);
%                     end
%                     denoms=dot(dummy,dummy,2);
%                     ccprior(nin+1)=(dot(dummy(1,:),dummy(2,:))/sqrt(denoms(1)*denoms(2))+dot(dummy(2,:),dummy(3,:))/sqrt(denoms(2)*denoms(3))+ ...
%                                     dot(dummy(3,:),dummy(1,:))/sqrt(denoms(3)*denoms(1)))/3;
%                     clear dummy
                        ccprior(nin+1)=0;                    
                    
%                     %CC in following 2.5*1.25 seconds, with 2*offset
%                     if iediff+2.5*cncntr+2*offset <= size(STAopt,2) && ...
%                        iediff-imaxSTA12wr+2.5*cncntr+2*offset <= size(STAopt,2) && ...
%                        iediff-imaxSTA13wr+2.5*cncntr+2*offset <= size(STAopt,2)
%                         dummy(1,:)=STAopt(1,iediff+1+2*offset:iediff+round(2.5*cncntr)+2*offset);  % 1.25s win after 0.5s after 2*offset
%                         dummy(2,:)=STAopt(2,iediff+1-imaxSTA12wr+2*offset:iediff-imaxSTA12wr+round(2.5*cncntr)+2*offset);
%                         dummy(3,:)=STAopt(3,iediff+1-imaxSTA13wr+2*offset:iediff-imaxSTA13wr+round(2.5*cncntr)+2*offset);
%                         denoms=dot(dummy,dummy,2);
%                         ccpost125(nin+1)=(dot(dummy(1,:),dummy(2,:))/sqrt(denoms(1)*denoms(2))+dot(dummy(2,:),dummy(3,:))/sqrt(denoms(2)*denoms(3))+ ...
%                                       dot(dummy(3,:),dummy(1,:))/sqrt(denoms(3)*denoms(1)))/3;
%                     else
%                         ccpost125(nin+1)=0;
%                     end    
%                     clear dummy
                        ccpost125(nin+1)=0;
                    
%                   %CC in following 4 seconds
%                   if iediff < winbig-winlen-mshift
%                       dummy(1,:)=STAopt(1,iediff+1:iediff+winlen);
%                       dummy(2,:)=STAopt(2,iediff+1-imaxSTA12wr:iediff-imaxSTA12wr+winlen);
%                       dummy(3,:)=STAopt(3,iediff+1-imaxSTA13wr:iediff-imaxSTA13wr+winlen);
%                   else
%                       dummy(1,:)=STAopt(1,iediff+1:winbig-mshift);
%                       dummy(2,:)=STAopt(2,iediff+1-imaxSTA12wr:winbig-mshift-imaxSTA12wr);
%                       dummy(3,:)=STAopt(3,iediff+1-imaxSTA13wr:winbig-mshift-imaxSTA13wr);
%                   end
%                   denoms=dot(dummy,dummy,2);
%                   ccpost(nin+1)=(dot(dummy(1,:),dummy(2,:))/sqrt(denoms(1)*denoms(2))+dot(dummy(2,:),dummy(3,:))/sqrt(denoms(2)*denoms(3))+ ...
%                                 dot(dummy(3,:),dummy(1,:))/sqrt(denoms(3)*denoms(1)))/3;
%                   clear dummy

%                     %Energy in prior 15 seconds
%                     %%% STAauto is dot product, an array, why use median?
%                     if isdiff > 15*sps+mshift
%                         dummy=median(STAauto(1,isdiff-15*sps:isdiff))+ ...    % sum of median of 15-s win
%                               median(STAauto(2,isdiff-15*sps-imaxSTA12wr:isdiff-imaxSTA12wr))+ ...
%                               median(STAauto(3,isdiff-15*sps-imaxSTA13wr:isdiff-imaxSTA13wr));
%                     else
%                         dummy=median(STAauto(1,1:isdiff))+ ...
%                               median(STAauto(2,1:isdiff-imaxSTA12wr))+ ...
%                               median(STAauto(3,1:isdiff-imaxSTA13wr));
%                     end
%                     Prev15(nin+1)=dummy;
%                     clear dummy
                        Prev15(nin+1)=0;
                    
%                     %Energy in prior 30 seconds
%                     if isdiff > 30*sps+mshift
%                         dummy=median(STAauto(1,isdiff-30*sps:isdiff))+ ...
%                               median(STAauto(2,isdiff-30*sps-imaxSTA12wr:isdiff-imaxSTA12wr))+ ...
%                               median(STAauto(3,isdiff-30*sps-imaxSTA13wr:isdiff-imaxSTA13wr));
%                     else
%                         dummy=median(STAauto(1,1:isdiff))+ ...
%                               median(STAauto(2,1:isdiff-imaxSTA12wr))+ ...
%                               median(STAauto(3,1:isdiff-imaxSTA13wr));
%                     end
%                     Prev30(nin+1)=dummy;
%                     clear dummy
                        Prev30(nin+1)=0;
                    
%                     % energy (amp squared sum) in post 1.25s window
%                     dummy=STAopt(1,iediff:iediff+round(2.5*cncntr)).^2+ ...   % 1.25s after iediff
%                           STAopt(2,iediff-imaxSTA12wr:iediff+round(2.5*cncntr)-imaxSTA12wr).^2+ ...
%                           STAopt(3,iediff-imaxSTA13wr:iediff+round(2.5*cncntr)-imaxSTA13wr).^2;
%                     dum2=cumsum(dummy);
%                     Post(nin+1)=dum2(end);
%                     clear dummy
                        Post(nin+1)=0;
                    
%                     %%%%%%%For comparison, start
%                     %
%                     %%% isdiff=istart+idiff-1, so indx == isdiff+round(cncntr/2)
%                     %%% round(cncntr/2) == offset
%                     indx=istart+idiff-1+round(cncntr/2); %indx should be centered ~ on zero-crossing of main arrival
%                     traces(1,:)=STAopt(1,indx-tempbef:indx+tempaft);   % -59;+60
%                     traces(2,:)=STAopt(2,indx-tempbef-imaxSTA12wr:indx+tempaft-imaxSTA12wr);
%                     traces(3,:)=STAopt(3,indx-tempbef-imaxSTA13wr:indx+tempaft-imaxSTA13wr);
%                     for ista=1:nsta
%                         % tempxc == template x-corr, |max lag| == 0.25s == 10 samples
%                         % cc between template and the window with the main arrival
%                         tempxc(ista,:)=xcorr(traces(ista,:),STAtemps(ista,:),floor(cncntr/2),'coeff');
%                     end
%                     sumxc=sum(tempxc)/nsta;    % sum along column, average cc coeff of all 3 stats at each shift
%                     [match(nin+1,1), ioff]=max(sumxc);  % get the max coeff and the index of the shift
%                     %%% the follow line is to convert index to real shift value
%                     ioff=ioff-(floor(cncntr/2)+1); %shift STAtemps by match(nin+1) (shift right for positive values)
%                     %%% timstemp == times of template, 1: 120, 121: 240
%                     timstemp(nin*templen+1:(nin+1)*templen)=timsSTA(indx-tempbef+ioff:indx+tempaft+ioff);
%                     for ista=1:nsta
%                         % tempxcneg means x-corr with negative template
%                         tempxcneg(ista,:)=xcorr(traces(ista,:),-STAtemps(ista,:),floor(cncntr/2),'coeff'); 
%                     end
%                     sumxcneg=sum(tempxcneg)/nsta;    
%                     [match(nin+1,2), ioff]=max(sumxcneg);
%                     ioff=ioff-(floor(cncntr/2)+1); %shift "snips" by match(nin+1) (shift right for positive values), what is snips ???
%                     timstempneg(nin*templen+1:(nin+1)*templen)=timsSTA(indx-tempbef+ioff:indx+tempaft+ioff);
%                     %
%                     %%%%%%%For comparison, end
                        match(nin+1,1)=0;
                        match(nin+1,2)=0;

                        % STA1file, STA2file, STA3file are back-to-back windows for the day (plus sample times)
                        % 1st column: time
                        % 2nd column: data
                        STA1file(nin*winlen+1:(nin+1)*winlen,1:2)=[timsSTA(istart:iend)' ...
                            STAopt(1,istart:iend)'];
                        STA2file(nin*winlen+1:(nin+1)*winlen,1:2)=[timsSTA(istart:iend)' ...
                            STAopt(2,istart-imaxSTA12wr:iend-imaxSTA12wr)'];
                        STA3file(nin*winlen+1:(nin+1)*winlen,1:2)=[timsSTA(istart:iend)' ...
                            STAopt(3,istart-imaxSTA13wr:iend-imaxSTA13wr)'];
                        STAamp(nin+1,1)=prctile(abs(STAopt(1,istart:iend)),80);  % prctile,  Percentiles of a sample, SEE NOTES
                        STAamp(nin+1,2)=prctile(abs(STAopt(2,istart-imaxSTA12wr:iend-imaxSTA12wr)),80);
                        STAamp(nin+1,3)=prctile(abs(STAopt(3,istart-imaxSTA13wr:iend-imaxSTA13wr)),80);
                        STAamp(nin+1,:)=STAamp(nin+1,:)/STAamp(nin+1,1);   % normalization relative to 1st stat
                        STA1bbfile(nin*winlen+1:(nin+1)*winlen,1:2)=[timsSTA(istart:iend)' ...
                            STAoptbb(1,istart:iend)'];  % bb is broader band
                        STA2bbfile(nin*winlen+1:(nin+1)*winlen,1:2)=[timsSTA(istart:iend)' ...
                            STAoptbb(2,istart-imaxSTA12wr:iend-imaxSTA12wr)'];
                        STA3bbfile(nin*winlen+1:(nin+1)*winlen,1:2)=[timsSTA(istart:iend)' ...
                            STAoptbb(3,istart-imaxSTA13wr:iend-imaxSTA13wr)'];
                        %%% imaxSTA12wr     -> cc shift of qualified win only, so remain no changed until the next qualified win;
                        %%% xcmaxSTA12n(n)  -> max cc coef in all 4s win
                        %%% cumsumtrdiff/cumsumtr(winlen) -> normalized cumsumtrdiff
                        STA12file(nin+1,1:2)=[imaxSTA12wr xcmaxSTA12n(n)];  
                        STA13file(nin+1,1:2)=[imaxSTA13wr xcmaxSTA13n(n)];
                        STA32file(nin+1,1:3)=[cumsumtrdiff/cumsumtr(winlen) xcmaxSTA32n(n) idiff];

                        nin=nin+1;
                        istartkeep(nin)=istart; %For adding other stations later, keep the istart of each win
                        aSTA12keep(nin,:)=[timswin(n) aSTA12(n)];   % timswin == time at the center of each win 
                        aSTA13keep(nin,:)=[timswin(n) aSTA13(n)];   
                        aSTA32keep(nin,:)=[timswin(n) aSTA32(n)];
                        loopoffkeep(nin,:)=[timswin(n) loopoff(n)];
                    
                        %%% 1. Without interpolation, xmaxSTA13ntmp == imaxSTA13wr
                        %%% 2. xcmaxAVEn is arithmetic average, == x-corr coeff max average normalized in each win 
                        %%% 3. Ampsq -> amplitude square == squared sum in each strongest 0.5s win
                        %%% 4. cumsumtrdiff -> the maximum squared sum of among all 0.5s windows
                        %%% 5. timswin(n)-winlensec/2+idiff/sps -> start time in sec of the strongest 0.5s win
                        %%% 6. cumsumtrdiff/cumsumtr(winlen) -> normalized cumsumtrdiff
                        %%% 7. match(1) -> max cc coeff between template and the window with the main arrival, 2s win, 120 samples
                        %%% 8. match(2) -> max cc coeff between negtive template and the window with the main arrival, 2s win, 120 samples
                        %%% 9. Prev -> amp squared sum (energy) in previous 1.25s window before the strongest window
                        %%% 10. Post -> amp squared sum (energy) in post 1.25s window before the strongest window
                        %%% 11. ccprior -> cc coeff in previous 4 secs window
                        %%% 12. STAamp -> normalized 80 percentile of amplitude of data
                        %%% 13. xcmaxSTA12n -> max cc coeff in each 4s window
                        %%% 14. sumsSTA12n(n,ceil(iSTA12bang/iup)) -> the CC coeff at detection timing
                        %%% 15. sigma, sigma12, sigma13 -> the different standard deviation estiamte
                    
%                     mapfile(nin,:)=[timswin(n) xmaxSTA13ntmp(n) xmaxSTA12ntmp(n) ...
%                         xcmaxAVEn(n) loopoff(n) Ampsq(nin) cumsumtrdiff timswin(n)-winlensec/2+idiff/sps cumsumtrdiff/cumsumtr(winlen) ...
%                         match(nin,1) match(nin,2) Prev(nin) Post(nin) Prev15(nin) Prev30(nin) ...
%                         ccprior125(nin) ccprior(nin) ccpost125(nin) STAamp(nin,2) STAamp(nin,3) xcmaxSTA12n(n) xcmaxSTA13n(n) xcmaxSTA32n(n) ];
                    
                        %%% NOTICE: add by Chao, save xcmaxAVEnbang instead of
                        %%% xcmaxAVEn, and the other 3 ccmax's
                    
                        %%% 26 cols, 2020/08/27
                        mapfile(nin,:)=[timswin(n) xmaxSTA13ntmp(n) xmaxSTA12ntmp(n) xcmaxAVEnbang(nin)...
                                        loopoff(n) Ampsq(nin) cumsumtrdiff ...
                                        timswin(n)-winlensec/2+(idiff+1)/sps ...
                                        cumsumtrdiff/cumsumtr(winlen) ...
                                        match(nin,1) match(nin,2) Prev(nin) Post(nin) Prev15(nin) ...
                                        Prev30(nin) ccprior125(nin) ccprior(nin) ccpost125(nin) ...
                                        STAamp(nin,2) STAamp(nin,3) sumsSTA12n(n,iSTA12bang) ...
                                        sumsSTA13n(n,iSTA13bang) sumsSTA32n(n,iSTA32bang) ...
                                        sigma(nin) sigma12(nin) sigma13(nin)];
                    else
                        xmaxSTA12ntmp(n)=mshift+1; xmaxSTA13ntmp(n)=mshift+1; xmaxSTA32ntmp(n)=mshift+1;
                    end                
                else
                    % 20 == mshift+1
                    xmaxSTA12ntmp(n)=mshift+1; xmaxSTA13ntmp(n)=mshift+1; xmaxSTA32ntmp(n)=mshift+1;
                end
            else
                xmaxSTA12ntmp(n)=mshift+1; xmaxSTA13ntmp(n)=mshift+1; xmaxSTA32ntmp(n)=mshift+1; 
            end
        end
    end
    
    nin
    
    if nin == 0
        fprintf('No detection is found in fam %s in day %s %s, no plots&files would be saved.\n', ...
                fam, YEAR, JDAY);
    else
        
    %% Write results into files
    % save analytics of all detections (nin), contain double counting
    fid = fopen([rstpath,'/MAPS/map',IDENTIF,'_',num2str(lo),'-',num2str(hi),'_',...
                 num2str(winlen/sps),'s',num2str(sps),'sps'],'w+');
    fprintf(fid,'%9.1f %6.2f %6.2f %8.3f %7.2f %10.3e %10.3e %10.3f %7.3f %7.3f %7.3f %10.3e %10.3e %10.3e %10.3e %7.3f %7.3f %7.3f %5.2f %5.2f %7.3f %7.3f %7.3f %6.2f %6.2f %6.2f \n',mapfile(1:nin,:)');
    fclose(fid);
    
    % save time traces of detections (nin), contain double counting
    fid = fopen(strcat(rstpath, '/MAPS/trace_',IDENTIF,'_',num2str(lo),'-',num2str(hi),'_',...
                       num2str(winlen/sps),'s',num2str(sps),'sps'),'w+');
    tracefile = [STA1file(1:nin*winlen,:) STA2file(1:nin*winlen,2) STA3file(1:nin*winlen,2)];
    fprintf(fid,'%.4f %.6f %.6f %.6f \n',tracefile');
    fclose(fid);
    
    % save time traces of all divided wins (nwin), original rotated
    % velocity seismogram, after bandpassing
    fid = fopen(strcat(rstpath, '/MAPS/trace1day_',IDENTIF,'_',num2str(lo),'-',num2str(hi),...
                       '_',num2str(winlen/sps),'s',num2str(sps),'sps'),'w+');
    oritracefile = [timsSTA' STAopt'];
    fprintf(fid,'%.4f %.6f %.6f %.6f \n',oritracefile');
    fclose(fid);
    
    
    %% Plot
    %%% 1st figure, all detections represented by offsets against time
    % The follow 4 subplots are basically the same, dividing the time axis
    % into 4 parts
    figure(101)
    subplot(4,1,1,'align'); 
    hold on
    plot(timswin,xcmaxAVEnmin*mshift*40/sps+zeros(nwin,1),'k:');
    plot(timsSTA(winlen:2*winlen),7+zeros(winlen+1,1),'k','linewidth',2);
    plot(timswin,zeros(nwin,1),'k:');
    plot(timswin,xcmaxAVEn*mshift*40/sps,'g');
    plot(timswin,xmaxSTA12ntmp*40/sps,'bs','MarkerSize',2);
    plot(timswin,xmaxSTA13ntmp*40/sps,'ro','MarkerSize',2);
%     plot(timswin,xmaxSTA32ntmp,'k*','MarkerSize',2);
    axis([0 86400/4 -40 40]);
    ylabel('samples')
    title([IDENTIF,'_{',num2str(lo),'-',num2str(hi),'}'])
    box on
    
    subplot(4,1,2,'align'); 
    hold on
    plot(timswin,xcmaxAVEnmin*mshift*40/sps+zeros(nwin,1),'k:');
    plot(timswin,zeros(nwin,1),'k:');
    plot(timswin,xcmaxAVEn*mshift*40/sps,'g');
    plot(timswin,xmaxSTA12ntmp*40/sps,'bs','MarkerSize',2);
    plot(timswin,xmaxSTA13ntmp*40/sps,'ro','MarkerSize',2);
%     plot(timswin,xmaxSTA32ntmp,'k*','MarkerSize',2);
    axis([86400/4 86400/2 -40 40]);
    ylabel('samples')
    box on
    
    subplot(4,1,3,'align'); 
    hold on
    plot(timswin,xcmaxAVEnmin*mshift*40/sps+zeros(nwin,1),'k:');
    plot(timswin,zeros(nwin,1),'k:');
    plot(timswin,xcmaxAVEn*mshift*40/sps,'g');
    plot(timswin,xmaxSTA12ntmp*40/sps,'bs','MarkerSize',2);
    plot(timswin,xmaxSTA13ntmp*40/sps,'ro','MarkerSize',2);
%     plot(timswin,xmaxSTA32ntmp,'k*','MarkerSize',2);
    axis([86400/2 3*86400/4 -40 40]);
    ylabel('samples')
    box on
    
    subplot(4,1,4,'align'); 
    hold on
    plot(timswin,xcmaxAVEnmin*mshift*40/sps+zeros(nwin,1),'k:');
    plot(timswin,zeros(nwin,1),'k:');
    plot(timswin,xcmaxAVEn*mshift*40/sps,'g');
    plot(timswin,xmaxSTA12ntmp*40/sps,'bs','MarkerSize',2);
    plot(timswin,xmaxSTA13ntmp*40/sps,'ro','MarkerSize',2);
%     plot(timswin,xmaxSTA32ntmp,'k*','MarkerSize',2);
    axis([3*86400/4 86400 -40 40]);
    xlabel('sec')
    ylabel('samples')
    box on
    orient landscape    % is used to set up the paper orientation of a Figure or Model window for printing
    print('-depsc',[rstpath,'/FIGS/',IDENTIF,'_',num2str(winlen/sps),'s_',num2str(sps),'sps_',...
                    num2str(lo),'-',num2str(hi),'b.eps'])
    close(101)
    
%     %%% Try to figure out which coef is more proper to use in figure 102,
%     %%% it seems that ymaxSTA12n is problematic, though it is more
%     %%% reasonable to use
%     %%% xcmaxSTA... is the max cc coef before parabolic fit
%     %%% ymaxSTA... is the max cc coef after parabolic fit
%     figure
% %     plot(xcmaxSTA12n, 'k-'); hold on
% %     plot(ymaxSTA12n,'r-');
% %     plot(xcmaxSTA13n, 'k:');hold on
% %     plot(ymaxSTA13n,'r:');
% %     plot(xcmaxSTA32n, 'k--');hold on
% %     plot(ymaxSTA32n,'r--');    
%     plot((xcmaxSTA12n+xcmaxSTA13n+xcmaxSTA32n)/3, 'k-.'); hold on
%     plot((ymaxSTA12n+ymaxSTA13n+ymaxSTA32n)/3, 'r-.');

%     %%% Try to figure which CC coef is more appropriate to save in mapfile
%     %%% xcmaxAVEn(nin), i.e. mapfile(:,4), the max cc coef hasn't met the
%     %%% criteria that offset 13+32-12=0.
%     %%% xcmaxAVEnbang, the cc coef of that offset when meet the criteria.
%     figure
%     plot(mapfile(:,4), 'k--'); hold on
%     plot(xcmaxAVEnbang, 'r--'); 
%     box on
%     %%% they are almost the same, but when meet the criteria, cc coef is
%     %%% not necessarily the max cc coef, <=, and very close
%     %%% AND, xcmaxAVEnbang is more reasonable to save 
    
    %%% 2nd figure, distribution of cc coefs and offsets
    figure(102)
    colormap(jet)
    %%% loopoff == xmaxSTA13n-xmaxSTA12n+xmaxSTA32n
    %%% AmpComp -> difference of (sum of squared amplitude in 4s win) between each win and its previous 4th win 
    %%% scatter, 3 is the marker size, AmpComp is the color of the marker
    scatter(xmaxSTA13n-xmaxSTA12n+xmaxSTA32n,xcmaxAVEn,3,AmpComp)
    hold on 
    plot(-50:50,xcmaxAVEnmin+zeros(101,1),'k:');
    axis([min(-5,-2.5*loopoffmax) max(5,2.5*loopoffmax) -0.2 1.0])
    hrf = plotreflinesr(gca,-loopoffmax,'x','k');colorbar   % plots dotted lines at the locations specified
    hrf = plotreflinesr(gca,loopoffmax,'x','k');colorbar
    box on
    %%% Tex syntax to use subscript, '_{}'
    title([IDENTIF,'_{',num2str(lo),'-',num2str(hi),'}'])
    print('-depsc',[rstpath,'/FIGS/',IDENTIF,'-',num2str(winlen/sps),'s_',num2str(sps),'sps_',...
                    num2str(lo),'-',num2str(hi),'e.eps'])
    close(102)    

    
    %% METHOD 1 to add additional stations for checking 
    %%%%%%%%%%%%%%%%%%%%%%%%%% METHOD 1, Begin %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    clear STAopt    % clear these two variables is just to save memory space
    clear STAort
    clear STAnzeros
%     stasnew=['LZB  '
%              'TWKB '
%              'MGCB '];  % twkb lzb mgcb
    stasnew=['PGC  '
             'SSIB '
             'SILB '
             'KLNB '];  % twkb lzb mgcb   
         
    if year==2003   % in 2003, station KLNB is named KELB, so use KELB to replace KLNB in 2003
        stasnew(4,:)='KELB ';
    else
        stasnew(4,:)='KLNB ';  % remember to change it back
    end
    
    nstanew=size(stasnew,1);

    fileflag = ones(nstanew,1);     % file flag to indicate if file exists
    STAopt=zeros(nstanew,nin*winlen);
    STAort=zeros(nstanew,nin*winlen);
    for istanew=1:nstanew
        found=0;
        [LIA,idx]=ismember(stasnew(istanew,:),PERMSTA,'rows');
        if LIA
            found=found+LIA;
            if strcmp(PERMSTA(idx,1:3),'PGC')
                fact=1.0e-3;
            elseif strcmp(PERMSTA(idx,1:3),'LZB')
                fact=1.8e-3;
            end
            fname = strcat(prename,'.',PERMSTA(idx,1:3),'..BHE.D.SAC');
            if isfile(fname)    % if have the data file
%                 [opt,ort,nzeros,timsSTAnew]=readperms(prename,PERMSTA,PERMROTS,idx,sps,lo,hi,npo,npa,fact,nwin,winlen,winoff,igstart);
                [opt,ort,nzeros,timsSTAnew]=readpermsv2(prename,PERMSTA,PERMROTS,idx,sps,lo,hi,...
                                                        npo,npa,fact,nwin,winlen,winoff,igstart);
            else
                fileflag(istanew) = 0;   % change the file flag to 0, meaning abnormal
                fprintf('No data for additional station %s in day %s / %s, will be omitted. \n',...
                        PERMSTA(idx,1:3), YEAR, JDAY);
                continue   % continue to next staion
            end         
        end
        
        [LIA,idx]=ismember(stasnew(istanew,:),POLSTA,'rows');
        if LIA
            found=found+LIA; %better be 1
            if year==2003 && jday<213
                fact=7.5e-3;
            else
                fact=1.5e-3; 
            end
            fname = strcat(prename,'.',POLSTA(idx,1:4),'..HHE.D.SAC');
            if isfile(fname)    % if have the data file
%                 [opt,ort,nzeros]=readpols(prename,POLSTA,POLROTS,idx,sps,lo,hi,npo,npa,fact,nwin,winlen,winoff,igstart);
                [opt,ort,nzeros]=readpolsv2(prename,POLSTA,POLROTS,idx,sps,lo,hi,npo,npa,fact,...
                                            nwin,winlen,winoff,igstart);
            else
                fileflag(istanew) = 0;   % change the file flag to 0, meaning abnormal
                fprintf('No data for additional station %s in day %s / %s, will be omitted. \n',...
                        POLSTA(idx,1:4), YEAR, JDAY);
                continue   % continue to next staion
            end 
        end
%         found=found
        if fileflag(istanew)
            for n=1:nin
                istart=istartkeep(n);   % saved istart of each successful window
                iend=istart+winlen-1;
                STAopt(istanew,(n-1)*winlen+1:n*winlen)=opt(istart:iend)';   % no need to divide scaleseisms here because they are all 1
                STAort(istanew,(n-1)*winlen+1:n*winlen)=ort(istart:iend)';
            end
        end
    end
    
    stasnew=stasnew(fileflag == 1,:);      
    nstanew=size(stasnew,1);
    STAopt=STAopt(fileflag == 1,:); 
    
    in=zeros(nstanew,nin);   % flag to indicate whether this win of this new sta is qualified (1)
    ioff=zeros(nstanew,nin);    % some kind of weighted offset according to the ccmax coef
    ccmaxave=zeros(nstanew,nin);
    loff=zeros(nstanew,nin);    % average sum of abs offset difference
    icc12=zeros(nstanew,nin);
    icc13=zeros(nstanew,nin);
    icc23=zeros(nstanew,nin);
    count=zeros(nstanew,1);
    %%% NOTE: for mshiftnew 's performance, 2.2 > 2.0 > 1.5 > 1.0,
    %%% meaningless to increase more
    mshiftnew=ceil(2.2*mshift); %1.8 was used in first-year report   % ceil(a) gives integer >=a, 
    ccavethres = 0.8*xcmaxAVEnmin;  %% average cc threshold
    ccindthres = 0.6*xcmaxAVEnmin;  %% individual cc threshold
    offavethres = loopoffmax;    % average differential offset threshold
    offindthres = 1.5*loopoffmax;   
    
    for istanew=1:nstanew   % 1st, for each additional station
        
        for n=1:nin     % 2nd, for each trio detection 
            
            istart=(n-1)*winlen+1;  % start from 1
            iend=n*winlen;
            % recall that STA1/2/3file saves the time and shifted seismogram of each station
            %%% do xcorr between 1 new station to each of the 3-sta trio
            [cc1, lag1] = xcorr(STA1file(istart:iend,2),STAopt(istanew,istart:iend),mshiftnew,'coeff');
            [cc2, lag2] = xcorr(STA2file(istart:iend,2),STAopt(istanew,istart:iend),mshiftnew,'coeff');
            [cc3, lag3] = xcorr(STA3file(istart:iend,2),STAopt(istanew,istart:iend),mshiftnew,'coeff');
            % get the max cc coef and index
            [cc1max,icc1]=max(cc1);
            [cc2max,icc2]=max(cc2);
            [cc3max,icc3]=max(cc3);
            ccmaxave(istanew,n)=(cc1max+cc2max+cc3max)/3;
            
            % NOTE:
            % In each deteciton, since sta 1,2,3 is already aligned to
            % denote one point, thus if it can be checked by the 4th station
            % the shifts between 4th station and
            % original 3 should be close enough, i.e. all differential shifts
            % should be close to 0 theoratically. That is why it is a
            % constraint to filter suspicious detections
            icc12(istanew,n) = lag1(icc1)-lag2(icc2);
            icc13(istanew,n) = lag1(icc1)-lag3(icc3);
            icc23(istanew,n) = lag2(icc2)-lag3(icc3);
            
            % loff and lofftest are exactly the same, although different syntax, depends on user habit
            % loff the absolute sum of the offset difference            
            loff(istanew,n)=(abs(icc12(istanew,n))+abs(icc13(istanew,n))+abs(icc23(istanew,n)))/3;
            
            ioff(istanew,n)=round((cc1max * lag1(icc1) + ...
                                   cc2max * lag2(icc2) + ...
                                   cc3max * lag3(icc3)) / ...
                                  (cc1max+cc2max+cc3max));
            %%% NOTE:
            %%% in first year report, the following criteria was used, so did the hf counterpart 
            %%% if abs(icc12)<=2 && abs(icc13)<=2 && abs(icc23)<=2 && ccmaxave(istanew,n)>1*xcmaxAVEnmin
            if abs(icc12(istanew,n))<=offindthres && abs(icc13(istanew,n))<=offindthres && ...
                    abs(icc23(istanew,n))<=offindthres && ...
                    ccmaxave(istanew,n)>=ccavethres && ... 
                    loff(istanew,n)<=offavethres && ...
                    cc1max >= ccindthres && ...
                    cc2max >= ccindthres && ...
                    cc3max >= ccindthres 
                
                in(istanew,n)=1; 
                count(istanew) = count(istanew)+1;                                
                                  
                %Now shift new station (for plotting purposes)
                if ioff(istanew,n)>0    % shift to right -->
                    STAopt(istanew,istart+ioff(istanew,n):iend)=STAopt(istanew,istart:iend-...
                                                                       ioff(istanew,n));
                    STAopt(istanew,istart:istart+ioff(istanew,n)-1)=0;
                elseif ioff(istanew,n)<0    % shift to left <--
                    STAopt(istanew,istart:iend+ioff(istanew,n))=STAopt(istanew,istart-...
                                                                       ioff(istanew,n):iend);
                    STAopt(istanew,iend+ioff(istanew,n)+1:iend)=0;                    
                end
                
            end
            
        end
        
    end
    
%     %% plot relationship between ave cc of 3-trio and ave cc with the 4th station
%     figure('Position',[scrsz(3)/10 1.5*scrsz(4)/10 1.5*scrsz(3)/5 4*scrsz(4)/5]);
%     ax=NaN(nstanew,1);
%     for istanew=1:nstanew
%         ax(istanew)=subplot(nstanew, 1, istanew,'align');
%         ind=find(loff(istanew,:)<=offavethres);
% %         scatter(ccmaxave(istanew,ind),xcmaxAVEnbang(ind), 6, loff(istanew,ind), 'filled');
%         scatter(ccmaxave(istanew,1:nin),xcmaxAVEnbang(1:nin), 6, abs(ioff(istanew,1:nin)), 'filled');
%         oldcmap = colormap(ax(istanew), hot);
%         colormap(ax(istanew), flipud(oldcmap) );
%         colorbar
%         caxis([0, 5]);
%         box on;
%         grid on;
%         axis equal
% %         xlim([-30, 30]);
% %         ylim([0.1, 0.8]);
% %         xticks(-10:2:10);
%         title(strcat(IDENTIF,'_Additional Station',{' '},stasnew(istanew,:)));
%         xlabel('Average CC max between 4th sta and trio pairs');
%         ylabel('Average CC max between trio pairs');
%     end    
%         
%     %% plot CC results between 3-sta trio and additional stations
%     
%     %%% median may be proper for ccmax, but not for icc
%     medcc = median(ccmaxave,2);     % median is a fair measure, maybe 0.8*xcmaxAVEnmin is similar
%     medicc12 = prctile(abs(icc12),30,2);  % 2 is a proper cutoff for all 3 additional stas
%     medicc13 = prctile(abs(icc13),30,2);
%     medicc23 = prctile(abs(icc23),30,2);
% 
%     figure('Position',[scrsz(3)/10 scrsz(4)/10 4*scrsz(3)/5 4*scrsz(4)/5]);
%     for istanew=1:nstanew
%         subplot(nstanew, nstanew, istanew,'align');
% %         plot(icc1(istanew,:)-mshiftnew-1, cc1max(istanew,:), 'ro','Markersize',8); hold on
% %         plot(icc2(istanew,:)-mshiftnew-1, cc2max(istanew,:), 'bs','Markersize',8);
% %         plot(icc3(istanew,:)-mshiftnew-1, cc3max(istanew,:), 'k+','Markersize',8);
%         plot(icc12(istanew,:), ccmaxave(istanew,:), 'ro','Markersize',8); hold on
%         plot([-30 30],[medcc(istanew) medcc(istanew)],'g--', 'linewidth',2);
%         plot([-medicc12(istanew) -medicc12(istanew)], [0.1 0.8], 'c--', 'linewidth',2);
%         plot([medicc12(istanew) medicc12(istanew)], [0.1 0.8], 'c--', 'linewidth',2);
%         xlim([-30, 30]);
%         ylim([0.1, 0.8]);
%         xticks(-10:2:10);
%         title(strcat(IDENTIF,'_Additional Station',{' '},stasnew(istanew,:)));
%         xlabel('Differential offset index 1-2');
%         ylabel('Average cc max');
% 
% 
%         subplot(nstanew, nstanew, nstanew+istanew,'align');
%         plot(icc13(istanew,:), ccmaxave(istanew,:), 'bs','Markersize',8); hold on
%         plot([-30 30],[medcc(istanew) medcc(istanew)],'g--', 'linewidth',2);
%         plot([-medicc13(istanew) -medicc13(istanew)], [0.1 0.8], 'c--', 'linewidth',2);
%         plot([medicc13(istanew) medicc13(istanew)], [0.1 0.8], 'c--', 'linewidth',2);
%         xlim([-30, 30]);
%         ylim([0.1, 0.8]);
%         xticks(-10:2:10);
%         xlabel('Differential offset index 1-3');
%         ylabel('Average cc max');
% 
%         subplot(nstanew, nstanew, 2*nstanew+istanew,'align');
%         plot(icc23(istanew,:), ccmaxave(istanew,:), 'k+','Markersize',8); hold on
%         plot([-30 30],[medcc(istanew) medcc(istanew)],'g--', 'linewidth',2);
%         plot([-medicc23(istanew) -medicc23(istanew)], [0.1 0.8], 'c--', 'linewidth',2);
%         plot([medicc23(istanew) medicc23(istanew)], [0.1 0.8], 'c--', 'linewidth',2);
%         xlim([-30, 30]);
%         ylim([0.1, 0.8]);
%         xticks(-10:2:10);
%         xlabel('Differential offset index 2-3');
%         ylabel('Average cc max');
%     end
    %%%%%%%%%%%%%%%%%%%%%%%%%% METHOD 1, End %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    %% Write checking result into files
    addrstfile = [in(:,1:nin); loff(:,1:nin) ; ioff(:,1:nin); ccmaxave(:,1:nin)]';      % number of columns should be 4*nstanew 
    fid = fopen([rstpath,'/MAPS/mapadd',IDENTIF,'_',num2str(lo),'-',num2str(hi),'_',num2str(winlen/sps),'s',num2str(sps),'sps', num2str(nstanew), 'nsta'],'w+');
    if nstanew == 3  
        fprintf(fid,'%d %d %d %.4f %.4f %.4f %d %d %d %.3f %.3f %.3f \n',addrstfile');  % 12 col
    elseif nstanew ==4
        fprintf(fid,'%d %d %d %d %.4f %.4f %.4f %.4f %d %d %d %d %.3f %.3f %.3f %.3f \n',addrstfile');   % 16 col
    end
    fclose(fid);
    
    fid = fopen(strcat(rstpath, '/MAPS/traceadd_',IDENTIF,'_',num2str(lo),'-',num2str(hi),'_',num2str(winlen/sps),'s',num2str(sps),'sps', num2str(nstanew), 'nsta'),'w+');
    addtracefile = [STA1file(1:nin*winlen,1) STAopt(:,1:nin*winlen)'];
    if nstanew == 3
        fprintf(fid,'%.4f %.6f %.6f %.6f \n',addtracefile');    % 4 col
    elseif nstanew == 4
        fprintf(fid,'%.4f %.6f %.6f %.6f %.6f \n',addtracefile');    % 5 col
    end
    fclose(fid);
    
    
    %% WAY 2 to avoid double counting, consider all additional stations together
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% WAY 2 start %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% construct a more complete result matrix of all nin windows
    %%% structure of allrst: 13+4*nstanew cols, if 4 new stas, then will be 32 cols
    %%% UPDATED at 2020/08/27    
    %%%   1:timswin(n) 2:xmaxSTA13ntmp(n) 3:xmaxSTA12ntmp(n) 4:xcmaxAVEnbang(nin) 5:loopoff(n)
    %%%   6:cumsumtrdiff 7:timswin(n)-winlensec/2+idiff/sps 8:cumsumtrdiff/cumsumtr(winlen)
    %%%   9:match(nin,1) 10:match(nin,2) 11:sumsSTA12n(n,iSTA12bang) 12:sumsSTA13n(n,iSTA13bang)
    %%%   13:sumsSTA32n(n,iSTA32bang) in(1:nstanew) loff(1:nstanew) ioff(1:nstanew) ccmaxave(1:nstanew)
    %%%   sigma(nin) sigma12(nin) sigma13(nin)  
    allrst = [mapfile(1:nin, 1:5) mapfile(1:nin, 7:11) mapfile(1:nin, 21:23) addrstfile ...
              mapfile(1:nin, 24:26)];
    dtmin = concentration;      % min time during which only 1 detection is retained
    indext = min(ceil(winlensec/winoffsec)+4, size(allrst,1));   % index extension, the index range that needs to be checked
    if nstanew == 3
        colnum = [7 4 23 24 25];    % col numbers that are timing and cc coefs
    elseif nstanew == 4
        colnum = [7 4 26 27 28 29];
    end
    allrstsort = sortrows(allrst, [7,1]);
    [allrsttemp1,~,~] = RemoveDoubleCounting(allrstsort,dtmin,indext,colnum);
    % do this twice to avoid extreme case like 20.2 20.5 20.5 20.6 20.9,
    % when first round saves the 2nd 20.5 and 20.9, but 20.5 and 20.9 are
    % still duplicates
    [allrsttemp2,~,~] = RemoveDoubleCounting(allrsttemp1,dtmin,indext,colnum); 
    bbb = allrsttemp2(2:end,7)-allrsttemp2(1:end-1,7);
    ccc = allrsttemp2(3:end,7)-allrsttemp2(1:end-2,7);
    if ~isempty(find(bbb<=dtmin,1)) || ~isempty(find(ccc<=dtmin,1))
        disp('not enough');
        [allrst_new,~,~] = RemoveDoubleCounting(allrsttemp2,dtmin,indext,colnum); 
    else
        allrst_new = allrsttemp2;
    end
    [~, isave, ~] = intersect(allrst, allrst_new, 'row', 'stable');     % saved indexes
    idis = setdiff(1:nin,isave);    % discarded indexes
    
    %%% re-assign those variables that need to be saved
    ind_ori=find(xmaxSTA12ntmp ~= mshift+1);
    idis_ori=ind_ori(idis);
    xmaxSTA12ntmp(idis_ori) = mshift+1;
    xmaxSTA13ntmp(idis_ori) = mshift+1;
    xmaxSTA32ntmp(idis_ori) = mshift+1;
    
    nin_new = length(isave);
    for i = 1: nin_new
        STA1file_new((i-1)*winlen+1: i*winlen,:)= STA1file((isave(i)-1)*winlen+1 : isave(i)*winlen, :);
        STA2file_new((i-1)*winlen+1: i*winlen,:)= STA2file((isave(i)-1)*winlen+1 : isave(i)*winlen, :);
        STA3file_new((i-1)*winlen+1: i*winlen,:)= STA3file((isave(i)-1)*winlen+1 : isave(i)*winlen, :);
        STAopt_new(1:nstanew,(i-1)*winlen+1: i*winlen)= STAopt(1:nstanew, (isave(i)-1)*winlen+1 : isave(i)*winlen);
        
        STAamp_new(i,:) = STAamp(isave(i),:);
        
        STA1bbfile_new((i-1)*winlen+1: i*winlen,:)= STA1bbfile((isave(i)-1)*winlen+1 : isave(i)*winlen, :);
        STA2bbfile_new((i-1)*winlen+1: i*winlen,:)= STA2bbfile((isave(i)-1)*winlen+1 : isave(i)*winlen, :);
        STA3bbfile_new((i-1)*winlen+1: i*winlen,:)= STA3bbfile((isave(i)-1)*winlen+1 : isave(i)*winlen, :);
        
        
        STA12file_new(i,:) = STA12file(isave(i), :);
        STA13file_new(i,:) = STA13file(isave(i), :);
        STA32file_new(i,:) = STA32file(isave(i), :);
        
        istartkeep_new(i) = istartkeep(isave(i));
        aSTA12keep_new(i,:) = aSTA12keep(isave(i), :);
        aSTA13keep_new(i,:) = aSTA13keep(isave(i), :);
        aSTA32keep_new(i,:) = aSTA32keep(isave(i), :);
        loopoffkeep_new(i,:) = loopoffkeep(isave(i), :);
        
        
    end
    
    
    %%% Write results into files
    fid = fopen([rstpath,'/MAPS/mapallrst',IDENTIF,'_',num2str(lo),'-',num2str(hi),'_',num2str(winlen/sps),'s',num2str(sps),'sps',num2str(nstanew), 'nsta'],'w+');
    if nstanew == 3
        fprintf(fid,'%9.1f %6.2f %6.2f %8.3f %7.2f %10.3e %10.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %d %d %d %.4f %.4f %.4f %d %d %d %.3f %.3f %.3f %6.2f %6.2f %6.2f \n',allrst_new(1:nin_new,:)');
    elseif nstanew == 4
        fprintf(fid,'%9.1f %6.2f %6.2f %8.3f %7.2f %10.3e %10.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %d %d %d %d %.4f %.4f %.4f %.4f %d %d %d %d %.3f %.3f %.3f %.3f %6.2f %6.2f %6.2f \n',allrst_new(1:nin_new,:)');
    end
    fclose(fid);
    
    fid = fopen(strcat(rstpath, '/MAPS/traceall_',IDENTIF,'_',num2str(lo),'-',num2str(hi),'_',num2str(winlen/sps),'s',num2str(sps),'sps',num2str(nstanew), 'nsta'),'w+');
    tracefile = [STA1file_new(1:nin_new*winlen,:) STA2file_new(1:nin_new*winlen,2) STA3file_new(1:nin_new*winlen,2) STAopt_new(:, 1:nin_new*winlen)'];
    if nstanew == 3
        fprintf(fid,'%.4f %.6f %.6f %.6f %.6f %.6f %.6f \n',tracefile');
    elseif nstanew == 4
        fprintf(fid,'%.4f %.6f %.6f %.6f %.6f %.6f %.6f %.6f \n',tracefile');
    end
    fclose(fid);
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% WAY 2 end %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    toc
    
        
%     %% Plot
%     close all
%     nt=0;
%     nrow=4;
%     if winlen<200
%         mcol=6;
%     elseif winlen<400
%         mcol=4;
%     else
%         mcol=2;
%     end
%     
%     for ifig=1:floor(nin_new/(nrow*mcol))+1
%         figure('Position',[scrsz(3)/10 scrsz(4)/10 4*scrsz(3)/5 9*scrsz(4)/10]);
%         fign(ifig) = gcf;
%         if ifig > 5
%             close(fign(ifig-5))
%         end
%         for n = 1:nrow
%             for m = 1:mcol
%                 nt=nt+1;
%                 if nt <= nin_new
%                      %if STA12file(nt,1) >= 10 && STA12file(nt,1) <= 16 && STA13file(nt,1) >= 2 && STA13file(nt,1) <= 8
%                      is = STA1file_new(winlen*(nt-1)+1,1);
%                      ien= STA1file_new(winlen*nt,1);
%                      %yma=0.4;
%                      yma=max(max([STA1file_new(winlen*(nt-1)+1:winlen*nt,2) STA2file_new(winlen*(nt-1)+1:winlen*nt,2) ...
%                          STA3file_new(winlen*(nt-1)+1:winlen*nt,2)]));
%                      ymi=min(min([STA1file_new(winlen*(nt-1)+1:winlen*nt,2) STA2file_new(winlen*(nt-1)+1:winlen*nt,2) ...
%                          STA3file_new(winlen*(nt-1)+1:winlen*nt,2)]));
%                      xvect=[is is+2*(yma-ymi)*(winlen/160.)]; %amplitude bar originally scaled for 4-s window
%                      yma=2.4*max(yma,-ymi);
%                      yvect=[-0.9*yma -0.9*yma];
%                      
%                      % plot checking result of all additional stations 
%                      for istanew = 1: nstanew
%                          subplot((nstanew+1)* nrow, mcol,((nstanew+1)*(n-1)+istanew-1)*mcol+m,'align');
%                          plot(STA1file_new(winlen*(nt-1)+1:winlen*nt,1),STA1file_new(winlen*(nt-1)+1:winlen*nt,2),'r');
%                          hold on
%                          plot(STA2file_new(winlen*(nt-1)+1:winlen*nt,1),STA2file_new(winlen*(nt-1)+1:winlen*nt,2),'b');
%                          plot(STA3file_new(winlen*(nt-1)+1:winlen*nt,1),STA3file_new(winlen*(nt-1)+1:winlen*nt,2),'k');
%                          plot(xvect,yvect,'r','linewidth',3);
%                          plot([is+1/hi is+1/lo],[-0.8*yma -0.8*yma],'k','linewidth',3);
%                          text(is+0.1*(ien-is), 0.66*yma, int2str(STA13file_new(nt,1)),'fontsize',6);
%                          text(ien-0.1*(ien-is), 0.66*yma, int2str(STA12file_new(nt,1)),'fontsize',6);
%                          if sum(in(:,nt))>=1
%                              plot(STA1file_new(winlen*(nt-1)+1:winlen*nt,1),STAopt_new(istanew,winlen*(nt-1)+1:winlen*nt),'c', 'linewidth',2);
%                              text(0.5*(is+ien), 0.66*yma, int2str(ioff(istanew,nt)),'fontsize',6);
%                              text(0.5*(is+ien)-0.1*(ien-is), -0.75*yma, num2str(ccmaxave(istanew,nt),3),'fontsize',6);
%                              text(0.5*(is+ien)+0.1*(ien-is), -0.75*yma, num2str(loff(istanew,nt),3),'fontsize',6);
%                          end
%                          if in(istanew,nt)==1
%                              plot(is+0.3*(ien-is), 0.66*yma, 'p');
%                          end
%                          %                           for istanew=1:nstanew
%                          %                               if (in(istanew,nt))
%                          %                                   plot(STA1file(winlen*(nt-1)+1:winlen*nt,1),STAoptnew(istanew,winlen*(nt-1)+1:winlen*nt),'c')
%                          %                                   text(0.5*(is+ien)-0.2, 0.66*yma, int2str(ioff(istanew,nt)),'fontsize',6);
%                          %                                   text(0.5*(is+ien)-0.2, -0.8*yma, num2str(ccmax(istanew,nt),3),'fontsize',6);
%                          %                               end
%                          %                           end
%                          box on
%                          axis([is ien -yma yma])
%                          set(gca,'XTick',[is (is+ien)/2],'fontsize',6);
%                      end
%                      
%                      % plot the original 3-station trio seimogram counterparts in broader frequency 
%                      subplot((nstanew+1)* nrow, mcol,((nstanew+1)*(n-1)+nstanew)*mcol+m,'align');
%                      plot(STA1bbfile_new(winlen*(nt-1)+1:winlen*nt,1),STA1bbfile_new(winlen*(nt-1)+1:winlen*nt,2),'r')
%                      hold on
%                      plot(STA2bbfile_new(winlen*(nt-1)+1:winlen*nt,1),STA2bbfile_new(winlen*(nt-1)+1:winlen*nt,2),'b')
%                      plot(STA3bbfile_new(winlen*(nt-1)+1:winlen*nt,1),STA3bbfile_new(winlen*(nt-1)+1:winlen*nt,2),'k')
%                      is = STA1bbfile_new(winlen*(nt-1)+1,1);
%                      ien= STA1bbfile_new(winlen*nt,1);
%                      plot([is+1/hibb is+1/lobb],[-0.8*yma -0.8*yma],'k','linewidth',3)
%                      box on
%                      axis([is ien -yma yma])
%                      set(gca,'XTick',[is (is+ien)/2],'fontsize',6);
%                      
%                 end
%             end
%         end
%         orient landscape
%         if ifig <= 9
%             print('-depsc',[rstpath,'/WIGS/',IDENTIF,'_','add','_',num2str(lo),'-',num2str(hi),'_',num2str(winlen/sps),'s_',num2str(sps),'sps','.',int2str(0),int2str(0),int2str(ifig),'.eps'])
%         elseif ifig <= 99
%             print('-depsc',[rstpath,'/WIGS/',IDENTIF,'_','add','_',num2str(lo),'-',num2str(hi),'_',num2str(winlen/sps),'s_',num2str(sps),'sps','.',int2str(0),int2str(ifig),'.eps'])
%         else
%             print('-depsc',[rstpath,'/WIGS/',IDENTIF,'_','add','_',num2str(lo),'-',num2str(hi),'_',num2str(winlen/sps),'s_',num2str(sps),'sps','.',int2str(ifig),'.eps'])
%         end
% 
%     end
    end
end

medlok=median(abs(loopoffkeep_new))     % median is actually effective to the 2nd column, which is loopoff
medaSTA12=median(aSTA12keep_new)        % 2nd column is aSTA12, which results from PARABOLIC FIT
medaSTA13=median(aSTA13keep_new)
medaSTA32=median(aSTA32keep_new)


% keyboard
end








