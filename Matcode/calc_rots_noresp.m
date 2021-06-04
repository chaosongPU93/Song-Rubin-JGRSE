function [ROTS]=calc_rots_noresp(fam,lo,hi,sps,bef,aft,splitoff,staoff,mainsta)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This is the script to get the shear-splitting parameter and rotation
% angles from scratch, following the method illustrated in NOTES. First to
% deal with the family 002 as a test, and will add other families in the
% future. 
%
% INPUT:
%   fam: fam, string
%   lo: lower corner freq. of passband
%   hi: higher corner freq. of passband
%   sps: sampling rate, use 40 if no resampling needed
%   bef: num of samples before the main dipole for estimating rots
%   aft: num of samples after the main dipole for estimating rots
%   splitoff: offset threshold in splitting
%   staoff: offset threshold between stations
%   mainsta: main station number, integer, ie. the reference sta when getting
%            the offsets between different sta
% OUTPUT:
%   ROTS: rotation parameters
%
% 
% Try to use the correct sign convention as Yajun did
% Make sure to understand this:
%
%   1. assume you have the N and E component data, and you can have the
%       complex data from this component. The shear wave spliting rotation
%       is try to find the two axes (slow and fast direction), on which the
%       projected record is the slow and fast component. So in fact, the
%       rotation is to ROTATE THE AXES (coordinate system), not the record,
%       counter-clockwise to get this two directions. Slow direction is defined 
%       as the rotated angle from E direction, whereas the fast direction
%       is the rotated angle from N direction. This angle is then accepted
%       as the 2nd column in ROTs parameters, ranges from [0,180) 
%   2. Projected on this axes, the fast and slow components should reach to
%       the maximum cc. Note the offset between slow and fast component, 
%       slow -fast (always positive) as the 1st column in ROTs parameter.
%       Only when the offset > 2 is recognized as obvious splitting,
%       otherwise it is forced to 0, as well as the angle.
%   3. After correcting the shear splitting effect (shift the slow one to left
%       by the offset), again, ROTATE THE AXES to find the optimal/orthogonal
%       direction. Rotate the axes by an angle counter-clockwise to get this 
%       two directions. Then the polarization angle is then recognized as
%       the angle from E again, and is the 3rd column in ROTs parameters,
%       ranges from [0,180)
%   4. After being rotated to optimal direction, do cc between different
%       stations to get the travel time differences. The overall difference
%       is recognized as the 4th column in ROTs parameter. The choice of
%       main station matters, since the quality of dipole affects the alignment
%       in CC.
%   5. MOST important thing to remember is, rotating the axes
%       counter-clockwise (using rotation matrix) by one angle then 
%       finding the new coordinates is numerially the same as,
%       rotating the record clockwise (using the complex rotation) by the 
%       same angle then get the new real and imaginary part. Check the code
%       'testrotation' to see the details. The second implemention is used
%        in Allan's code, and is more elegant.  
%       
%      
% USE data removed from station response
%
% NOTES:
%   1. there are fams that the stacked templates are pretty low-amplitude,
%   then the maximum offset allowed in CC would be necessary, if the result
%   reaches the boundary, then it is necessary to decrease the threshold.
%   However, check the N and E component to see if the splitting should be
%   small or large. For example, fam 010 should allow splitthres <=8, but
%   others could works fine when splitthres <=12, and usually this offset
%   ranges from 0~9.
%
%
% By Chao Song, chaosong@princeton.edu
% First created date:   2019/09/09
% Last modified date:   2019/11/10
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% default value for easy debugging
defval('fam', '002');
defval('lo', 0.5);
defval('hi', 6.5);
defval('sps', 40);
defval('bef', 25);
defval('aft', 35);
defval('splitoff', 12);
defval('staoff', 8);
defval('mainsta', 4);

%% Initialization
format short e   % Set the format to 5-digit floating point
% clear
% close all

% set(0,'DefaultFigureVisible','on');
% set(0,'DefaultFigureVisible','off');   % switch to show the plots or not
scrsz=get(0,'ScreenSize');   % for example, scrsz=[1,1,1680,1050]
wid=scrsz(3);       % width
hite=scrsz(4);       % height

% set path
workpath = getenv('ALLAN');

% datapath = workpath;
datapath = strcat(workpath,'/data-no-resp');

% old LFE family pool, from Yajun's selections
ofampool = ['002';  % 002,246
            '013';  % 043,152
            '025';  % 141
            '028';  % 047
            '056';  % 010
            '084';  % 144,102,121,58,257
            '099';  % 099,006
            '115';  % 068
            '125';  % 125
            '147';  % 147,17
            '149']; % 017,047

% the newest families from Bostock are different, so select the closest
% one as new family pool, also selected family have the most LFEs
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

nfam = length(nfampool);


%%% corresponds to PERMROTS
PERMSTA=['PGC'        % permanent station names
    'LZB'];
POLSTA=['SSIB '           % polaris station names
    'SILB '
    'KLNB '
    'MGCB '
    'TWKB '];

%  station trio used, 1st station is the main station
%     stas=['LZB  '
%         'PGC  '
%         'SSIB '
%         'SILB '
%         'TWKB '
%         'MGCB '
%         'KLNB '];
stas=['PGC  '
    'SSIB '
    'SILB '
    'LZB  '
    'TWKB '
    'MGCB '
    'KLNB '];

[timoffrot,tempoffs] = GetDays4Stack(fam);
tempoffs = round(tempoffs*sps/40);

% number of used stations
nsta=size(stas,1);         %  number of stations



%% stack at bostock's detection timings on N and E components

%%% data parameters
% sps=40;     % samples per second
lenofday = 24*3600;
STAE=zeros(nsta,sps * lenofday);
STAN=zeros(nsta,sps * lenofday);

%%% desired template parameters
templensec = 60;
templen = templensec * sps;
% before = templen/2-1;
% after = templen/2;
before = templen/8;
after = 7*templen/8-1;
extra = 6*sps;    % use extra more samples than the original window.
stackexE = zeros(nsta, templen+ 2*extra);
stackexN = zeros(nsta, templen+ 2*extra);

fam
nday = size(timoffrot, 1)
% timoffrot=timoffrot(5:9, :);
nday = size(timoffrot, 1);

if nday == 0
    disp("No day with detections found in this family of BOSTOCK's catalog");
    return
end


% get the all catalog LFEs in that family
bostname = ('/BOSTOCK/total_mag_detect_0000_cull_NEW.txt');
catalog = load(strcat(workpath, bostname));
famnum = str2double(fam);
dateall = catalog(famnum == catalog(:, 1), :);

%%% loop for day
nstack = 0;
nLFE = 0;
fprintf('Part 1: Raw Template \n');

for id = 1: nday
    %     nd = 1;
    year = timoffrot(id,1);
    YEAR = int2str(year);
    jday = timoffrot(id,2);
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
    nLFE = nLFE + size(bostsamp, 1);
    
    if jday <= 9
        JDAY=['00',int2str(jday)];
    elseif jday<= 99
        JDAY=['0',int2str(jday)];
    else
        JDAY=int2str(jday);
    end
    %     MO = 'SEP';
    MO=day2month(jday,year);     % EXTERNAL function, day2month, get the month of one particular date
    direc=[datapath, '/', YEAR,'/',MO,'/'];     % directory name
    fprintf('%s \n',direc);
    datafnm = [direc, YEAR,'.',JDAY,'.00.00.00.0000.CN'];
    fprintf('Processing date: %s / %s \n',YEAR, JDAY);
    nlfe1day = size(bostocks, 1)
    
    if year==2003   % in 2003, station KLNB is named KELB, so use KELB to replace KLNB in 2003
        POLSTA(3,:)='KELB ';
    else
        POLSTA(3,:)='KLNB ';  % remember to change it back
    end
    
    fileflag = 1;   % 1 means all files exist, the file status is normal
    
    %%% loop for each station for reading data
    for ista = 1: nsta
        found = 0;
        [LIA, idx] = ismember(stas(ista, :), PERMSTA, 'rows');
        
        %%% if station is a permanent one
        if LIA
            found = found+ LIA;
            if strcmp(PERMSTA(idx, 1:3),'PGC')
                fact = 1.0e-3;
            elseif strcmp(PERMSTA(idx, 1:3),'LZB')
                fact = 1.8e-3;
            end
            fname = strcat(datafnm,'.',PERMSTA(idx,1:3),'..BHE.D.SAC');
            if isfile(fname)    % if have the data file
                % 1. this is for data without removing station response
                %   [opt,ort,timeperm]=readperm_1daynofilterv2(datafnm, PERMSTA, PERMROTS, idx, sps, fact);
                % 2. this is for data with no response
                [dataE, dataN, ~] = readperms_norotsv2(datafnm, PERMSTA, idx, sps, fact);
            else
                fileflag = 0;   % change the file flag to 0, meaning abnormal
                fprintf('No data for station %s in day %s %s, this day will be omitted. \n',...
                        PERMSTA(idx,1:3), YEAR, JDAY);
                break   % break the entire station loop
            end
            
        end
        
        %%% if station is a polaris one
        [LIA, idx] = ismember(stas(ista, :), POLSTA, 'rows');
        if LIA
            found = found+ LIA;
            if year == 2003 && jday < 213
                fact = 7.5e-3;
            else
                fact = 1.5e-3;
            end
            fname = strcat(datafnm,'.',POLSTA(idx,1:4),'..HHE.D.SAC');
            if isfile(fname)    % if have the data file
                %                 % 1. this is for data without removing station response
                %                 [opt,ort,nzeros]=readpols(prename,POLSTA,POLROTS,idx,sps,lo,hi,npo,npa,fact,nwin,winlen,winoff,igstart);
                % 2. this is for data with no response
                [dataE, dataN, ~] = readpols_norotsv2(datafnm, POLSTA, idx, sps, fact);
            else
                fileflag = 0;   % change the file flag to 0, meaning abnormal
                fprintf('No data for station %s in day %s / %s. \n',POLSTA(idx,1:4), YEAR, JDAY);
                break   % break the entire station loop
            end
            
        end
        
        STAE(ista, :) = dataE;
        STAN(ista, :) = dataN;
    end
    
    if fileflag == 0    % means there are missing files
        fprintf('Day %s / %s will be omitted because of missing files. \n', YEAR, JDAY);
        continue    % continue to the next day
    end
    
    len = size(STAE,2);
    for n=1: size(bostsamp, 1)
        % set the timing be at the center of the window
        if (bostsamp(n)+ after+ extra <= len) && ...
                (bostsamp(n)- before- extra >= 1)
            
            %%% it is better to use a longer window than templen to make
            %%% sure the shifting in step2 would preserve the frequency
            %%% content as much as possible
            windata = STAE(:, bostsamp(n)- before- extra: ...
                bostsamp(n)+ after+ extra);
            % remove mean and linear trend
            tmp = detrend(windata');
            stackexE = stackexE+ tmp';
            
            windata2 = STAN(:, bostsamp(n)- before- extra: ...
                bostsamp(n)+ after+ extra);
            tmp2 = detrend(windata2');
            stackexN = stackexN+ tmp2';
            
            nstack = nstack + 1;
            
            % datamat is 3D, for storing all qualified wins for stacking
            % size is (nsta, nstack, templen)
            dmatE(:, nstack, :) = STAE(:, bostsamp(n)- before- extra: ...
                bostsamp(n)+ after+ extra);
            dmatN(:, nstack, :) = STAN(:, bostsamp(n)- before- extra: ...
                bostsamp(n)+ after+ extra);
            
        end
    end
    
end     % loop for days

stackE = stackexE(:, 1+extra: end-extra)/nstack;
stackN = stackexN(:, 1+extra: end-extra)/nstack;

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% filterring before all operations
% lo = 0.5;
% hi = 6.5;
npo = 2;
npa = 2;
for ista =1: nsta
    stackE(ista,:) = Bandpass(stackE(ista,:), sps, lo, hi, npo, npa, 'butter');
    stackN(ista,:) = Bandpass(stackN(ista,:), sps, lo, hi, npo, npa, 'butter');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


pltbef = 10*sps;
pltaft = 10*sps;
figure  % plot original stack, fig 1
subplot(2,1,1)  % plot E component 
for ista = 1:nsta
    plot(stackE(ista,tempoffs(ista)-pltbef+1: tempoffs(ista)+pltaft)*5 + 1* ista); hold on
end

subplot(2,1,2)  % plot N component
for ista = 1:nsta
    plot(stackN(ista,tempoffs(ista)-pltbef+1: tempoffs(ista)+pltaft)*5 + 1* ista); hold on
end

vertical_cursors;

% keyboard

%% rotate to get PERMROTS and POLAROTS parameters

%%%%%%%% Convention of each column in PERMROTS/POLROTS %%%%%%%%
%%% 1. offset in samples between fast/slow direction, unit in 40 sps, same
%%%     station, slow-fast, positive 
%%% 2. rotation angle to get fast/slow direction, same station, rotate this
%%%     angle from E counter-clockwise to get SLOW direction, from N cc to
%%%     get FAST direction, 0-180 deg
%%% 3. rotation angle to maximize the energy/particle motion/polarization,
%%%     same station, rotate this angle from E counter-clockwise to get
%%%     optimal component, 0-180 deg
%%% 4. offset in samples between the arrival times at different stations, unit in 40 sps
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

splitrot = 0:5:175;     % shear splitting rotation angle array, starting from E c-c, rotate the axes by this angle
polarrot = 0:1:179;     % polarization rotation angle array, starting from E c-c, rotate the axes by this angle

splitrotrad=splitrot*pi/180;     % convert angle to rad
polarrotrad=polarrot*pi/180;

ROTS = zeros(nsta, 4);  % rots parameter,

slowfastlagmat = zeros(nsta, length(splitrotrad));
slowfastcoefmat = zeros(nsta, length(splitrotrad));
optortratmat = zeros(nsta, length(polarrotrad));
for ista = 1:nsta
    
%     ista=2;
%     bef = 25;
%     aft = 35;
    seglen = bef+aft;
    
    % STA will be the same length, construct the complex record
    STA = stackE(ista,tempoffs(ista)-bef+1: tempoffs(ista)+aft) + ...
          1i * stackN(ista,tempoffs(ista)-bef+1: tempoffs(ista)+aft);
    
    % remove trend and mean
    STA = detrend(STA(:));
    
    for irot = 1: length(splitrotrad)
        %%% rotate the axes counter-clockwise by splitrotrad(irot), is equal to
        %%% rotate the complex record counter-clockwise by the same angle  
        STAfastslow = STA * exp(-1i * splitrotrad(irot));  % rotate the complex record counter-clockwise    
        STA1 = real(STAfastslow);
        STA2 = imag(STAfastslow);
        
        % remove trend and mean
        STA1 = detrend(STA1(:));
        STA2 = detrend(STA2(:));
        
%         if isequal(fam,'010')
%             off = round(8*sps/40);
%         else
%             off = round(12*sps/40);
%         end
        
        [coef, lag] = xcorr(STA1, STA2, 'coeff', splitoff);  % cc STA1 relative to STA2 
        [maxcoef, idx] = max(coef);
        lagsamp = lag(idx);     % if lagsamp >0, shift STA1 to left; <0, shift STA1 to right to align 
        coefmat(irot) = maxcoef;
        lagmat(irot) = lagsamp;
    end
    
    slowfastlag(ista) = lagmat(coefmat==max(coefmat));
    maxsplit(ista) = splitrot(coefmat==max(coefmat));
    slowfastlagmat(ista, :) = lagmat;
    slowfastcoefmat(ista, :) = coefmat;
    
    if abs(slowfastlag(ista))>=2   % accept the solution only if the split time is >=2 samples, from Peng etal (2015)
        if slowfastlag(ista) > 0   % i.e. STA1 is the slow component
            ROTS(ista, 1) =  slowfastlag(ista);
            ROTS(ista, 2) = maxsplit(ista);    % slow direction is the rotation angle
        else    % i.e. STA2 is the slow component, STA1 is the fast one
            ROTS(ista, 1) = -slowfastlag(ista);    % make sure this offset is always slow-fast, positive 
            if maxsplit(ista)< 90  
                ROTS(ista, 2) = maxsplit(ista)+90;  % +/- 90 degree to get slow direction
            else
                ROTS(ista, 2) = maxsplit(ista)-90;  % make sure that all direction is within [0,180)
            end
                
        end
    else    % <2, consider as no obvious splitting
        ROTS(ista, 1) = 0;
        ROTS(ista, 2) = 0;
    end
    
    figure   % plot individual correction results, fig 2~nsta+1
    
    STAfastslow = STA * exp(-1i * ROTS(ista, 2)*pi/180);     % rotate counter-clockwise
    STAslow = real(STAfastslow);    % now the real component must be the slow one
    STAfast = imag(STAfastslow);
    
    % remove trend and mean
    STAslow = detrend(STAslow(:));
    STAfast = detrend(STAfast(:));
    
    % plot
    plot(STAslow/max(STAslow) + 1, 'linewidth', 2); hold on
    plot(STAfast/max(STAfast) + 2, 'linewidth', 2); hold on
    
    %%% 1st column of PERMROTS == offset in samples between fast/slow direction
    STAslow(splitoff: seglen-splitoff) = STAslow(splitoff+ ROTS(ista, 1): seglen- splitoff+ ....
                                                 ROTS(ista, 1));    % shift the slow one to left (positive) to align
    STAsplitcorrected = (STAslow+ 1i * STAfast) * exp(1i* ROTS(ista, 2)*pi/180);   % rotate back
    
    % plot
    plot(STAslow/max(STAslow) + 3, 'linewidth', 2); hold on
    plot(STAfast/max(STAfast) + 4, 'linewidth', 2); hold on
    
    
    for irot = 1: length(polarrotrad)
        %%% Again, rotate the axes counter-clockwise by polarrotrad(irot), is equal to
        %%% rotate the complex record counter-clockwise by the same angle
        STAscrot = STAsplitcorrected * exp(-1i * polarrotrad(irot));   
        STAopt = real(STAscrot);    % assume the real component is optimal direction    
        STAort = imag(STAscrot);    
        
        % remove trend and mean
        STAopt = detrend(STAopt(:));
        STAort = detrend(STAort(:));
        engopt = sum(STAopt.^2);
        engort = sum(STAort.^2);
        engratio(irot) = engopt/engort;
    end
    maxratio(ista) = max(engratio);
    maxpolar = polarrot(engratio==max(engratio));
    optortratmat(ista, :) = engratio;
    
    %%% usually this maxpolar angle can make sure real axis is optimal
    %%% direction
    ROTS(ista, 3) = maxpolar;
    
    STAscrot = STAsplitcorrected* exp(-1i* maxpolar*pi/180);
    STAopt = real(STAscrot);
    STAort = imag(STAscrot);
    
    % remove trend and mean
    STAopt = detrend(STAopt(:));
    STAort = detrend(STAort(:));
    
    % plot
    plot(STAopt + 5, 'linewidth', 2); hold on
    plot(STAort + 6, 'linewidth', 2); hold on
    vertical_cursors
    
    STAoptmat(ista,:) = STAopt;
    STAortmat(ista,:) = STAort;
    
    %%% this is for saving the result in a longer segment
    STAlong = stackE(ista,:) + 1i * stackN(ista,:);
    STAfastslowlong = STAlong * exp(-1i * ROTS(ista, 2)*pi/180);
    STAslowlong = real(STAfastslowlong);
    STAfastlong = imag(STAfastslowlong);
    STAslowlong(splitoff: templen-splitoff) = STAslowlong(splitoff+ ROTS(ista, 1): ...
                                                          templen- splitoff+ ROTS(ista, 1));
    STAsplitcorrectedlong = (STAslowlong+ 1i * STAfastlong) * exp(1i* ROTS(ista, 2)*pi/180);
    STAscrotlong = STAsplitcorrectedlong* exp(-1i* ROTS(ista, 3)*pi/180);
    STAscrotlongmat(ista,:) = STAscrotlong;
    
end

%%%% just for plotting purpose 
for ista = 1: nsta
    ampmax = max(STAoptmat(ista, :));
    ampmin = min(STAoptmat(ista, :));
    norm = max(ampmax, -ampmin);
    dummy = real(STAscrotlongmat(ista, :))/norm;
%     %%% bandpass function encodes remove trend and mean
%     STAoptnorm(ista,:) = Bandpass(dummy, sps, lo, hi, npo, npa, 'butter');
    STAoptnormlong(ista,:) = dummy;
end

figure  % plot the longer normalized split corrected seismogram before shifting between different stations
for ista = 1:nsta
    plot(STAoptnormlong(ista,:) + 1* ista, 'linewidth', 1); hold on
    text(2/3*size(STAoptnormlong,2),ista,stas(ista,:));
end
%%%% just for plotting purpose 


for ista = 1: nsta
    ampmax = max(STAoptmat(ista, :));
    ampmin = min(STAoptmat(ista, :));
    norm = max(ampmax, -ampmin);
    dummy = STAoptmat(ista, :)/norm;
%     %%% bandpass function encodes remove trend and mean
%     STAoptnorm(ista,:) = Bandpass(dummy, sps, lo, hi, npo, npa, 'butter');
    STAoptnorm(ista,:) = dummy;
end

figure  % plot the shorter normalized split corrected seismogram before shifting between different stations
for ista = 1:nsta
    plot(STAoptnorm(ista, :) + 1* ista, 'linewidth', 2); hold on
    text(2/3*size(STAoptnorm,2),ista,stas(ista,:));
end

% mainsta = 4;    % index of main station, such as 4 for LZB
% staoff=8;
for ista = 1: nsta
    [coef, lag] = xcorr(STAoptnorm(ista,:), STAoptnorm(mainsta,:), 'coeff', staoff); % relative to LZB
    [maxcoef, idx] = max(coef);
    lagsamp = lag(idx);
    STAoptcoefmat(ista) = maxcoef;
    STAoptlagmat(ista) = lagsamp;
    ROTS(ista,4) = STAoptlagmat(ista)+tempoffs(ista)-tempoffs(mainsta);
end

figure      % plot normalized final optimal results, fig 9
for ista = 1: nsta
    plot(STAoptnorm(ista, staoff+STAoptlagmat(ista): seglen-staoff+STAoptlagmat(ista)) + ...
         1* ista, 'linewidth', 2); hold on
    text(2/3*size(STAoptnorm,2),ista,stas(ista,:));
end
vertical_cursors

% check if the alignment can be improved, lag should approach to 0
for ista = 1: nsta
    [coef, lag] = xcorr(STAoptnorm(ista, staoff+STAoptlagmat(ista): seglen-staoff+...
                        STAoptlagmat(ista)), STAoptnorm(mainsta, staoff+STAoptlagmat(mainsta)...
                        : seglen-staoff+STAoptlagmat(mainsta)), 'coeff');                         
    [maxcoef, idx] = max(coef);
    lagsamp = lag(idx);
    coefmat22(ista) = maxcoef;
    lagmat22(ista) = lagsamp;
end

%%% check if the orthogonal component is indeed very small, fig 10
figure('Position',[scrsz(3)/10 scrsz(4)/10 4*scrsz(3)/5 2*scrsz(4)/5]);  % plot individual optimal and orthogonal components
for ista = 1: nsta
    subplot(2,round(nsta/2),ista)
    plot(STAoptmat(ista, staoff+STAoptlagmat(ista): seglen-staoff+STAoptlagmat(ista)), 'k',...
         'linewidth', 2); hold on
    plot(STAortmat(ista, staoff+STAoptlagmat(ista): seglen-staoff+STAoptlagmat(ista)), 'r',...
          'linewidth', 2); hold on
    text(2/3*size(STAoptmat,2),0,stas(ista,:));
end


%%% check the overall one-day data, to see if the final operations are
%%% correct
for ista = 1: nsta
    if ROTS(ista,4) > -1    % or, >=0
        STAscrotlongmat(ista, 1:templen-ROTS(ista,4))=STAscrotlongmat(ista, ROTS(ista,4)+1:templen);
        STAscrotlongmat(ista, templen-ROTS(ista,4)+1:templen)=0;
    else
        STAscrotlongmat(ista, -ROTS(ista,4)+1:templen)=STAscrotlongmat(ista, 1:templen+ROTS(ista,4));
        STAscrotlongmat(ista, 1:-ROTS(ista,4))=0;
    end
    STAoptlongmat(ista,:)=real(STAscrotlongmat(ista,:));
    STAortlongmat(ista,:)=imag(STAscrotlongmat(ista,:));
end

% % fig 11
% figure('Position',[scrsz(3)/10 scrsz(4)/10 4*scrsz(3)/5 2*scrsz(4)/5]);  % plot individual optimal and orthogonal components
% for ista = 1: nsta
%     subplot(2,round(nsta/2),ista)
%     plot(STAoptlongmat(ista, :), 'k', 'linewidth', 2); hold on
%     plot(STAortlongmat(ista, :), 'r', 'linewidth', 2); hold on
%     text(2/3*size(STAortlongmat,2),0,stas(ista,:));
% end

% fig 12, plot the cc coef between fast and slow comp.
figure
for ista = 1: nsta
    plot(slowfastcoefmat(ista,:),'linewidth',1.5); hold on
end
legend(stas(:,:),'location','best');
text(1/3*size(slowfastcoefmat,2),0.5,'cc coef between fast and slow comp');

% fig 13, plot the enengy ratio between opt and ort comp.
figure
for ista = 1: nsta
    plot(optortratmat(ista,:),'linewidth',1.5); hold on
end
legend(stas(:,:),'location','best');
text(1/3*size(optortratmat,2),5,'enengy ratio between opt and ort comp');


% % save results for each family
% PREFIX = strcat(fam,'rot_para',num2str(lo),'-',num2str(hi),'_',num2str(bef),'_',num2str(aft),stas(mainsta,:));
% fid = fopen(strcat(datapath, '/split_chao/',PREFIX,'_',num2str(sps),'sps',num2str(nsta),'stas'),'w+');
% fprintf(fid,'%d %d %d %d \n',ROTS');
% fclose(fid);
    

keyboard

%     % convert the offset in 40 sps unit to the sps specified.
%     %%% 1st column of PERMROTS == offset in samples between fast and slow components
%     fastslowoff = round(POLROTS(idx,1)*sps/40);  % this offset is based on 40sps data, even for polaris stations
%     %%% 4th column of PERMROTS == offset in samples between the arrival times at different stations
%     STAsoff = round(POLROTS(idx,4)*sps/40);  %input offsets given at 40 sps.
%
%     %Rotate & split-correct
%     STA = STAEd+ 1i* STANd;
%     STAfastslow = STA* exp(-1i* POLROTS(idx,2));   %%% 2th column of PERMROTS, == rotation angle to get fast/slow direction
%     STAslow = real(STAfastslow);
%     STAfast = imag(STAfastslow);
%     len = length(STA);
%     off = round(10*sps/40);
%     %%% 1st column of PERMROTS == offset in samples between fast/slow direction
%     STAslow(off: len-off) = STAslow(off+ fastslowoff: len- off+ fastslowoff);
%     STAsplitcorrected = (STAslow+ 1i* STAfast)* exp(1i* POLROTS(idx,2));
%     STAscrot = STAsplitcorrected* exp(-1i* POLROTS(idx,3));
%
%     %Timeshift
%     if STAsoff > -1
%         STAscrot(1: tracelen- STAsoff) = STAscrot(STAsoff+ 1: tracelen);
%         STAscrot(tracelen- STAsoff+ 1: tracelen) = 0;
%     else
%         STAscrot(-STAsoff+ 1: tracelen) = STAscrot(1: tracelen+ STAsoff);
%         STAscrot(1: -STAsoff) = 0;
%     end
%
%     dataE = real(STAscrot);
%     dataN = imag(STAscrot);

