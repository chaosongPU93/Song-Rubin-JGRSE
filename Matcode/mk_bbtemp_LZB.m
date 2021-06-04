function [dstack, ccstack] = mk_bbtemp_LZB(fam,sps,templensec,ccmethod,plflag)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This is the script to generate the broadband LFE family template from Bostock's
% catalog, use the shear splitting and polarization parameters
% (theoratically) either from Yajun's or Chao's calculations for LZB trio.
% Saved for future use. 
%
% INPUT:
%   fam: family string
%   sps: desired sampling rate (original is uniformed to 40)
%   ccmethod: choice number to indicate which criteria to use in CC, 1 with threshold, 2 without
%   plflag: plotting flag whether to plot the templates (1) or not (0)
%
% OUTPUT:
%   dstack: direct stacks
%   ccstack: cross-coorelating corrected stacks
%
% USE data removed from station response
%
% Assume now you already have a catalog by hand, from others' LFE
% detections, or from your own rough detections, this script trys to make a
% template based on the catalog.
% TWO-STEP Precedure:
%   1. stack directly according to the timing of Bostock's catalog with 
%      cutted window. (stack at each station in each LFE family)
%   2. use the stacked template in 1 to cc with every cutted window to
%      make a time shift, then stack again to get a finer template
%
% UPDATE:
%   1. 2021/05/12.  insteading using a fixed before and after, now try a fixed,
%       then find the zero crossing, then update to center the window at zero
%       crossing
%
% By Chao Song, chaosong@princeton.edu
% First created date:   2019/08/30
% Last modified date:   2021/05/12
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% default value for easy debugging
defval('fam', '043');
defval('sps', 40);
defval('templensec', 60);
defval('ccmethod', 2);
defval('plflag', 1);

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
temppath = strcat(datapath, '/templates/LZBtrio/');

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

% number of used stations
nsta=size(stas,1);         %  number of stations

% ifam=2;
% for ifam = 1: length(nfampool)
% close all
% fam = nfampool(ifam, :)
fprintf('Family ID is %s \n', fam);

[timoffrot,~] = GetDays4Stack(fam);

% get permanent and polaris station rotation parameters
sft2=0;     % centroid shift of station 2
sft3=0;     % centroid shift of station 3
[PERMROTS, POLROTS] = GetRotsChao(fam,datapath,sft2,sft3);
reftime = POLROTS(5,4);     % reftime is the reference time at the 1st station,depends on the choice of 1st station
PERMROTS(:,4) = PERMROTS(:,4)-reftime;    % to make sure that 1st station is 0
POLROTS(:,4) = POLROTS(:,4)-reftime;

% convert angles to rads
PERMROTS(:,2:3)=pi*PERMROTS(:,2:3)/180.;     % convert 2-3 columns to rad
POLROTS(:,2:3)=pi*POLROTS(:,2:3)/180.;
% if isequal(fam,'243')
%     scaleseisms=[1.0 0.76 0.95 1.0 1.0 1.0 1.0];       % scaleseisms scales seismograms
% else
%     scaleseisms=[1.0 1.0 1.0 1.0 1.0 1.0 1.0];    
% end

%% stack to templates
%%% data parameters
% sps=40;     % samples per second
lenofday = 24*3600;
STAopt=zeros(nsta,sps * lenofday);
STAort=zeros(nsta,sps * lenofday);

%%% desired template parameters
% templensec = 60;
templen = templensec * sps;
before = 1*templen;
after = 1*templen-1;
switch fam
    case '001'
        offsec = 22.6;
    case '002'
        offsec = 22.7;
    case '006'
        offsec = 22.6;
    case '010'
        offsec = 22.625;
    case '017'
        offsec = 22.6;
    case '043'
        offsec = 22.75;
    case '047'
        offsec = 22.775;
    case '068'
        offsec = 22.55;
    case '099'
        offsec = 22.625;
    case '125'
        offsec = 22.6;
    case '141'
        offsec = 22.625;
    case '144'
        offsec = 23;
    case '147'
        offsec = 22.65;
    otherwise
        offsec = 0;
end
before = before - offsec*sps;
after = after + offsec*sps;

extra = 6*sps;    % use extra more samples than the original window.
stacktmp = zeros(nsta, 2*templen+ 2*extra);
stacktmport = zeros(nsta, 2*templen+ 2*extra);
% stacktmp = zeros(nsta, templen+ 2*extra);
% stacktmport = zeros(nsta, templen+ 2*extra);

nday = size(timoffrot, 1)

if nday == 0
    disp("No day with detections found in this family of BOSTOCK's catalog");
    return
end

% get the all catalog LFEs in that family
bostname = ('/BOSTOCK/total_mag_detect_0000_cull_NEW.txt');
catalog = load(strcat(workpath, bostname));
famnum = str2double(fam);
dateall = catalog(famnum == catalog(:, 1), :);


%% STEP 1: In a selected family, loop for every day, read data & catalog, stack
%%% loop for day
nstack = 0;
nLFE = 0;
fprintf('Part 1: Raw Template \n');

for id = 1: nday
    %     id = 1;
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
    prename = [direc, YEAR,'.',JDAY,'.00.00.00.0000.CN'];
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
            fname = strcat(prename,'.',PERMSTA(idx,1:3),'..BHE.D.SAC');
            if isfile(fname)    % if have the data file
                % 1. this is for data without removing station response
                %   [opt,ort,timeperm]=readperm_1daynofilterv2(datafnm, PERMSTA, PERMROTS, idx, sps, fact);
                % 2. this is for data with no response
                [opt, ort, timeperm] = readperm_nofilterv2(prename, ...
                    PERMSTA, PERMROTS, idx, sps, fact);
            else
                fileflag = 0;   % change the file flag to 0, meaning abnormal
                fprintf('No data for station %s in day %s / %s, this day will be omitted. \n',...
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
            fname = strcat(prename,'.',POLSTA(idx,1:4),'..HHE.D.SAC');
            if isfile(fname)    % if have the data file
                %                 % 1. this is for data without removing station response
                %                 [opt,ort,nzeros]=readpols(prename,POLSTA,POLROTS,idx,sps,lo,hi,npo,npa,fact,nwin,winlen,winoff,igstart);
                % 2. this is for data with no response
                [opt, ort, ~] = readpola_nofilterv2(prename, ...
                    POLSTA, POLROTS, idx, sps, fact);
            else
                fileflag = 0;   % change the file flag to 0, meaning abnormal
                fprintf('No data for station %s in day %s / %s, this day will be omitted. \n',...
                        POLSTA(idx,1:4), YEAR, JDAY);
                break   % break the entire station loop
            end
            
        end
        
        STAopt(ista, :) = opt;      % different from PGC trio since lack of 'scaleseisms' info
        STAort(ista, :) = ort;
    end
    
    if fileflag == 0    % means there are missing files
        fprintf('Day %s / %s will be omitted because of missing files. \n', YEAR, JDAY);
        continue    % continue to the next day
    end
    
    len = size(STAopt,2);
    for n=1: size(bostsamp, 1)
        % set the timing be at the center of the window
        if (bostsamp(n)+ after+ extra <= len) && ...
                (bostsamp(n)- before- extra >= 1)
            
            %%% it is better to use a longer window than templen to make
            %%% sure the shifting in step2 would preserve the frequency
            %%% content as much as possible
            windata = STAopt(:, bostsamp(n)- before- extra: ...
                bostsamp(n)+ after+ extra);
            % remove mean and linear trend
            tmp = detrend(windata');
            stacktmp = stacktmp+ tmp';
            
            windata2 = STAort(:, bostsamp(n)- before- extra: ...
                bostsamp(n)+ after+ extra);
            tmp2 = detrend(windata2');
            stacktmport = stacktmport+ tmp2';
            
            nstack = nstack + 1;
            
            % datamat is 3D, for storing all qualified wins for stacking
            % size is (nsta, nstack, templen)
            datamat(:, nstack, :) = STAopt(:, bostsamp(n)- before- extra: ...
                bostsamp(n)+ after+ extra);
        end
    end
    
end     % loop for days

%% find the location of zero-crossing, and adjust the window to be centered at zero-crossing
%this step is replaced by direct correction in time, original step can be found in 'testmk_bbtemp'
%now the indzero is basically the center of current segment, no need to find it any more
% [~,indmax] = max(stacktmp(:, :),[],2);
% indmaxmed = median(indmax);
% [~,indmin] = min(stacktmp(:, :),[],2);
% indminmed = median(indmin);
% indzero = round((indminmed+indmaxmed)/2);
indzero = size(stacktmp, 2)/2;
indst = indzero - templen/2 -extra;
inded = indzero + templen/2 +extra-1;

windata = [];
stackex = zeros(nsta, templen+ 2*extra);
for i = 1: nstack
    windata(:,:) = datamat(:, i, indst: inded);
    % remove mean and linear trend
    tmp = detrend(windata');
    stackex = stackex+ tmp';
    
end


%% Averaging
%%% cut off extra points, stack should be templen long
% dstack = stacktmp(:, 1+extra: end-extra);
dstack = stackex(:, 1+extra: end-extra);
% dstackort = stackexort(:, 1+extra: end-extra);

%%% averaging
dstack = dstack/ nstack;
% dstackort = dstackort/nstack;

%%% detrend
tmp = detrend(dstack');
dstack = tmp';
% tmp2 = detrend(dstackort');
% dstackort = tmp2';


%% normalization
mid = templen/2;
ampmax = max(dstack(:, mid-sps+1: mid+sps),[],2);
ampmin = min(dstack(:, mid-sps+1: mid+sps),[],2);   % /2
for ista = 1: nsta
    if ampmax(ista) >= -ampmin(ista)
        fprintf('At station %s positive is bigger \n', strtrim(stas(ista, :)));
    else
        fprintf('At station %s negative is bigger \n', strtrim(stas(ista, :)));
    end
end
norm = max(ampmax, -ampmin);
dstackno = dstack./ norm;
% dstackortno = dstackort./ norm;


%% plotting if desired
if plflag
    %%% plot the 1-step stacked template & write into file    
    %%% figure 1
    figure('Position',[0.1*wid 0.1*hite 0.3*wid 0.3*hite]);
    for ista=1: nsta
        plot(dstack(ista, :) + 0.1* ista, 'linewidth', 2); hold on
        title('direct stack');
    end
    box on
    grid on
    
    %%% plot the normalized 1-step template
    %%% figure 2
    figure('Position',[0.1*wid 0.2*hite 0.3*wid 0.3*hite]);
    for ista = 1: nsta
        plot((1:10*sps+1)/sps,dstackno(ista, mid-5*sps: mid+5*sps) + 1* ista, 'linewidth', 1);hold on
        %     plot(stackno(ista, mid-5*sps: mid+5*sps) + 1* ista, 'linewidth', 1,'color',...
        %          color(ista)); hold on
        text(1,1*ista,stas(ista,:),'fontsize',12);
    end
    title('normalized direct stack');
    xlabel('Time (s)');
    % xlabel('Samples (on 40 sps)');
    ylabel('Normalized Amp.');
    box on
    grid on
    
    %%% figure 3
    figure('Position',[0.1*wid 0.3*hite 0.3*wid 0.3*hite]);
    for ista = 1: nsta
        plot(dstack(ista, mid-5*sps: mid+5*sps) + 0.1* ista, 'linewidth', 2); hold on
        title('direct stack-short window');
    end
    box on
    grid on
    
end


%% PART 2: Multi-CC the 1-step template and each window to get a finer template
% ccmethod = 1;   % which method is using, 1 with threshold; 2 without

if ccmethod == 1    % start for CCMETHOD choice
    %%%%%%%%%%%%%%%% WAY 1, stack only the ones that pass the cc and offset thresholds %%%%%%%%%%%%%%%%%
    %%% the previous 'extra' is tO avoid adding zeros to the array due to shifting,
    nstack2 = 0;   % here nstack2 can vary w/ different stations
    % according to template, the duration of dipole seems to be ~0.5 sec (10
    % samples), could even as small as 5 samples
    offsetmax = sps/4;
    % 0.25/ 0.3 could be a reasonable cc coef threshold
    coefmin = 0.25;     % should vary according to different trio & fam
    ccstack = zeros(nsta, templen);
    %%% SEE NOTES for proper freq. band choices
    lo = 2;     % empirical 
    hi = 8;
    npo = 2;
    npa = 2;
    coefmat = zeros(nsta, nstack);
    lagmat = zeros(nsta, nstack);
    fprintf('Part 2: Precise Template \n');
    % figure
    for ista = 1: nsta    % loop for station
        %     ista = 1;
%         dummy(:,:) = datamat(ista,:,:);
        dummy(:,:) = datamat(ista,:, indst: inded);
        windata = dummy';
        dummy2 = stackex(ista, :)/ nstack;
        template = dummy2';
        for istack = 1: nstack
            fwindata(:,istack) = Bandpass(windata(:,istack), sps,lo,hi,npo,npa,'butter');  % also contains detrend
        end
        ftemplate = Bandpass(template, sps,lo,hi,npo,npa,'butter');
        %%% to SEE how bad is our windowed data before/after filtering
        %         subplot(4,1,1)
        %         plot(windata(templen/2-4*sps+1:templen/2+4*sps), 'b');
        %         subplot(4,1,2)
        %         plot(fwindata(templen/2-4*sps+1:templen/2+4*sps), 'r');
        %         subplot(4,1,3)
        %         plot(template(templen/2-4*sps+1:templen/2+4*sps), 'b');
        %         subplot(4,1,4)
        %         plot(ftemplate(templen/2-4*sps+1:templen/2+4*sps), 'r');
        %         pause(2)
        %%%
        
        %%% e.g. if lag = argmax(xcorr(A, B)) = 4, A needs to shift 4 to
        %%% the left to align with B, (i.e. A-B=4), B is the reference station
        win = fwindata(templen/2+extra-2*sps+1:templen/2+extra+2*sps,:);  % the length is 4*sps
        tem = ftemplate(templen/2+extra-2*sps+1:templen/2+extra+2*sps,:);
        win = detrend(win);
        tem = detrend(tem);
        for istack = 1: nstack
            [coef, lag] = xcorr(win(:,istack), tem, offsetmax, 'coeff');
            [maxcoef, idx] = max(coef);
            lagsamp = lag(idx);
            coefmat(ista, istack) = maxcoef;
            lagmat(ista, istack) = lagsamp;
        end
        
    end   % loop for station
            
    windata = [];
    avecoef = sum(coefmat,1);    
    for istack = 1: nstack
        
        if (abs(lagmat(:, istack)) <= offsetmax) && (avecoef(istack) >= coefmin) ...
                && (coefmat(:,istack) >= 2/3*coefmin)
%             windata = datamat(:,istack,:);
            windata = datamat(:,istack, indst: inded);
            newdata = windata(:,1+lagmat(ista, istack)+extra: ...
                                lagmat(ista, istack)+extra+templen);
            nstack2 = nstack2+ 1;
            winlag(:, nstack2) = lagmat(:, istack);
            tmp = detrend(newdata');
            ccstack = ccstack+ tmp';
        end
    end
            %%% This was the test for long/short cc comparison ,now in a
            %%% separate test file
            %         wincut = detrend(windata(templen+extra-offsetmax+1: templen+extra+offsetmax));
            %         temcut = detrend(template(templen+extra-offsetmax+1: templen+extra+offsetmax));
            %         [coef, lag] = xcorr(wincut, temcut, 'coeff');
            %         [maxcoef, idx] = max(coef);
            %         lagsamp = lag(idx);
            %
            %         if (abs(lagsamp) < offsetmax)  && (maxcoef >= coefmin)
            %             newdata(:) = windata(1+lagsamp+extra: lagsamp+extra+templen);
            %             nstack2(ista) = nstack2(ista)+ 1;
            %             winlag(ista, nstack2(ista)) = lagsamp;
            %             tmp = detrend(newdata');
            %             stack2(ista, :) = stack2(ista, :)+ tmp';
            %         end
            %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%% WAY 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
elseif ccmethod == 2
    %%%%%%%%%%%%%%%% WAY 2, stack all, only use cc to correct time %%%%%%%%%%%%%%%%%
    %%% the previous 'extra' is tO avoid adding zeros to the array due to shifting,
    nstack2 = nstack;
    % according to template, the duration of dipole seems to be ~0.5 sec (10
    % samples), could even as small as 5 samples
    offsetmax = sps/4;
    ccstack = zeros(nsta, templen);
    %%% SEE NOTES for proper freq. band choices
    lo = 2;     % empirical 
    hi = 8;
    npo = 2;
    npa = 2;
    coefmat = zeros(nsta, nstack);
    lagmat = zeros(nsta, nstack);
    fprintf('Part 2: Precise Template \n');
    % figure
    for ista = 1: nsta    % loop for station
        %     ista = 1;
%         dummy(:,:) = datamat(ista,:,:);
        dummy(:,:) = datamat(ista,:, indst: inded);
        windata = dummy';
        dummy2 = stackex(ista, :)/ nstack;
        template = dummy2';
        %             tic
        for istack = 1: nstack
            fwindata(:,istack) = Bandpass(windata(:,istack), sps, lo, hi, npo, npa, 'butter');  % also contains detrend
        end
        %             toc
        ftemplate = Bandpass(template, sps, lo, hi, npo, npa, 'butter');
        %%% to SEE how bad is our windowed data before/after filtering
        %         subplot(4,1,1)
        %         plot(windata(templen/2-4*sps+1:templen/2+4*sps), 'b');
        %         subplot(4,1,2)
        %         plot(fwindata(templen/2-4*sps+1:templen/2+4*sps), 'r');
        %         subplot(4,1,3)
        %         plot(template(templen/2-4*sps+1:templen/2+4*sps), 'b');
        %         subplot(4,1,4)
        %         plot(ftemplate(templen/2-4*sps+1:templen/2+4*sps), 'r');
        %         pause(2)
        %%%
        
        %%% e.g. if lag = argmax(xcorr(A, B)) = 4, A needs to shift 4 to
        %%% the left to align with B, (i.e. A-B=4), B is the reference station
        win = fwindata(templen/2+extra-2*sps+1:templen/2+extra+2*sps,:);  % the length is 4*sps
        tem = ftemplate(templen/2+extra-2*sps+1:templen/2+extra+2*sps,:);
        win = detrend(win);
        tem = detrend(tem);
        for istack = 1: nstack
            [coef, lag] = xcorr(win(:,istack), tem, offsetmax, 'coeff');
            [maxcoef, idx] = max(coef);
            lagsamp = lag(idx);
            coefmat(ista, istack) = maxcoef;
            lagmat(ista, istack) = lagsamp;
            
            newdata(:) = windata(1+lagsamp+extra: lagsamp+extra+templen,istack);
            winlag(ista, istack) = lagsamp;
            tmp = detrend(newdata');
            ccstack(ista, :) = ccstack(ista, :)+ tmp';
        end
            
            %%% This was the test for long/short cc comparison ,now in a
            %%% separate test file
            %         wincut = detrend(windata(templen+extra-offsetmax+1: templen+extra+offsetmax));
            %         temcut = detrend(template(templen+extra-offsetmax+1: templen+extra+offsetmax));
            %         [coef, lag] = xcorr(wincut, temcut, 'coeff');
            %         [maxcoef, idx] = max(coef);
            %         lagsamp = lag(idx);
            %
            %         if (abs(lagsamp) < offsetmax)  && (maxcoef >= coefmin)
            %             newdata(:) = windata(1+lagsamp+extra: lagsamp+extra+templen);
            %             nstack2(ista) = nstack2(ista)+ 1;
            %             winlag(ista, nstack2(ista)) = lagsamp;
            %             tmp = detrend(newdata');
            %             stack2(ista, :) = stack2(ista, :)+ tmp';
            %         end
            %%%
                    
    end     % loop for station
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%% WAY 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
end     % end for CCMETHOD choice

%% Averaging and normalization 
%%% averaging
for ista = 1: nsta
    ccstack(ista, :) = ccstack(ista, :)/ nstack2;
end

%%% detrend
tmp = detrend(ccstack');
ccstack = tmp';

%%% normalization
ampmax = max(ccstack(:, mid-sps+1: mid+sps),[],2);    % /2
ampmin = min(ccstack(:, mid-sps+1: mid+sps),[],2);
for ista = 1: nsta
    if ampmax(ista) >= -ampmin(ista)
        fprintf('At station %s positive is bigger \n', strtrim(stas(ista, :)));
    else
        fprintf('At station %s negative is bigger \n', strtrim(stas(ista, :)));
    end
end
norm = max(ampmax, -ampmin);
ccstackno = ccstack./ norm;

%% plotting if desired
if plflag
    %%% plot the 2-step stacked template 
    %%% figure 1
    figure('Position',[0.1*wid 0.1*hite 0.3*wid 0.3*hite]);
    for ista=1: nsta
        plot(ccstack(ista, :) + 0.1* ista, 'linewidth', 2); hold on
        title('CC stack');
    end
    box on
    grid on
    
    %%% plot the normalized 2-step template
    %%% figure 2
    figure('Position',[0.1*wid 0.2*hite 0.3*wid 0.3*hite]);
    for ista = 1: nsta
        plot((1:10*sps+1)/sps,ccstackno(ista, mid-5*sps: mid+5*sps) + 1* ista, 'linewidth', 1);hold on
        %     plot(stack2no(ista, mid-5*sps: mid+5*sps) + 1* ista, 'linewidth', 1,'color',...
        %          color(ista)); hold on
        text(1,1*ista,stas(ista,:),'fontsize',12);
    end
    title('normalized CC stack');
    xlabel('Time (s)');
    % xlabel('Samples (on 40 sps)');
    ylabel('Normalized Amp.');
    box on
    grid on
    
    if ccmethod == 1
        %%% figure 3
        figure
        supertit(gca,'max cc coef vs sample shift of all traces');
        for ista = 1: nsta
            subplot(nsta,1,ista)
            scatter(lagmat(ista, :), coefmat(ista, :)); hold on
            plot([-offsetmax, offsetmax], [coefmin, coefmin], 'r-', 'linewidth', 2);
        end
        box on
        grid on
    end
    
    %%% figure 4
    figure
    supertit(gca,'sample shift of preserved traces passing threshold');
    for ista = 1: nsta
        subplot(nsta,1,ista)
        scatter(1: nstack2, winlag(ista, 1:nstack2));
    end
    box on
    grid on

end


%% calculate the SNR
for ista = 1: nsta
    s = dstack(ista, mid-0.5*sps+1: mid+0.5*sps);
    Es(ista) = sum(s.^2);
    np1 = dstack(ista, mid-1.5*sps+1: mid-0.5*sps);
    np2 = dstack(ista, mid+0.5*sps+1: mid+1.5*sps);
    Enp1(ista) = sum(np1.^2);
    Enp2(ista) = sum(np2.^2);
    En(ista) = (Enp1(ista)+Enp2(ista))/2;
    snr(ista) = Es(ista)/En(ista);
end

for ista = 1: nsta
    s = ccstack(ista, mid-0.5*sps+1: mid+0.5*sps);
    Es2(ista) = sum(s.^2);
    np1 = ccstack(ista, mid-1.5*sps+1: mid-0.5*sps);
    np2 = ccstack(ista, mid+0.5*sps+1: mid+1.5*sps);
    Enp1(ista) = sum(np1.^2);
    Enp2(ista) = sum(np2.^2);
    En2(ista) = (Enp1(ista)+Enp2(ista))/2;
    snr2(ista) = Es2(ista)/En2(ista);
end

%% plot if desired
if plflag
    %%% plot the 2 kinds of template together
    f.fig = figure(111);
    set(f.fig,'Position',[scrsz(3)/10 scrsz(4)/10 3*scrsz(3)/5 4*scrsz(4)/5]);    
    for ista = 1: nsta
        subplot(nsta,2,2*ista-1);
        plot(dstack(ista, mid-5*sps: mid+5*sps), 'k', 'linewidth', 1); hold on
        plot(ccstack(ista, mid-5*sps: mid+5*sps), 'r', 'linewidth', 1);
        title(strcat('Stacked template with/without CC at station', {' '}, strtrim(stas(ista, :))));
        legend('Direct stack', 'CC stack');
        xlabel('Samples');
        ylabel('Amplitude');
        box on
        grid on
        hold off
    end
    
    for ista = 1: nsta
        subplot(nsta,2,2*ista);
        plot(dstackno(ista, mid-5*sps: mid+5*sps), 'k', 'linewidth', 1); hold on
        plot(ccstackno(ista, mid-5*sps: mid+5*sps), 'r', 'linewidth', 1);
        title(strcat('Normalized stacked template with/without CC at station', {' '}, ...
              strtrim(stas(ista, :))));
        legend('Direct stack', 'CC stack');
        xlabel('Samples');
        ylabel('Normalized amplitude');
        box on
        grid on
        hold off
    end
    
    print(f.fig,'-depsc2',strcat(temppath,'/',fam,'.broadbandstack.eps'));
    
    close(f.fig);
end

% keyboard


















