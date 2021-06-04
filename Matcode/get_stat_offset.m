%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% This is the script to get the templates at different stations in a 
% particular family and the offset between them (4th col. in PERMROTS and 
% POLROTS), assuming that all splitting and polarization angles obtained
% by Yajun are correct (first 3 col. needed in PERMROTS and POLROTS).
% 
% Features:
%   1. Based on the corrections by Yajun, assuming they are all right.
%   2. Obtain the templates using corrections at each station in one family
%   3. Obtain time offset between each template (4th col. needed by in 
%      PERMROTS and POLROTS). Then align them relative to the 1st station
%
%
% By Chao Song, chaosong@princeton.edu
% First created date:   2019/04/19
% Last modified date:   2019/04/19
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Initialization
format short e   % Set the format to 5-digit floating point
clear
close all

set(0,'DefaultFigureVisible','on');
% set(0,'DefaultFigureVisible','off');   % switch to show the plots or not
scrsz=get(0,'ScreenSize');   % for example, scrsz=[1,1,1680,1050]
wid=scrsz(3);       % width
hite=scrsz(4);       % height

% set path
workpath = getenv('ALLAN');
temppath = strcat(workpath, '/templates/');
if ~exist(temppath, 'dir')
  mkdir(temppath);
end
% datapath = workpath;
datapath = strcat(workpath,'/data-no-resp');
rstpath = strcat(datapath, '/PGCtrio');

ofampool = ['002';   
           '013';
           '025';
           '028';
           '056';
           '084';
           '099';
           '115';
           '125';
           '147';
           '149'];

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
       
% fam='010'     % family number
% for ifam = 1: length(nfampool)
    ifam=1;
    close all
    
    fam = nfampool(ifam, :)
    
    if isequal(fam,'002')  %DON'T FORGET additional family-specific delcarations around line 156
        
        %%% permanent stations rotation angles, NAMES can be found in Rubin &
        %%% Armbruster 2013, Armbruster et al., 2014
        
        %%%%%%%% Meanning of each column in PERMROTS/POLROTS %%%%%%%%
        %%% 1. offset in samples between fast/slow direction for SAME station,
        %      given at 40 sps
        %%% 2. rotation angle to get fast/slow direction for SAME station
        %%% 3. rotation angle to maximize the energy/particle
        %%%    motion/polarization, for SAME station
        %%% 4. offset in samples between the arrival times at different
        %%%    stations, given at 40 sps
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        PERMROTS=[0 90 32 0;  %PGC , Yajun's "0" changed to 90.  Don't think that matters, if no splitting.
            0 90 54 0]; %LZB
        %%% POLARIS stations rotation angles
        POLROTS=[6 85 33 0;  %SSIB from Yajun
            0 90 39 0;  %SILB
            0 90  7 0;  %KLNB
            4 70 48 0;  %MGCB
            4 75 38 0]; %TWKB
        
        % generate unique dates matrix that in that family
        timoffrot = Readbostock(fam);
        ind=[];
        timoffrot(ind, :)=[];
        
        tempoffs=[1221 1211 1297 1231 1207 1186]; %these are the zero crossings
        
    elseif isequal(fam, '043')  % use 043 instead of 013, but with error
        
        PERMROTS=[0  0 29  0;  %PGC
            6 95 41  0]; %LZB
        POLROTS =[0  0  0  0;  %SSIB
            0  0  0  0;  %SILB
            0  0 20  0;  %KLNB
            2 65 36  0;  %MGCB
            3 65 12  0]; %TWKB
        
        % generate unique dates matrix that in that family
        timoffrot = Readbostock(fam);
        ind=[];
        timoffrot(ind, :)=[];
        
        tempoffs=[1195 1289 1374 1349 1210 1208]; %these are the zero crossings
        
    elseif isequal(fam, '141')  % use 141 instead of 025, but with error
        
        PERMROTS=[0  0  0  0;  %PGC
            3 80 15  0]; %LZB
        POLROTS =[0  0  0  0;  %SSIB
            0  0  0  0;  %SILB
            0  0 27  0;  %KLNB
            0  0 41  0;  %MGCB
            3 45 13  0]; %TWKB
        
        % generate unique dates matrix that in that family
        timoffrot = Readbostock(fam);
        ind=[];
        timoffrot(ind, :)=[];
        
        tempoffs=[1182 1318 1398 1391 1204 1214]; %these are the zero crossings
        
    elseif isequal(fam, '047')  % use 047 instead of 028
        
        PERMROTS=[0  0 36  0;  %PGC
            5 60 39  0]; %LZB
        POLROTS =[5 90 44  0;  %SSIB
            0  0 42  0;  %SILB
            0  0  4  0;  %KLNB
            3 80 48  0;  %MGCB
            4 70 31  0]; %TWKB
        
        % generate unique dates matrix that in that family
        timoffrot = Readbostock(fam);
        ind=[];
        timoffrot(ind, :)=[];
        
        tempoffs=[1212 1240 1331 1274 04 1196]; %these are the zero crossings
        
    elseif isequal(fam, '010')  % use 010 instead of 056
        
        PERMROTS=[0  0  0  0;  %PGC
            2 55  2  0]; %LZB
        POLROTS =[0  0  0  0;  %SSIB
            0  0  0  0;  %SILB
            0  0  0  0;  %KLNB
            2 95 51  0;  %MGCB
            5 70 27  0]; %TWKB
        
        % generate unique dates matrix that in that family
        timoffrot = Readbostock(fam);
        ind=[7 10];
        timoffrot(ind, :)=[];
        
        %%% PROBLEMATIC
        tempoffs=[1179 1278 1391 1296 1204 1201]; %these are the zero crossings
        
    elseif isequal(fam, '144')  % use 144 instead of 084
        
        PERMROTS=[3 85 52  0;  %PGC
            4 90 29  0]; %LZB
        POLROTS =[0  0  0  0;  %SSIB
            0  0  0  0;  %SILB
            0  0  0  0;  %KLNB
            3 90 50  0;  %MGCB
            4 80 30  0]; %TWKB
        
        % generate unique dates matrix that in that family
        timoffrot = Readbostock(fam);
        ind=[];
        timoffrot(ind, :)=[];
        
        tempoffs=[1191 1341 1413 1419 1217 1226]; %these are the zero crossings
        
        
    elseif isequal(fam, '099')
        
        PERMROTS=[0  0  0  0;  %PGC
            2 80 13  0]; %LZB
        POLROTS =[0  0  0  0;  %SSIB
            0  0  0  0;  %SILB
            0  0  0  0;  %KLNB
            0  0 35  0;  %MGCB
            2 65 19  0]; %TWKB
        
        % generate unique dates matrix that in that family
        timoffrot = Readbostock(fam);
        ind=[8];
        timoffrot(ind, :)=[];
        
        tempoffs=[1179 1315 1400 1378 1205 1212]; %these are the zero crossings
        
    elseif isequal(fam, '068')  % use 068 instead of 115
        
        PERMROTS=[0  0 48  0;  %PGC
            8 65  6  0]; %LZB
        POLROTS =[0  0  0  0;  %SSIB
            0  0  0  0;  %SILB
            2 50 25  0;  %KLNB
            5 65 41  0;  %MGCB
            4 60 14  0]; %TWKB
        
        % generate unique dates matrix that in that family
        timoffrot = Readbostock(fam);
        ind=[];
        timoffrot(ind, :)=[];
        
        tempoffs=[1193 1281 1347 1352 1201 1201]; %these are the zero crossings
        
    elseif isequal(fam, '125')
        
        PERMROTS=[0  0 40  0;  %PGC
            8 75 26  0]; %LZB
        POLROTS =[0  0  0  0;  %SSIB
            0  0  0  0;  %SILB
            0  0 10  0;  %KLNB
            4 95 60  0;  %MGCB
            6 55 16  0]; %TWKB
        
        % generate unique dates matrix that in that family
        timoffrot = Readbostock(fam);
        ind=[];
        timoffrot(ind, :)=[];
        
        tempoffs=[1200 1262 1341 1316 1203 1198]; %these are the zero crossings for 002: PGC,SSIB,SILB.
        
    elseif isequal(fam, '147')
        
        PERMROTS=[0  0 39  0;  %PGC
            6 85 36  0]; %LZB
        POLROTS =[0  0  0  0;  %SSIB
            0  0  0  0;  %SILB
            0  0 24  0;  %KLNB
            3 80 51  0;  %MGCB
            4 75 27  0]; %TWKB
        
        % generate unique dates matrix that in that family
        timoffrot = Readbostock(fam);
        ind=[2 4];
        timoffrot(ind, :)=[];
        
        tempoffs=[1191 1309 1378 1383 1205 1210]; %these are the zero crossings for 002: PGC,SSIB,SILB.
        
    elseif isequal(fam, '017')    % use 017 instead of 149
        
        PERMROTS=[0   0 39  0;  %PGC
            6  90 35  0]; %LZB
        POLROTS =[0   0  0  0;  %SSIB
            6 115 57  0;  %SILB
            2  55 21  0;  %KLNB
            3  85 54  0;  %MGCB
            4  70 21  0]; %TWKB
        
        % generate unique dates matrix that in that family
        timoffrot = Readbostock(fam);
        ind=[2 7];
        timoffrot(ind, :)=[];
        
        tempoffs=[1186 1311 1379 1389 1203 1209]; %these are the zero crossings for 002: PGC,SSIB,SILB.
        
    end
    
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
    
    % convert angles to rads
    PERMROTS(:,2:3)=pi*PERMROTS(:,2:3)/180.;     % convert 2-3 columns to rad
    POLROTS(:,2:3)=pi*POLROTS(:,2:3)/180.;
    
    %%% data parameters
    sps=40;     % samples per second
    lenofday = 24*3600;
    STAopt=zeros(nsta,sps * lenofday);
    STAort=zeros(nsta,sps * lenofday);
    
    %%% desired template parameters
    templensec = 60;
    templen = templensec * sps;
    % before = templen/2-1;
    % after = templen/2;
    before = templen/8;
    after = 7*templen/8-1;
    extra = 6*sps;    % use extra more samples than the original window.
    stackex = zeros(nsta, templen+ 2*extra);
    stackexort = zeros(nsta, templen+ 2*extra);
    
    
    nday = size(timoffrot, 1)
    % timoffrot=timoffrot(5:9, :);
    nday = size(timoffrot, 1);
    
    if nday == 0
        disp("No day with detections found in this family of BOSTOCK's catalog");
    else
        
        % get the all catalog LFEs in that family
        bostname = ('/BOSTOCK/total_mag_detect_0000_cull_NEW.txt');
        catalog = load(strcat(workpath, bostname));
        famnum = str2double(fam);
        dateall = catalog(famnum == catalog(:, 1), :);
        
        %%
        %%% loop for day
        nstack = 0;
        nLFE = 0;
        fprintf('Part 1: Raw Template \n');
        
        for id = 1: nday
            %     nd = 1;
            close all
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
            bostsamp = round(bostsec*40);
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
                        [opt, ort, timeperm] = readperm_nofilterv2(datafnm, ...
                            PERMSTA, PERMROTS, idx, sps, fact); 
                    else
                        fileflag = 0;   % change the file flag to 0, meaning abnormal
                        fprintf('No data for station %s in day %s %s, this day will be omitted. \n',PERMSTA(idx,1:3), YEAR, JDAY);
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
                        [opt, ort, timepola] = readpola_nofilterv2(datafnm, ...
                            POLSTA, POLROTS, idx, sps, fact); 
                    else
                        fileflag = 0;   % change the file flag to 0, meaning abnormal
                        fprintf('No data for station %s in day %s / %s. \n',POLSTA(idx,1:4), YEAR, JDAY);
                        break   % break the entire station loop
                    end

                end
                
                STAopt(ista, :) = opt;
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
                    stackex = stackex+ tmp';
                    
                    windata2 = STAort(:, bostsamp(n)- before- extra: ...
                        bostsamp(n)+ after+ extra);
                    tmp2 = detrend(windata2');
                    stackexort = stackexort+ tmp2';
                    
                    nstack = nstack + 1;
                    
                    % datamat is 3D, for storing all qualified wins for stacking
                    % size is (nsta, nstack, templen)
                    datamat(:, nstack, :) = STAopt(:, bostsamp(n)- before- extra: ...
                        bostsamp(n)+ after+ extra);
                end
            end
            
        end     % loop for days
        
        %% Averaging and and plotting
        
        %%% cut off extra points, stack should be templen long
        stack = stackex(:, 1+extra: end-extra);
        stackort = stackexort(:, 1+extra: end-extra);
        
        %%% averaging
        stack = stack/ nstack;
        stackort = stackort/nstack;
        
        %%% detrend
        tmp = detrend(stack');
        stack = tmp';
        tmp2 = detrend(stackort');
        stackort = tmp2';
        
        %%% plot the 1-step stacked template & write into file
        %%% plot
        
        %%% figure 1
        figure('Position',[0.1*wid 0.1*hite 0.3*wid 0.3*hite]);
        for ista=1: nsta
            plot(stack(ista, :) + 0.1* ista, 'linewidth', 2); hold on
            title('direct stack');
        end
        box on
        
        
        %% normalization
        stackno = zeros(nsta, templen);
        stackortno = zeros(nsta, templen);
        tempoffs_sort = sort(tempoffs);
        
        if mod(nsta,2) == 0
            mid = floor((tempoffs_sort(nsta/2)+tempoffs_sort(nsta/2+1))/2);
        else
            mid = tempoffs_sort((nsta+1)/2);
        end
        
        for ista = 1: nsta
            
            ampmax = max(stack(ista, mid-5*sps+1: mid+5*sps));
            ampmin = min(stack(ista, mid-5*sps+1: mid+5*sps));
            if ampmax >= -ampmin
                fprintf('At station %s positive is bigger \n', strtrim(stas(ista, :)));
            else
                fprintf('At station %s negative is bigger \n', strtrim(stas(ista, :)));
            end
            norm = max(ampmax, -ampmin);
            stackno(ista, :) = stack(ista, :)/ norm;
            stackortno(ista, :) = stackort(ista, :)/ norm;
        end
        
        %%% plot the normalized 1-step template
        %%% figure 2
        figure('Position',[0.1*wid 0.2*hite 0.3*wid 0.3*hite]);
        for ista = 1: nsta
            plot(stackno(ista, :) + 2* ista, 'linewidth', 2); hold on
            title('normalized direct stack');
        end
        box on
        
        %%% figure 3
        figure('Position',[0.1*wid 0.3*hite 0.3*wid 0.3*hite]);
        for ista = 1: nsta
            plot(stack(ista, mid-5*sps: mid+5*sps) + 0.1* ista, 'linewidth', 2); hold on
            title('direct stack-short window');
        end
        box on
        
        %% cross-correlation
        %%% directly cross-correlate the template to get the offset between stations
        lohf = 0.5;
        hihf = 6.5;
        npo = 2;
        npa = 2;
        coefmat = zeros(nsta, 1);
        lagmat = zeros(nsta,1);
        aligntemp = zeros(nsta, templen);
        % template is with extra length
        template1 = stackex(1, :)/ nstack;
        % ftemp is filtered with template with extra length
        ftemp1 = Bandpass(template1, sps, lohf, hihf, npo, npa, 'butter');
        %     ftemp1 = template1;
        % temp is 4s win ready for cc
        temp1 = ftemp1(mid+extra-5*sps+1: mid+extra+5*sps);
        
        %%% figure 4
        figure
        for ista = 1: nsta
            
            template2 = stackex(ista, :)/ nstack;
            
            ftemp2 = Bandpass(template2, sps, lohf, hihf, npo, npa, 'butter');
            %         ftemp2 = template2;
            temp2 = ftemp2(mid+extra-5*sps+1: mid+extra+5*sps);
            
            plot(temp2 + 0.1* ista, 'linewidth', 2); hold on
            
            [coef, lag] = xcorr(temp2, temp1, 'coeff'); % 5*sps,
            [maxcoef, idx] = max(coef);
            lagsamp = lag(idx);
            coefmat(ista) = maxcoef;
            lagmat(ista) = lagsamp;
            
            % if lagsamp > 0, shift this station to left relative to the
            % reference station to align them, i.e. the arrival of this sta
            % is late than the reference, so: arr2 == arr1 +lagsamp
            % and vise versa.
            
            aligntemp(ista,:) = template2(1+lagsamp+extra: lagsamp+extra+templen);
        end
        
        figure('Position',[0.1*wid 0.4*hite 0.3*wid 0.3*hite]);
        for ista = 1: nsta
            plot(aligntemp(ista, :) + 0.1* ista, 'linewidth', 2); hold on
            title('aligned direct stack');
        end
        box on
        
        figure('Position',[0.1*wid 0.5*hite 0.3*wid 0.3*hite]);
        for ista = 1: nsta
            plot(aligntemp(ista, mid-5*sps: mid+5*sps) + 0.1* ista, 'linewidth', 2); hold on
            title('aligned direct stack-short window');
        end
        box on
        
        
        %% write into files
        for ista = 1: nsta
            fid = fopen(strcat(temppath,'new_', fam, '_', strtrim(stas(ista, :)), '_', ...
                num2str(sps), 'sps_', num2str(templensec), 'sec_', num2str(nstack), ...
                'DirectStacks_', 'opt_Nofilter_Nonorm.txt'), 'w+');
            fprintf(fid, '%f \n', aligntemp(ista, :)');
            fclose(fid);
        end
        
        %% bandpass and compare the cc results under different freq band
        lohf = 1.25;
        hihf = 6.5;
        npo = 2;
        npa = 2;
        spshf = 40;
        [num, denom] = rat(spshf/40);
        for ista = 1: nsta
            ampmax = max(aligntemp(ista, mid-5*spshf+1: mid+5*spshf));
            ampmin = min(aligntemp(ista, mid-5*spshf+1: mid+5*spshf));
            norm = max(ampmax, -ampmin);
            dummy=aligntemp(ista, :)/norm;
            temphf(ista,:) = Bandpass(dummy, 40, lohf, hihf, npo, npa, 'butter');
            temphf(ista,:)=resample(temphf(ista,:), num, denom); 
        end
        
        figure
        for ista = 1: nsta
            plot(temphf(ista, mid-5*spshf: mid+5*spshf) + 1* ista, 'linewidth', 2); hold on
            text(50,1*ista,stas(ista,:));
        end
        box on
        title('Direct stack-short window-bandpass 1.25-6.5 hz');

        figure
        for ista = 1: nsta
            [coefhf, laghf] = xcorr(temphf(ista,:), temphf(1,:), 'coeff'); % 5*sps,
            [maxcoef, idx] = max(coefhf);
            lagsamp = laghf(idx);
            coefmathf(ista) = maxcoef;
            lagmathf(ista) = lagsamp;
            
            % if lagsamp > 0, shift this station to left relative to the
            % reference station to align them, i.e. the arrival of this sta
            % is late than the reference, so: arr2 == arr1 +lagsamp
            % and vise versa.
            plot(temphf(ista, mid+lagsamp-5*spshf: mid+lagsamp+5*spshf) + 1* ista, 'linewidth', 2); hold on
            text(50,1*ista,stas(ista,:));
        end
        title('Aligned direct stack-short window-bandpass 1.25-6.5 hz');
        lagmathf=lagmathf';
        
    
        lolf = 0.5;
        hilf = 1.25;
        npo = 2;
        npa = 2;
        spslf = 40;
        [num, denom] = rat(spslf/40);
        for ista = 1: nsta
            ampmax = max(aligntemp(ista, mid-5*spshf+1: mid+5*spshf));
            ampmin = min(aligntemp(ista, mid-5*spshf+1: mid+5*spshf));
            norm = max(ampmax, -ampmin);
            dummy=aligntemp(ista, :)/norm;
            dummy1 = Bandpass(dummy, spshf, lolf, hilf, npo, npa, 'butter');
            templf(ista,:)=resample(dummy1, num, denom);

        end
        
%         midlf=608;
        midlf=mid;
        figure
        for ista = 1: nsta
            plot(templf(ista, midlf-5*spslf: midlf+5*spslf) + 1* ista, 'linewidth', 2); hold on
            text(50,1*ista,stas(ista,:));
        end
        box on
        title('Direct stack-short window-bandpass 0.5-1.25 hz');

        figure
        for ista = 1: nsta
            [coeflf, laglf] = xcorr(templf(ista,:), templf(1,:), 'coeff'); % 5*sps,
            [maxcoef, idx] = max(coeflf);
            lagsamp = laglf(idx);
            coefmatlf(ista) = maxcoef;
            lagmatlf(ista) = lagsamp;
            
            % if lagsamp > 0, shift this station to left relative to the
            % reference station to align them, i.e. the arrival of this sta
            % is late than the reference, so: arr2 == arr1 +lagsamp
            % and vise versa.
            plot(templf(ista, midlf+lagsamp-5*spslf: midlf+lagsamp+5*spslf) + 1* ista, 'linewidth', 2); hold on
            text(50,1*ista,stas(ista,:));
        end
        
        title('Aligned direct stack-short window-bandpass 0.5-1.25 hz');
    
    
    
    
    
    
    
    end
    
    
% end

















