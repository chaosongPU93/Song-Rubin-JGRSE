% function mig_linear_fit_LZB_addfam_autortm_v3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script is used to use the poropagation direction estimate fit the hf
% the migration linearly to get the propagation slope and error distribution.
% Use the new automatic time ranges of migrations that come from 'identify_RTMs_v3_requestfam.m'
%
%    2020/09/10, i am adding a fam 001 back again for last try, 13 fam in total now
%                this is for adding new fams 006, 001
%
% Chao Song, chaosong@princeton.edu
% First created date:   2021/01/18
% Last modified date:   2021/01/18
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Initialization
format short e   % Set the format to 5-digit floating point
clear
close all
% clc

% set(0,'DefaultFigureVisible','on');
set(0,'DefaultFigureVisible','off');   % switch to show the plots or not

[scrsz, res] = pixelperinch(1);

workpath = getenv('MHOME');
rstpath = strcat(workpath, '/Seisbasics/hypoinverse/lzbrst');

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
%             '017';
            '006';
            '001';
%             '158';      % 158, 20200916,testing purpose
            ];     % 2020/09/07, i am adding a new family 006 to the pool, final version
        
nfam = size(nfampool,1);
disp(nfam); 
        
famex = '099'

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

  
% load files
winlenhf = 4;

winlenlf = 16;
lofflf = 4;
ccminlf = 0.35;

SUFFIXhf = strcat('up.hf.time.',num2str(winlenhf),'_',num2str(winlenlf),'.',num2str(lofflf),...
                  '.',num2str(ccminlf));
% fname = strcat(rstpath, '/evtloc.requestfam.dcutnodou.',SUFFIXhf);
fname = strcat(rstpath, '/evtloc.ex',famex,'.dcutnodou.',SUFFIXhf);
hftime = load(fname);
% 25 cols, format is:
%   E(043) N(043) E(own) N(own) dist(own) lon lat dep off12 off13 off12sec off13sec date
%   main_arrival_time  cent_of_win  avecc ccadd1 ccadd2 ccadd3 ccadd4 ampmax ampmin eng normeng famnum
%


SUFFIXlf = strcat('lf.time.',num2str(winlenhf),'_',num2str(winlenlf),'.',num2str(lofflf),...
                  '.',num2str(ccminlf));
% fname = strcat(rstpath, '/evtloc.requestfam.dcutnodou.',SUFFIXlf);
fname = strcat(rstpath, '/evtloc.ex',famex,'.dcutnodou.',SUFFIXlf);
lftime = load(fname);

% % this is inverted from (0,0) of all fams, same order, location of control points
% loccont = [-123.492667 48.451500 38.1400; 
%            -123.772167 48.493000 35.5900; 
%            -123.863167 48.528167 35.2100;
%            -123.603333 48.440167 36.7100;
%            -123.800167 48.408833 34.5200;
%            -123.893333 48.536500 35.0700;
%            -123.864500 48.498667 34.8800;
%            -123.753333 48.525667 36.2000;
%            -123.703667 48.502667 36.4100;
%            -123.814333 48.538667 35.7900;
% %            -123.838500 48.544833 35.6600;
%            -123.908000 48.494167 34.5100;       % 006
%            -123.879667 48.446167 34.2600;       % 001
% %            -123.967000 48.458667 33.7800;       % 158, 20200916,testing purpose
%            ];
%        
% 
% relacont = loccont;
% [dx, dy] = absloc2relaloc(loccont(:,1),loccont(:,2), -123.772167, 48.493000);
% relacont(:,1) = dx;
% relacont(:,2) = dy;

%%% load the location resolution points, +-1 & 2 samples
reslzb1 = load(fullfile(rstpath,'evtloc.13fam_locres_hf1spl'));
reslzb2 = load(fullfile(rstpath,'evtloc.13fam_locres_hf2spl'));

% sort according to day, sec
hftime = sortrows(hftime, [13, 15]);
lftime = sortrows(lftime, [13, 15]);

% this is new time ranges that come from 'identify_RTMs_v3.m'
trange = [
    2004195   4.1773e+04   4.3362e+04   % combine, y
    2004197   4.0974e+04   4.2205e+04   % speed direct 90, y, larger pear
    2004197   5.4453e+04   5.5867e+04   % speed direct 205, y, but use rmse 220
    2004197   5.6270e+04   5.7227e+04   % speed direct 175, y, but use rmse 185
%     2004197   7.9949e+04   8.1453e+04   % speed direct 260, y, but use rmse 195
    2004197   8.3402e+04   8.4888e+04   % divided into 3, 1 discarded, y
    2004197   8.5320e+04   8.6659e+04   % divided into 3, 1 discarded, y
    2004198   5.9530e+03   9.7940e+03   % speed direct 270, y, but use rmse 235
    2004198   1.9151e+04   2.1616e+04   % acceptted, y
    2004198   4.2632e+04   4.3795e+04   % acceptted, y
    2004198   5.3506e+04   5.8050e+04   % acceptted, y
    2004198   6.2966e+04   6.3885e+04   % acceptted, y
    2004198   7.3540e+04   7.6029e+04   % acceptted, y 
    2004198   8.4222e+04   8.6400e+04   % combine, y
    2004199   4.1670e+03   6.2660e+03   % speed direct 280, y, larger pear
    2004199   4.2845e+04   4.3359e+04   % speed direct 270, y, but use rmse 250
    2004199   4.6340e+04   4.9126e+04   % combine, care speed direct, y, larger pear
    2004199   4.9744e+04   5.0614e+04   % speed direct 180, y, but use rmse 195
    2004199   8.0861e+04   8.2008e+04   % divided into 2, y, care speed direct, y, but use rmse 90
    2004199   8.2008e+04   8.3515e+04   % divided into 2, y, care speed direct, RECHECK PEAR, slightly larger pear    
    2004200   1.2600e+04   1.5799e+04   % modified time, y
    2004200   1.5914e+04   1.7125e+04   % acceptted, y
    2004200   1.9104e+04   1.9900e+04   % check OLD param. to determine, y, use OLD end time
    2004200   4.8109e+04   4.8666e+04   % acceptted, y
    2004201   4.3700e+02   1.6290e+03   % speed direct 225, y, but use rmse 260
    2004203   1.6586e+04   2.0776e+04   % combine them, y
%     2005255   3.4400e+04   3.5612e+04   % speed direct 80, y, larger pear
    2005255   6.0381e+04   6.1557e+04   % speed direct 150, y, but use rmse 110
    2005255   6.7268e+04   6.8535e+04   % speed direct 80, y, larger pear
    2005255   8.4296e+04   8.5634e+04   % combine them, y
%     2005256   3.6200e+03   4.9570e+03   % speed direct 265, y, larger pear
    2005256   8.0330e+03   9.8280e+03   % change end time, y
    2005256   1.3412e+04   1.4346e+04   % speed direct 65, y, but use rmse 80   
%     2005256   1.9069e+04   2.0203e+04   % speed direct 195, y, larger pear
    2005256   2.1492e+04   2.2294e+04   % change start time, y
    2005256   2.8080e+04   2.9172e+04   % change start time, y
    2005256   3.4776e+04   3.5238e+04   % change start time, y
    2005256   3.6900e+04   3.8520e+04   % change times, y
%     2005256   4.6440e+04   4.8421e+04   % change start time, y, check speed direc, y, a bit awful
    2005256   5.2020e+04   5.2596e+04   % change times, y
    2005256   6.2640e+04   6.4865e+04   % change start time, y, check speed direc, y, but use rmse 130
    2005256   7.2720e+04   7.3914e+04   % divided into 2, y, but check speed direc, y
    2005256   7.6381e+04   7.7940e+04   % change end time, y
    2005257   6.1200e+03   7.7460e+03   % change start time, y, but check speed direc, y
%     2005257   1.2240e+04   1.3248e+04   % CHANGED time, st-3.68, y, use speed direc, larger pear, SUSPICOUS
    2005257   1.6920e+04   1.8655e+04   % change start time, y
    2005257   2.1070e+04   2.5365e+04   % accepted, y 
    2005257   3.9000e+04   4.5418e+04   % change start time, y
    2005257   6.1449e+04   6.4012e+04   % accepted, y
    2005257   7.3440e+04   7.8840e+04   % change end time, y, but check speed direc, y, but use rmse 235
    2005257   8.2159e+04   8.3380e+04   % speed direct 350, y, larger pear
    2005258   3.5000e+04   3.6540e+04   % divided into 3, y, but check speed direc, RECHECK PEAR, slightly larger pear, SUSPICOUS
    2005258   3.6540e+04   3.8160e+04   % divided into 3, y
    2005258   3.8160e+04   4.0021e+04   % divided into 3, y, but check speed direc y, slightly larger pear 
%     2005259   1.5400e+02   1.6340e+03   % speed direct 75, y, larger pear
%     2005259   1.9560e+03   4.0460e+03   % speed direct 185, y, but use rmse 140
    2005259   5.2320e+03   7.7140e+03   % check OLD param. to determine, y, but check speed direc, larger pear, use old direc 255
    2005259   3.7000e+04   4.1293e+04   % speed direct 85, RECHECK PEAR, y, but use rmse 115
%     2005260   2.6510e+03   3.4380e+03   % speed direct 80, y, but use rmse 115
    2005260   3.6040e+03   4.5000e+03   % divided into 2, y
%     2005260   4.6080e+03   6.6780e+03   % divided into 2, y, but check speed direc, y, larger pear
    2005260   5.6300e+04   5.8200e+04   % use old times, y
%     2005260   9.4090e+03   1.0357e+04   % speed direct 55, y, but use rmse 115
    ]; 

%%%%% for cases of different fam excluded
if strcmp(famex, '017') 
    %%% 1. without fam 017
    angbest = [250,90,220,185,230,270,235,200,230,80,70,235,245,310,250,120,195,90,140,130,70,255,250,260,120,110,80,235,215,80,115,105,260,95,260,130,230,250,195,115,245,170,245,240,10,95,85,185,255,115,250,240];
    indspeed = [2,27,34,37,47];
    indlargediff = [2,7,11,14,24,28,29,32,34,36,39,40,43,49];
elseif strcmp(famex, '043')
    %%% 2. without fam 043
    angbest = [250,115,200,100,260,270,240,200,225,60,70,265,240,310,255,120,170,90,140,130,130,255,265,260,120,100,85,300,260,125,190,105,260,90,260,125,230,250,200,115,250,165,245,235,15,95,155,175,255,115,250,250];
    indspeed = [4,5,12,21,27,28,29,34,37,52];
    indlargediff = [4,5,11,13,14,24,27,28,29,32,34,39,40,43];
elseif strcmp(famex, '068')
    %%% 3. without fam 068
    angbest = [250,90,220,185,225,270,230,200,225,80,80,205,245,305,250,120,195,90,140,130,90,260,250,255,120,110,80,240,220,80,115,105,260,90,260,130,230,250,195,115,245,165,240,235,0,95,85,185,255,115,250,240];
    indspeed = [2,21,27,34,37,47,52];
    indlargediff = [2,9,32,36,40,44,52];
elseif strcmp(famex, '099')
    %%% 4. without fam 099
    angbest = [250,95,215,185,225,270,255,200,230,80,65,235,245,310,250,120,195,90,140,130,70,255,245,260,125,110,80,240,220,80,115,105,265,90,245,135,225,250,195,115,250,165,245,235,10,60,85,185,255,120,250,240];
    indspeed = [2,27,34,37,46,47];
    indlargediff = [2,11,14,24,32,34,36,40,43];
end


%% linear fitting and median distance 
resprophf = nan(length(trange),1200);
resproplf = nan(length(trange),500);

medallhf = nan(length(trange),4);
medalllf = nan(length(trange),4);
ranvechf = nan(length(trange),2);
ranveclf = nan(length(trange),2);
ranvechf95 = nan(length(trange),2);
ranveclf95 = nan(length(trange),2);
ranvechf98 = nan(length(trange),2);
ranveclf98 = nan(length(trange),2);
ranvechf99 = nan(length(trange),2);
ranveclf99 = nan(length(trange),2);

% array of vertical distance from LF from fitted HF line 
vertdistraw = nan(300 , length(trange));    % raw vertical distance
vertdistwt = nan(300 , length(trange));     % weighted vertical distance
weight = nan(300 , length(trange));     % weights from regression
indlfindpen = nan(300 , length(trange));    % store the index of indepedent LF detections
indltdiv = nan(300 , length(trange));    % store the index that is lower than the division
indgtdiv = nan(300 , length(trange));    % store the index that is greater than the division

% percentage of detections from specific fam relative to all in each RTM
famperchf = zeros(length(trange), nfam);
famperclf = zeros(length(trange), nfam);

xran = [-20 25];
yran = [-20 20];

angle = 0:5:355;

slopehf = zeros(size(trange,1),length(angle));
rmsehf = zeros(size(trange,1),length(angle));
angrmsehf = zeros(size(trange,1),1);
angslopehf = zeros(size(trange,1),1);

% create fit object with free constraints
fttpfree = fittype( @(a,b,x) a*x+b);

% a box around fam 002 region
reg002 = [ -5 -15;
            25 0
            35 -20;
            5 -35;
           -5 -15];
       
% index of migrations whose majority is NOT at 002 region but have some detections there  
indmig002 = [2,4,5,6,7,8,9,12,13,14,17,19,20,24,25,28,29,31,32,35,36,37,39,41,42,43,44,46,47,48,49,...
             50,52];
         
for i = 1: length(trange)
% for i = 1: 211
%     i=46;
    disp(trange(i,:));
    indhf = find(hftime(:,13)==trange(i,1) & hftime(:,15)>=trange(i,2) & ...
                 hftime(:,15)<=trange(i,3));
    mighf = hftime(indhf,:);
    
    indlf = find(lftime(:,13)==trange(i,1) & lftime(:,15)>=trange(i,2) & ...
                 lftime(:,15)<=trange(i,3));
    miglf = lftime(indlf,:); 
    
    %%% find best angle, now is mainly to get the variation of se and slope with the trial angle
    for iang = 1: length(angle)
        %%% propagation trial of hf
        mighfdum = mighf;
        for j = 1: size(mighf,1)
            x0 = mighfdum(j,1);
            y0 = mighfdum(j,2);
            [newx,newy] = coordinate_rot(x0,y0,-(angle(iang)-90),0,0);
            mighfdum(j,1) = newx;
            mighfdum(j,2) = newy;
        end        
        % linear robust least square
        [fitobj,gof,~] = fit(mighfdum(:,15)/3600, mighfdum(:,1),fttpfree,'Robust','Bisquare',...
                                  'StartPoint',[1 1]);
        % output fit parameters
        coef = coeffvalues(fitobj);
        slopehf(i,iang) = coef(1);
        rmsehf(i,iang) = gof.rmse;
            
    end

    %%% best angle estimate from hf
    ind = find(slopehf(i,:)>0);
    ind3 = find(rmsehf(i,ind)==min(rmsehf(i,ind)));     % rmse, Root Mean Squared Error, the smaller, the better
    if length(ind3) > 1
        disp(strcat({'multiple angles have same rmse: '},num2str(angle(ind3))));
    end
    angrmsehf(i) = angle(ind(ind3(1)));

    ind6 = find(slopehf(i,:)==max(slopehf(i,:))); % one with the largest slope, i.e., migrating speed
    if length(ind6) > 1
        disp(strcat({'multiple angles have same slope:'},num2str(angle(ind6))));
    end
    angslopehf(i) = angle(ind6(1));
    
    
    %%% Actually used is the pre-determined best prop direc to do the fitting
    mighfdum = mighf;
    miglfdum = miglf;
    
%     % for the migrations whose majority is occuring not at 002 region, leave out 
%     % the detections at 002 region
%     if ~isempty(intersect(indmig002,i))
%         % only use the detections ouside the bound
%         [is,ion] = inpolygon(mighfdum(:,1),mighfdum(:,2),reg002(:,1),reg002(:,2));
%         isinreg = is | ion;
%         mighfdum = mighfdum(isinreg ~= 1, :);
%         
%         [is,ion] = inpolygon(miglfdum(:,1),miglfdum(:,2),reg002(:,1),reg002(:,2));
%         isinreg = is | ion;
%         miglfdum = miglfdum(isinreg ~= 1, :);
%     end
    
    for j = 1: size(mighfdum,1)
        x0 = mighfdum(j,1);
        y0 = mighfdum(j,2);
        [newx,newy] = coordinate_rot(x0,y0,-(angbest(i)-90),0,0);
        mighfdum(j,1) = newx;
        mighfdum(j,2) = newy;
    end

    for j = 1: size(miglfdum,1)
        x0 = miglfdum(j,1);
        y0 = miglfdum(j,2);
        [newx,newy] = coordinate_rot(x0,y0,-(angbest(i)-90),0,0);
        miglfdum(j,1) = newx;
        miglfdum(j,2) = newy;
    end
    
    % detections involved into fitting 
    mintime = max(min(mighfdum(:,15)), min(miglfdum(:,15)));
    maxtime = min(max(mighfdum(:,15)), max(miglfdum(:,15)));
    mighfdum2 = mighfdum(mighfdum(:,15)>=mintime & mighfdum(:,15)<=maxtime, :);
    miglfdum2 = miglfdum(miglfdum(:,15)>=mintime & miglfdum(:,15)<=maxtime, :);

    % get the percentage of detections from specific fam
    for j = 1: nfam
        famperchf(i,j) = sum(mighfdum2(:,end) == str2double(nfampool(j,:))) / size(mighfdum2,1);
        famperclf(i,j) = sum(miglfdum2(:,end) == str2double(nfampool(j,:))) / size(miglfdum2,1);
    end
    
    % assign the location resolution from one fam
%     ifam = find(famperchf(i,:) == max(famperchf(i,:))); % choose fam contributes most detections
    ifam = 2;   % choose fam 043 as the representative
    tmp1 = reslzb1((ifam-1)*4+1: ifam*4, :);
    [dx, dy] = absloc2relaloc(tmp1(:,1),tmp1(:,2),-123.772167, 48.493000);
    vecreslzb1 = [dx dy];
    tmp2 = reslzb2((ifam-1)*4+1: ifam*4, :);
    [dx, dy] = absloc2relaloc(tmp2(:,1),tmp2(:,2),-123.772167, 48.493000);
    vecreslzb2 = [dx dy];

    %%% define and position the figure frame and axes of each plot
    f.fig=figure;
    widin = 8;  % maximum width allowed is 8.5 inches
    htin = 9;   % maximum height allowed is 11 inches
    set(f.fig,'Position',[scrsz(1)+1*scrsz(3)/10 scrsz(2)+scrsz(4)/20 widin*res htin*res]);
    nrow = 2;
    ncol = 2;    
    for isub = 1:nrow*ncol
        f.ax(isub) = subplot(nrow,ncol,isub);
    end

    %%% reposition
    set(f.ax(1), 'position', [ 0.08, 0.64, 0.36, 0.32]);
    set(f.ax(2), 'position', [ 0.52, 0.64, 0.36, 0.32]);
    set(f.ax(3), 'position', [ 0.08, 0.35, 0.36, 0.22]);
    set(f.ax(4), 'position', [ 0.52, 0.35, 0.36, 0.22]);
%     set(f.ax(5), 'position', [ 0.1, 0.078, 0.32, 0.22]);
%     set(f.ax(6), 'position', [ 0.5, 0.078, 0.32, 0.22]);
    
    % subplot 1 of figure i
    hold(f.ax(1),'on');
    plot(f.ax(1),[-100 100],[0 0],'k--');
    plot(f.ax(1),[0 0],[-100 100],'k--');
    f.ax(1).FontSize = 9;
    mighf = sortrows(mighf,-15);
    scatter(f.ax(1),mighf(:,1),mighf(:,2), 20, mighf(:,15)/3600, 'filled','o');  %, 'MarkerEdgeColor', 'w')
    colormap(f.ax(1),'jet');
    c=colorbar(f.ax(1),'SouthOutside');
    pos = f.ax(1).Position;
    c.Position = [pos(1), pos(2)-0.018, pos(3), 0.02];
%     c.TickLabels=[];
    juldate = num2str(trange(i,1));
    yr = str2double(juldate(1:4));
    date = str2double(juldate(5:end));
    a = jul2dat(yr,date);
    mo = a(1);
    if mo == 9 
        mo = {' Sep. '};
    elseif mo == 7
        mo = {' Jul. '};
    else
        mo = {' Mar. '};
    end
    day = num2str(a(2));
    yr = num2str(a(3));
    c.Label.String = strcat({'Time (hr) on '}, day, mo, yr);
%     c.Label.String = strcat(num2str(trange(i,1)),' of HF',' (hr)');
    c.Label.FontSize = 11;
    caxis(f.ax(1),[trange(i,2)/3600 trange(i,3)/3600])
    text(f.ax(1),0.85,0.1,'HF','FontSize',12,'unit','normalized');
    text(f.ax(1),0.04,0.93,'a','FontSize',11,'unit','normalized','EdgeColor','k','Margin',2);
    text(f.ax(1),0.04,0.1,'TWKB','FontSize',11,'unit','normalized','EdgeColor','k','Margin',2);
    f.ax(1).Box = 'on';
    grid(f.ax(1), 'on');
    axis(f.ax(1), 'equal');
    f.ax(1).GridLineStyle = '--';
    f.ax(1).XAxisLocation = 'top';
    medxhf = median(mighf(:,1));
    medyhf = median(mighf(:,2));
    [rotx, roty] = complex_rot(0,5,-angbest(i));
    xvect = [medxhf-rotx medxhf+rotx];
    yvect = [medyhf-roty medyhf+roty];
    drawArrow(f.ax(1),xvect,yvect,xran,yran,'linewidth',1.5);
    [rotx, roty] = complex_rot(-5,0,-angbest(i));
    xvect = [medxhf-rotx medxhf+rotx];
    yvect = [medyhf-roty medyhf+roty];
    drawArrow(f.ax(1),xvect,yvect,xran,yran,'linewidth',1.5,'linestyle','--','color',[0.4 0.4 0.4]);
    x = [16+ vecreslzb1(1,1), 16+ vecreslzb1(2,1)];
    y = [6+ vecreslzb1(1,2), 6+ vecreslzb1(2,2)];  % loc resolution indicated by +-1 sample cross
    plot(f.ax(1),x,y,'color','r','linewidth',1);
    x = [16+ vecreslzb1(3,1), 16+ vecreslzb1(4,1)];
    y = [6+ vecreslzb1(3,2), 6+ vecreslzb1(4,2)];  % loc resolution indicated by +-1 sample cross
    plot(f.ax(1),x,y,'color','r','linewidth',1);
    xticks(f.ax(1),xran(1):5:xran(2));
    yticks(f.ax(1),yran(1):5:yran(2));    
    xlabel(f.ax(1),'E (km)','fontsize',11);
    ylabel(f.ax(1),'N (km)','fontsize',11);
    text(f.ax(1),0.5,0.1,num2str(i),'FontSize',12,'unit','normalized');
    text(f.ax(1),0.84,0.95,strcat(num2str(size(mighf,1)),{' detections'}),'FontSize',8,...
         'unit','normalized','horizontalalignment','center');
    text(f.ax(1),0.84,0.90,strcat({'in '},num2str(trange(i,3)-trange(i,2)),{' s'}),'FontSize',8,...
         'unit','normalized','horizontalalignment','center');
    rate = sprintf('%.3f',size(mighf,1)/(trange(i,3)-trange(i,2)));
    text(f.ax(1),0.84,0.85,strcat({'rate: '},rate),'FontSize',8,...
         'unit','normalized','horizontalalignment','center');
    text(f.ax(1),0.84,0.80,strcat(num2str(angbest(i)),{'{\circ}'}),'FontSize',10,...
         'unit','normalized','horizontalalignment','center'); 
    hold(f.ax(1),'off');

    % subplot 2 of figure i
    hold(f.ax(2),'on');
    plot(f.ax(2),[-100 100],[0 0],'k--');
    plot(f.ax(2),[0 0],[-100 100],'k--');
    f.ax(2).FontSize = 9;
    miglf = sortrows(miglf,-15);
    scatter(f.ax(2),miglf(:,1),miglf(:,2), 20, miglf(:,15)/3600, 'filled','o');  %, 'MarkerEdgeColor', 'w')
    colormap(f.ax(2),'jet');
    c=colorbar(f.ax(2),'SouthOutside');
    pos = f.ax(2).Position;
    c.Position = [pos(1), pos(2)-0.018, pos(3), 0.02];
%     c.TickLabels=[];
    c.Label.String = strcat({'Time (hr) on '}, day, mo, yr);
%     c.Label.String = strcat(num2str(trange(i,1)),' of LF',' (hr)');
    c.Label.FontSize = 11;
    caxis(f.ax(2),[trange(i,2)/3600 trange(i,3)/3600])
    text(f.ax(2),0.85,0.1,'LF','FontSize',12,'unit','normalized');
    text(f.ax(2),0.04,0.93,'b','FontSize',11,'unit','normalized','EdgeColor','k','Margin',2);
    text(f.ax(2),0.04,0.1,'TWKB','FontSize',11,'unit','normalized','EdgeColor','k','Margin',2);
    f.ax(2).Box = 'on';
    grid(f.ax(2), 'on');
    axis(f.ax(2), 'equal');
    f.ax(2).GridLineStyle = '--';
    f.ax(2).XAxisLocation = 'top';
    medxlf = median(miglf(:,1));
    medylf = median(miglf(:,2));
    [rotx, roty] = complex_rot(0,5,-angbest(i));
    xvect = [medxlf-rotx medxlf+rotx];
    yvect = [medylf-roty medylf+roty];
    drawArrow(f.ax(2),xvect,yvect,xran,yran,'linewidth',1.5);
    [rotx, roty] = complex_rot(-5,0,-angbest(i));
    xvect = [medxlf-rotx medxlf+rotx];
    yvect = [medylf-roty medylf+roty];
    drawArrow(f.ax(2),xvect,yvect,xran,yran,'linewidth',1.5,'linestyle','--','color',[0.4 0.4 0.4]);
    x = [16+ vecreslzb2(1,1), 16+ vecreslzb2(2,1)];
    y = [6+ vecreslzb2(1,2), 6+ vecreslzb2(2,2)];  % loc resolution indicated by +-1 sample cross
    plot(f.ax(2),x,y,'color','r','linewidth',1);
    x = [16+ vecreslzb2(3,1), 16+ vecreslzb2(4,1)];
    y = [6+ vecreslzb2(3,2), 6+ vecreslzb2(4,2)];  % loc resolution indicated by +-1 sample cross
    plot(f.ax(2),x,y,'color','r','linewidth',1);
    xticks(f.ax(2),xran(1):5:xran(2));
    yticks(f.ax(2),yran(1):5:yran(2));    
    xlabel(f.ax(2),'E (km)','fontsize',11);
    ylabel(f.ax(2),'N (km)','fontsize',11);
    text(f.ax(2),0.5,0.1,num2str(i),'FontSize',12,'unit','normalized');
    text(f.ax(2),0.84,0.95,strcat(num2str(size(miglf,1)),{' detections'}),'FontSize',8,...
         'unit','normalized','horizontalalignment','center');
    text(f.ax(2),0.84,0.9,strcat({'in '},num2str(trange(i,3)-trange(i,2)),{' s'}),'FontSize',8,...
         'unit','normalized','horizontalalignment','center');
    rate = sprintf('%.3f',size(miglf,1)/(trange(i,3)-trange(i,2)));
    text(f.ax(2),0.84,0.85,strcat({'rate: '},rate),'FontSize',8,...
         'unit','normalized','horizontalalignment','center');
    text(f.ax(2),0.84,0.80,strcat(num2str(angbest(i)),{'{\circ}'}),'FontSize',10,...
         'unit','normalized','horizontalalignment','center');
    hold(f.ax(2),'off');
    
    % subplot 3 of figure i
    hold(f.ax(3),'on');
    f.ax(3).FontSize = 9;
    scatter(f.ax(3),mighfdum(:,15)/3600,mighfdum(:,1),20,[0.6 0.6 0.6],'filled','o',...
            'MarkerEdgeColor','k');    
    scatter(f.ax(3),miglfdum(:,15)/3600,miglfdum(:,1),20,[0.6 1 1],'filled','o',...
            'MarkerEdgeColor','k');
    % create fit object
        
    [fitobjhfprop,gofhfprop,outphfprop] = fit(mighfdum2(:,15)/3600, mighfdum2(:,1),fttpfree,'Robust',...
                                              'Bisquare','StartPoint',[1 1]);
    % output fit parameters
    coefprop = coeffvalues(fitobjhfprop);
    slopeprophf(i) = coefprop(1);
    intcptprophf(i) = coefprop(2);
    fitprophf = feval(fitobjhfprop,mighfdum2(:,15)/3600);
    plot(f.ax(3),mighfdum2(:,15)/3600,fitprophf,'-','linewidth',2,'color',[0.6 0.6 0.6]);
        
%     intcptproplf = ones(size(miglfdum(:,end)))\(miglfdum(:,1)-slopeprophf*miglfdum(:,end)/3600); % direct LS fit
%     fitproplf = slopeprophf*miglfdum(:,end)/3600+intcptproplf;
    a=slopeprophf(i);
    fttpfix = fittype( @(b,x) a*x+b);
    [fitobjlfprop,goflfprop,outplfprop] = fit(miglfdum2(:,15)/3600, miglfdum2(:,1),fttpfix,'Robust',...
                                              'Bisquare','StartPoint',1);
%     [fitobjlfprop,goflfprop,outplfprop] = fit(miglfdum(:,10)/3600, miglfdum(:,1),fttpfix,'Robust',...
%                                               'Bisquare','Weights',miglfdum(:,11));
    fitproplf = feval(fitobjlfprop,miglfdum2(:,15)/3600);
    intcptproplf(i) = coeffvalues(fitobjlfprop);
    offset(i) = intcptproplf(i)-intcptprophf(i);
    plot(f.ax(3),miglfdum2(:,15)/3600,fitproplf,'-','linewidth',2,'color',[0.6 1 1]);
    f.ax(3).Box = 'on';
    grid(f.ax(3), 'on');
    f.ax(3).GridLineStyle = '--';
    xran1 = [trange(i,2)/3600 trange(i,3)/3600];
%     yran1 = [round(min([mighfdum(:,1);miglfdum(:,1)]))-1 ...
%              round(max([mighfdum(:,1);miglfdum(:,1)]))+1];
    aa = round(prctile([mighfdum(:,1);miglfdum(:,1)], 98));
    bb = round(prctile([mighfdum(:,1);miglfdum(:,1)], 2));
    aap = aa + ceil((aa-bb)/6);
    bbp = bb - ceil((aa-bb)/3);
    yran1 = [bbp aap];
    xlim(f.ax(3),xran1);
    ylim(f.ax(3),yran1);
    text(f.ax(3),0.04,0.91,'c','FontSize',11,'unit','normalized','EdgeColor','k','Margin',2);
    xlabel(f.ax(3),'Time (hr)','fontsize',11);
    ylabel(f.ax(3),'Dist. along prop. (km)','fontsize',11); 
    
    % compute the HF weights in robust linear regression, see NOTES above
    rhf = outphfprop.residuals;   % usual residuals
    x = [mighfdum2(:,15)/3600];
    hatmat = x*inv(x'*x)*x';
    h = zeros(size(hatmat,1),1);    % leverage of least square
    for jj = 1 : size(hatmat,1)
        h(jj) = hatmat(jj,jj);
    end    
    radj = rhf./sqrt(1-h);      % adjusted residuals
    K = 4.685;
    s = mad(rhf,1)/0.6745;
    u = radj/(K*s);
    wthf = zeros(length(u),1);    % rubust weight of next iteration
    for jj = 1 : length(u)
        if abs(u(jj)) < 1
            wthf(jj) = (1-(u(jj))^2)^2;
        else
            wthf(jj) = 0;
        end
    end
    
    % get the standard error of the estimated parameters, may indicate the compare the quality
    % of fitting cross different RTMs, NOTE that rmse, mse are definitely not a good measure
    % the obtained CI is comfirmed to be correct by comparing it with the 'fitobjhfprop'
    slopesehf = gofhfprop.rmse./sqrt(sum((x-mean(x)).^2));
    est = slopeprophf(i);
    slopeprophfCI(i,:) = confidence_interval_general(est,slopesehf,length(x)-2,95);
    interceptsehf = slopesehf.*sqrt(sum(x.^2)./length(x));
    est = intcptprophf(i);
    intcptprophfCI(i,:) = confidence_interval_general(est,interceptsehf,length(x)-2,95);
    sehf(i, :) = [gofhfprop.rmse slopesehf interceptsehf];     
    
    x = miglfdum2(:,15)/3600;
    slopeself = goflfprop.rmse./sqrt(sum((x-mean(x)).^2));
    interceptself = slopeself.*sqrt(sum(x.^2)./length(x));
    self(i, :) = [goflfprop.rmse slopeself interceptself];
    
    x = mighfdum2(:,15)/3600;
    y = mighfdum2(:,1);
    x_bar = wt_mean(x,wthf);
    y_bar = wt_mean(y,wthf);
    x_var = sum(wthf.*(x-x_bar).^2) / sum(wthf);
    y_var = sum(wthf.*(y-y_bar).^2) / sum(wthf);
    xy_cov = sum(wthf.*(x-x_bar).*(y-y_bar)) / sum(wthf);
    pearwthf(i) = xy_cov / sqrt(x_var*y_var);
    
    % compute the LF weights in robust linear regression, see NOTES above
    rlf = outplfprop.residuals;   % usual residuals
    x = [miglfdum2(:,15)/3600];
    hatmat = x*inv(x'*x)*x';
    h = zeros(size(hatmat,1),1);    % leverage of least square
    for jj = 1 : size(hatmat,1)
        h(jj) = hatmat(jj,jj);
    end    
    radj = rlf./sqrt(1-h);      % adjusted residuals
    K = 4.685;
    s = mad(rlf,1)/0.6745;
    u = radj/(K*s);
    wtlf = zeros(length(u),1);    % rubust weight of next iteration
    for jj = 1 : length(u)
        if abs(u(jj)) < 1
            wtlf(jj) = (1-(u(jj))^2)^2;
        else
            wtlf(jj) = 0;
        end
    end
    
    % get the vertical distance from LF detections to the fitted HF line
    fittemp = feval(fitobjhfprop,miglfdum2(:,15)/3600);
    verttemp = miglfdum2(:,1) - fittemp;     % can be + or -
    vertdistraw(1: length(verttemp), i) = verttemp;     % raw vertical distance
    verttempwt = wtlf.*verttemp;      % weighted vertical distance by multiplying the weight
    vertdistwt(1: length(verttempwt), i) = verttempwt;
    weight(1: length(verttempwt), i) =  wtlf;
       
    % get the information of propagation vector indicating range etc.
    [tmpx,tmpy] = coordinate_rot(medxhf,medyhf,-(angbest(i)-90),0,0);
    medallhf(i,1:4) = [medxhf medyhf tmpx tmpy];
    [tmpx,tmpy] = coordinate_rot(medxlf,medylf,-(angbest(i)-90),0,0);
    medalllf(i,1:4) = [medxlf medylf tmpx tmpy];
    
    ranvechf(i,1:2) = [min(mighfdum(:,1)) max(mighfdum(:,1))];
    ranveclf(i,1:2) = [min(miglfdum(:,1)) max(miglfdum(:,1))];
    
    ranvechf95(i,1:2) = [prctile(mighfdum(:,1),5) prctile(mighfdum(:,1),95)];
    ranveclf95(i,1:2) = [prctile(miglfdum(:,1),5) prctile(miglfdum(:,1),95)];
    
    ranvechf98(i,1:2) = [prctile(mighfdum(:,1),2) prctile(mighfdum(:,1),98)];
    ranveclf98(i,1:2) = [prctile(miglfdum(:,1),2) prctile(miglfdum(:,1),98)];
    
    nlt = sum(mighfdum(:,1)< ranvechf98(i,1));      % number that lower than the current range
    ngt = sum(mighfdum(:,1)> ranvechf98(i,2));      % number that greater than the current range
    
    mighfsort = sort(mighfdum(:,1));    % sort in descend order
    if nlt < 2
        ranvechf98(i,1) = mighfsort(3);     % so at least throw away 2 abnormals
    end
    
    if ngt < 2
        ranvechf98(i,2) = mighfsort(end-2);     % so at least throw away 2 abnormals
    end
    
    ranvechf99(i,1:2) = [prctile(mighfdum(:,1),1) prctile(mighfdum(:,1),99)];
    ranveclf99(i,1:2) = [prctile(miglfdum(:,1),1) prctile(miglfdum(:,1),99)];
                
    %%% find the independent LF detections, i.e. non-overlapping time window
    [miglfindpen,indindpen] = find_independent_detection(miglfdum2, 0);
    inumlf(i) = length(indindpen);
    indlfindpen(1: inumlf(i), i) =  indindpen;
    
    text(f.ax(3),0.45,0.2,sprintf('Slope: %.1f km/h',slopeprophf(i)),'FontSize',8,...
         'unit','normalized','horizontalalignment','left');
    text(f.ax(3),0.97,0.2,sprintf('SE: %.2f',slopesehf),'FontSize',8,...
         'unit','normalized','horizontalalignment','right'); 
    text(f.ax(3),0.45,0.13,sprintf('Pearson: %.3f',pearwthf(i)),'FontSize',8,...
         'unit','normalized','horizontalalignment','left');
    text(f.ax(3),0.97,0.12,sprintf('SE_{LF}: %.2f',slopeself),'FontSize',8,...
         'unit','normalized','horizontalalignment','right');  
    text(f.ax(3),0.45,0.05,sprintf('LF-HF: %.2f km',offset(i)),'FontSize',10,'unit','normalized'); 
    hold(f.ax(3),'off');    
    
    % subplot 4 of figure i
    hold(f.ax(4),'on');
    f.ax(4).FontSize = 9;
    f.ax(4).Box = 'on';
    yyaxis(f.ax(4),'left');
    plot(f.ax(4),angle,rmsehf(i,:),'o-','linew',1,'color','k','markers',4);
    yranl = f.ax(4).YLim;    
    if angbest(i) == angrmsehf(i)
        plot(f.ax(4),[angrmsehf(i) angrmsehf(i)],[0 100],'--','linew',1.5,'color','k');
    else
        plot(f.ax(4),[angrmsehf(i) angrmsehf(i)],[0 100],'--','linew',0.5,'color','k');
    end
    f.ax(4).YColor = 'k';
    f.ax(4).YLim = yranl;
    ylabel(f.ax(4),'RMSE of HF','fontsize',11);
    
    yyaxis(f.ax(4),'right');
    plot(f.ax(4),angle,slopehf(i,:),'^:','linew',1,'color',[0.12 0.56 1],'markers',4);
    yranr = f.ax(4).YLim;
    if angbest(i) == angslopehf(i)
        plot(f.ax(4),[angslopehf(i) angslopehf(i)],f.ax(4).YLim,'--','linew',1.5,...
             'color',[0.12 0.56 1]);    %[0 0 0.5];
    else
        plot(f.ax(4),[angslopehf(i) angslopehf(i)],f.ax(4).YLim,'--','linew',0.5,...
             'color',[0.12 0.56 1]);
    end    
    f.ax(4).YColor = [0.12 0.56 1]; %[0 0 0.5];
    f.ax(4).YLim = yranr;
    ylabel(f.ax(4),'Slope of HF (km/h)','fontsize',11);
    text(f.ax(4),0.04,0.91,'d','FontSize',11,'unit','normalized','EdgeColor','k','Margin',2);
    xlim(f.ax(4),[0 360]);
%     ymax = f.ax(4).YLim(2)+0.1;
%     ylim(f.ax(4),[0 ymax]);
%     ylim(f.ax(4),[0 1]);
    xlabel(f.ax(4),'Trial propagation direction (deg)','fontsize',11);
    hold(f.ax(4),'off');
    
%     % subplot 4 of figure i
%     hold(f.ax(4),'on');
%     f.ax(4).FontSize = 9;
    numhf(i) = length(mighfdum2(:,1));
    numlf(i) = length(miglfdum2(:,1));
    resprophf(i,1:numhf(i)) = outphfprop.residuals;
    resproplf(i,1:numlf(i)) = outplfprop.residuals+(intcptproplf(i)-intcptprophf(i));
%     [xlochf, a, b, pdfvalhf] = weighted_bar_pdf(resprophf(i,1:numhf(i)), wthf, 0.5, 'int');
%     hfHdl = bar(f.ax(4),xlochf, pdfvalhf,1,'stacked');
%     hfHdl(1).FaceColor = [0.6 0.6 0.6];
%     [xloclf, ~, ~, pdfvallf] = weighted_bar_pdf(resproplf(i,1:numlf(i)), wtlf, 0.5, 'int');
%     lfHdl = bar(f.ax(4),xloclf, pdfvallf,1,'stacked','facea',0.6);
%     lfHdl(1).FaceColor = [0.6 1 1];


% returns estimates of normal distribution parameters, mu and sigma
% under 95% confidence interval
    [muhf,sigmahf,muhfCI,sigmahfCI] = normfit((resprophf(i,1:numhf(i)))', 0.05, [], wthf);
    [mulf,sigmalf,mulfCI,sigmalfCI] = normfit((resproplf(i,1:numlf(i)))', 0.05, [], wtlf);
%     text(f.ax(4),0.04,0.91,'d','FontSize',11,'unit','normalized','EdgeColor','k','Margin',2);
%     text(f.ax(4),0.4,0.92,sprintf('LF-HF = %.2f km',offset(i)),'FontSize',11,'unit','normalized');
%     f.ax(4).Box = 'on';
%     grid(f.ax(4), 'on');
%     f.ax(4).GridLineStyle = '--';
%     ymax = f.ax(4).YLim(2)+0.1;
%     ylim(f.ax(4),[0 ymax]);
% %     ylim(f.ax(4),[0 1]);
%     xlabel(f.ax(4),'Residual in prop. (km)','fontsize',11);
%     ylabel(f.ax(4),'PDF estimate','fontsize',11);
%     hold(f.ax(4),'off');
    
    
%     %%% save figure
%     print(f.fig,'-dpdf',strcat(rstpath,'/LZB.',num2str(nfam),'autortmv3.mig.proj.lfit',...
%           num2str(trange(i,1)),'_',num2str(trange(i,2)),'-',num2str(trange(i,3)),...
%           '.',num2str(winlenhf),'_',num2str(winlenlf),'.',num2str(lofflf),'.',...
%           num2str(ccminlf),'.pdf'));
%     print(f.fig,'-dpdf',strcat('/home/data2/chaosong/CurrentResearch/Song_Rubin_2020/Figures/LZB',...
%           num2str(i),'.pdf'));

end

[tmp,tmpind] = sort(sehf(:,2),'descend');
sortsehf = [tmpind tmp];
[tmp,tmpind] = sort(pearwthf(:),'ascend');
sortpearwthf = [tmpind tmp];
[tmp,tmpind] = sort(self(:,2),'descend');
sortself = [tmpind tmp];


%% obtain the vertical distance from lf points to hf line in certain regions

% NOTE here, 2020/09/22, u can choose which vertical distance to use
% vertdist = vertdistwt;
vertdist = vertdistraw;

% SW migs
% ones within angle range and excluding those at distant fams, eg. 002, 010
indsw = setdiff(find(angbest>=220 & angbest<=280),[1]); % [3,5,6,7,9,12,13,15,22,23,24,28,29,33,35,37,38,41,43,44,49,51,52]
% ones have poor lf constraints
indpoorlf = sortself(sortself(:,2)>=3, 1); %[40,4,5,33,3,15,11,23,35,21,45,17,2,31,30,24,32,39,6,9,27,37,51,19,48,26,46,22,8,29,18,47]

% so you can have several options, comment others:  USE OPTION 4!
% % 1. all within angles
% indsw = indsw;  % [3,5,6,7,9,12,13,15,22,23,24,28,29,33,35,37,38,41,43,44,49,51,52]
% % 2. all excluding EITHER poor LF constraints OR there will be very different offsets
% tmp = union(indpoorlf, indlargediff);
% indsw = setdiff(indsw, tmp); % [7,12,13,28,38,41]
% % 3. all excluding ones have large difference in offsets
% indsw = setdiff(indsw, indlargediff);  % [5,6,7,9,12,13,15,22,23,28,33,35,37,38,41,51]
% 4. all excluding ones have large difference in offsets due to LF constraint
tmp = intersect(indlargediff, indpoorlf);
indsw = setdiff(indsw, tmp);   % [5,6,7,9,12,13,15,22,23,28,33,35,37,38,41,43,44,49,51,52]

vdistsw = [];
wtsw = [];
ivdistsw = [];  % dist and weight of independent detections
iwtsw = [];

for i = 1: length(indsw)
    iind = indlfindpen(1:inumlf(indsw(i)), indsw(i));  % index of independent detection in that migration indsw(i)

        temp1 = vertdist(1:numlf(indsw(i)), indsw(i));
%         temp1 = vertdist(1:numlf(indsw(i)), indsw(i))*cos(deg2rad(angbest(indsw(i)) - 250));
        vdistsw = [vdistsw; temp1];
        ivdistsw = [ivdistsw; temp1(iind)];
        
        temp1 = weight(1:numlf(indsw(i)), indsw(i));
        wtsw = [wtsw; temp1];
        iwtsw = [iwtsw; temp1(iind)];
end


% NE migs
% ones within angle range and excluding those at distant fams, eg. 002, 010
indne = setdiff(find(angbest>=40 & angbest<=100),[26,27]); % [2,10,11,18,21,30,34,46,47]
% ones have poor lf constraints
indpoorlf = sortself(sortself(:,2)>=3, 1); %[40,4,5,33,3,15,11,23,35,21,45,17,2,31,30,24,32,39,6,9,27,37,51,19,48,26,46,22,8,29,18,47]
% ones have large difference in offsets due to a different selection of propagation direction
% indlargediff = [2,3,11,14,24,29,32,34,36,40,43,44,49,52];  %
% so you can have several options, comment out others:  USE OPTION 4!
% % 1. all within angles
% indne = indne;  % [2,10,11,18,21,30,34,46,47]
% % 2. all excluding EITHER poor LF constraints OR there will be very different offsets
% tmp = union(indpoorlf, indlargediff);
% indne = setdiff(indne, tmp); % [10]
% % 3. all excluding ones have large difference in offsets
% indne = setdiff(indne, indlargediff);  % [10,18,21,30,46,47]
% 4. all excluding ones have large difference in offsets due to LF constraint
tmp = intersect(indlargediff, indpoorlf);
indne = setdiff(indne, tmp);   % [10,18,21,30,34,46,47]

vdistne = [];
wtne = [];
ivdistne = [];  % dist and weight of independent detections
iwtne = [];
for i = 1: length(indne)
    iind = indlfindpen(1:inumlf(indne(i)), indne(i));  % index of independent detection in that migration indsw(i)
    temp1 = vertdist(1:numlf(indne(i)), indne(i));
%     temp1 = vertdist(1:numlf(indne(i)), indne(i))*cos(deg2rad(angbest(indne(i)) - 70));
    vdistne = [vdistne; temp1];
    ivdistne = [ivdistne; temp1(iind)];
    
    temp1 = weight(1:numlf(indne(i)), indne(i));
    wtne = [wtne; temp1];
    iwtne = [iwtne; temp1(iind)];

end



%%% Change sign of distance to geographical direction
%%% NOTE, 2020/10/27, recall that this 'vertical' distance is LF-HF in the propagation direction,
%%% in which positive means LF is leading, negative is lagging. But this is limited to the
%%% propagation. So if change to absolute geographycal direction, say WSW-ENW, if look at wsw
%%% propagating migrations, if lf leads, then distance should be +, but you need to flip the sign to
%%% convert to geo direction, since lf is in fact at the wsw of hf, - in geo direction. Same if lf
%%% lags. On the other hand, if look at ene group, if lf lags, the raw distance should be -, in
%%% fact, it means lf is at the wsw of hf, still - in the chosen geo coordinates, thus sign of this
%%% group remains unchanged.
%%%
%%% even divided spatially into two parts, for the same propagationn, sign of the offset is the same 

%%% currently focus on main LZB region and 002 region only

% main LZB region 
vdistsw = -vdistsw;     % flip the sign for SW group, check note above
ivdistsw = -ivdistsw;
vdistne = vdistne;     % keep the sign for NE group, check note above
ivdistne = ivdistne;



%%% test the NULL hypothesis that SW and NE group has the same mean, without assuming same variance
%%% NOTE:
%%% test if the possibility of this null hypothesis is higher than a significance level (alpha). If
%%% so, it means the sample mean difference is NOT significant, then fail to reject this null 
%%% hypothesis. Otherwise, if p-value is smaller than alpha, then reject it and regard the
%%% difference is significant.

%%% 1. for SW and NE group at main LZB region, i.e. the NW part of the study region
% a. use significance level 0.05; frequency weights (2), and unequal variance (0)
[testmaina,probmaina,statsmaina] = ttest2_wt(vdistsw, wtsw, vdistne, wtne, 0.05, 2, 0)

% % b. use significance level 0.05; reliability weights (1), and unequal variance (0)
% [testmainb,probmainb,statsmainb] = ttest2_wt(vdistsw, wtsw, vdistne, wtne, 0.05, 1, 0)
% 
% % c. try equal variance (1) as well, since the variances in fact are very close 
% [testmainc,probmainc,statsmainc] = ttest2_wt(vdistsw, wtsw, vdistne, wtne, 0.05, 2, 1)



%%% plot the combined migrations, unshifted all
%%% 1. for SW and NE group at main LZB region, i.e. the NW part of the study region
% indplt = union(indsw,indne);    % index to plot
xran = [-15 10];
yran = [-15 15];
binw = 0.5;
conf = 99;
[f,barsw,pdfxlocsw,pdfvalsw,barne,pdfxlocne,pdfvalne,lgd] = ...
    plt_hist_combined_RTMs(indne,indsw,ranvechf98,medallhf,angbest,xran,yran,vdistsw,wtsw,...
                           vdistne,wtne,binw,conf,'int');
text(f.ax(1),0.5,0.06,'All LF detections','FontSize',12,'unit','normalized',...
     'horizontalalignment','center');
text(f.ax(1),0.5,0.12,strcat({'Excluding '},famex),'FontSize',12,'unit','normalized',...
     'horizontalalignment','center');
xlim(f.ax(2), [-12.5 12.5]); 
xvect = [-7 -12.5+0.5];
yvect = [0.015 0.015];
drawArrow(f.ax(2),xvect,yvect,f.ax(2).XLim,f.ax(2).YLim,'linewidth',1.5,'linestyle','-',...
          'color',[0.5 0.5 0.5]);
xvect = [7 12.5-0.5];
yvect = [0.015 0.015];
drawArrow(f.ax(2),xvect,yvect,f.ax(2).XLim,f.ax(2).YLim,'linewidth',1.5,'linestyle','-',...
          'color',[0.5 0.5 0.5]);
text(f.ax(2),0.05,0.1,strcat({'WSW'}),'fontsize',11,'unit','normalized');
text(f.ax(2),0.95,0.1,strcat({'ENE'}),'fontsize',11,'unit','normalized',...
     'horizontalalignment','right');
lgd.String = [{'ENE propagation'},{'WSW propagation'}];
text(f.ax(3),0.95,0.86,strcat({'10X'}),'fontsize',10,'unit','normalized','EdgeColor','k',...
     'Margin',2,'horizontalalignment','right');
text(f.ax(3),0.05,0.15,strcat({'99% CI'}),'fontsize',12,'unit','normalized','horizontalalignment','left'); 

xlimit = f.ax(2).XLim/10;
xlim(f.ax(3),xlimit); 
% convert data units to global figure units
[xs,ys] = ds2nfu(f.ax(3),xlimit(1),0.2);
[xe,ye] = ds2nfu(f.ax(2),xlimit(1),0.3);
% plot two dashed lines denoting the zoom-in effect
annotation('line',[xs,xe],[ys,ye],'color','k','linestyle','--');
% convert data units to global figure units
[xs,ys] = ds2nfu(f.ax(3),xlimit(2),0.2);
[xe,ye] = ds2nfu(f.ax(2),xlimit(2),0.3);
% plot two dashed lines denoting the zoom-in effect
annotation('line',[xs,xe],[ys,ye],'color','k','linestyle','--');

%%% save figure
print(f.fig,'-dpdf',strcat(rstpath,'/autortmv3.ex',famex,'.vdist_mainSW+NE_all',num2str(nfam),'.pdf'));



%%% test the NULL hypothesis to those independent detections
%%% NOTE:
%%% test if the possibility of this null hypothesis is higher than a significance level (alpha). If
%%% so, it means the sample mean difference is NOT significant, then fail to reject this null 
%%% hypothesis. Otherwise, if p-value is smaller than alpha, then reject it and regard the
%%% difference is significant.

%%% 1. for SW and NE group at main LZB region, i.e. the NW part of the study region
% a. use significance level 0.05; frequency weights (2), and unequal variance (0)
[testmaina,probmaina,statsmaina] = ttest2_wt(ivdistsw, iwtsw, ivdistne, iwtne, 0.05, 2, 0)

% % b. use significance level 0.05; reliability weights (1), and unequal variance (0)
% [testmainb,probmainb,statsmainb] = ttest2_wt(ivdistsw, iwtsw, ivdistne, iwtne, 0.05, 1, 0)
% 
% % c. try equal variance (1) as well, since the variances in fact are very close
% [testmainc,probmainc,statsmainc] = ttest2_wt(ivdistsw, iwtsw, ivdistne, iwtne, 0.05, 2, 1);

% %%% 2. for SW and NE group at fam 002 region
% % a. use significance level 0.05; frequency weights (2), and unequal variance (0)
% [test002a,prob002a,stats002a] = ttest2_wt(ivdist002sw, iwt002sw, ivdist002ne, iwt002ne, 0.05, 2, 0)
% 
% % b. use significance level 0.05; reliability weights (1), and unequal variance (0)
% [test002b,prob002b,stats002b] = ttest2_wt(ivdist002sw, iwt002sw, ivdist002ne, iwt002ne, 0.05, 1, 0)
% 
% % c. try equal variance (1) as well, since the variances in fact are very close
% [test002c,prob002c,stats002c] = ttest2_wt(ivdist002sw, iwt002sw, ivdist002ne, iwt002ne, 0.05, 2, 1)


%%% plot the combined migrations, unshifted independent detections
%%% 1. for SW and NE group at main LZB region, i.e. the NW part of the study region
% indplt = union(indsw,indne);    % index to plot
xran = [-15 10];
yran = [-15 15];
binw = 0.5;
conf = 99;
[f,barisw,pdfxlocisw,pdfvalisw,barine,pdfxlocine,pdfvaline,lgd] = ...
    plt_hist_combined_RTMs(indne,indsw,ranvechf98,medallhf,angbest,xran,yran,ivdistsw,iwtsw,...
                           ivdistne,iwtne,binw,conf,'int');
text(f.ax(1),0.5,0.06,'Independent LF detections','FontSize',12,'unit','normalized',...
     'horizontalalignment','center');
text(f.ax(1),0.5,0.12,strcat({'Excluding '},famex),'FontSize',12,'unit','normalized',...
     'horizontalalignment','center');
xlim(f.ax(2), [-12.5 12.5]); 
xvect = [-7 -12.5+0.5];
yvect = [0.015 0.015];
drawArrow(f.ax(2),xvect,yvect,f.ax(2).XLim,f.ax(2).YLim,'linewidth',1.5,'linestyle','-',...
          'color',[0.5 0.5 0.5]);
% arrow1 = annotation('arrow','linewidth',1.5,'linestyle','-','color',[0.5 0.5 0.5]);
% arrow1.Parent = f.ax(2);
% arrow1.Position = [-7, 0.015, -5, 0] ;
xvect = [7 12.5-0.5];
yvect = [0.015 0.015];
drawArrow(f.ax(2),xvect,yvect,f.ax(2).XLim,f.ax(2).YLim,'linewidth',1.5,'linestyle','-',...
          'color',[0.5 0.5 0.5]);
text(f.ax(2),0.05,0.1,strcat({'WSW'}),'fontsize',11,'unit','normalized');
text(f.ax(2),0.95,0.1,strcat({'ENE'}),'fontsize',11,'unit','normalized',...
     'horizontalalignment','right');
lgd.String = [{'ENE propagation'},{'WSW propagation'}];
text(f.ax(3),0.95,0.86,strcat({'10X'}),'fontsize',10,'unit','normalized','EdgeColor','k',...
     'Margin',2,'horizontalalignment','right');
text(f.ax(3),0.05,0.15,strcat({'99% CI'}),'fontsize',12,'unit','normalized','horizontalalignment','left'); 

xlimit = f.ax(2).XLim/10;
xlim(f.ax(3),xlimit); 
% convert data units to global figure units
[xs,ys] = ds2nfu(f.ax(3),xlimit(1),0.2);
[xe,ye] = ds2nfu(f.ax(2),xlimit(1),0.3);
% plot two dashed lines denoting the zoom-in effect
annotation('line',[xs,xe],[ys,ye],'color','k','linestyle','--');
% convert data units to global figure units
[xs,ys] = ds2nfu(f.ax(3),xlimit(2),0.2);
[xe,ye] = ds2nfu(f.ax(2),xlimit(2),0.3);
% plot two dashed lines denoting the zoom-in effect
annotation('line',[xs,xe],[ys,ye],'color','k','linestyle','--');

print(f.fig,'-dpdf',strcat(rstpath,'/autortmv3.ex',famex,'.vdist_mainSW+NE_ind',num2str(nfam),...
      '.pdf'));
                                      

%% summary of results with error bars from lf fitting
CIFcn = @(x,p)std(x(:),'omitnan')/sqrt(sum(~isnan(x(:)))) * tinv(abs([0,1]-(1-p/100)/2),...
              sum(~isnan(x(:)))-1) + mean(x(:),'omitnan');
f5.fig=figure;
f5.fig.Renderer='Painters';
widin = 8;  % maximum width allowed is 8.5 inches
htin = 6.5;   % maximum height allowed is 11 inches
set(f5.fig,'Position',[scrsz(1)+1*scrsz(3)/10 scrsz(2)+scrsz(4)/20 widin*res htin*res]);
nrow = 1;
ncol = 2;
for isub = 1:nrow*ncol
    f5.ax(isub) = subplot(nrow,ncol,isub);
end

set(f5.ax(1), 'position', [ 0.1, 0.1, 0.4, 0.85]);
set(f5.ax(2), 'position', [ 0.53, 0.1, 0.4, 0.85]);

% index of main LZB region RTMs
% indne = [10,18,21,30,34,46,47];
% indsw = [5,6,7,9,12,13,15,22,23,28,33,35,37,38,41,43,44,49,51,52];
indlzb = [fliplr(reshape(indne,1,[])), fliplr(reshape(indsw,1,[]))];

% index of 002 region RTMs
ind002ne = [21,22];
ind002sw = [4,25];
ind002 = [fliplr(reshape(ind002ne,1,[])), fliplr(reshape(ind002sw,1,[]))];
ind002 = [];

% index of all other RTMs
indall = 1: 1: size(trange,1);
tmp = union(indlzb,ind002);
indelse = setdiff(indall, tmp);
indelse = fliplr(indelse);

ax = f5.ax(1);
hold(ax,'on');
patarea = [0 -100;
           -100 -100;
           -100 100;
            0 100;
            0 -100];
% patch(ax,patarea(:,1),patarea(:,2),'k','Facealpha',0.05,'edgecolor','none');
plot(ax,[0 0],[-100 100],'--','color',[0.6 0.6 0.6],'linew',0.7);

% main LZB region RTMs 
ii = length(indlzb)+ length(ind002);
indplt = indlzb;
for i = 1: length(indplt)
    hfN = numhf(indplt(i));
    
    lfN = numlf(indplt(i));
    lf=resproplf(indplt(i),1:lfN);
    lfm = mean(lf);
    lfstd = std(lf);
    lfsem = lfstd/sqrt(lfN);
    ci95 = tinv([0.025 0.975], lfN-1);
    offcor = bsxfun(@times, lfsem, ci95(:));
    lfci95 = lfm+offcor;
    lfCI95 = CIFcn(lf,95);
    
    ypos = ii-i+1;
    if i > 7
        ypos = ii-i;
    end
    if ~isempty(intersect(indne,indplt(i)))   % NE
        xpos = offset(indplt(i));   % sign unchanged
        e1 = errorbar(ax,xpos,ypos,offcor(1),offcor(2),'horizontal','o','markersize',4.5,'color',...
         'k','linewidth',0.8,'MarkerEdgeColor','k','MarkerFaceColor',[0 1 1],'CapSize',5);
    elseif ~isempty(intersect(indsw,indplt(i))) % SW
        xpos = -offset(indplt(i));   % sign flipped
        e3 = errorbar(ax,xpos,ypos,offcor(1),offcor(2),'horizontal','o','markersize',4.5,'color',...
         'k','linewidth',0.8,'MarkerEdgeColor','k','MarkerFaceColor',[1 96/255 0],'CapSize',5);
    end
    text(ax,xpos+0.1,ypos+0.4,num2str(indplt(i)),'fontsize',7,...
        'HorizontalAlignment','left');
    text(ax,4.9,ypos,strcat(num2str(hfN),'; ',num2str(lfN)),'fontsize',7,...
        'HorizontalAlignment','right');
end
plot(ax,[-100 100],[ii-length(indne) ii-length(indne)],'--','color',[0 0 0],'linew',0.5);
text(ax,0.04,0.96,'a','FontSize',11,'unit','normalized','EdgeColor','k','Margin',2);
text(ax,0.96,0.96,'Oppositely-propagating RTMs','FontSize',10,'unit','normalized',...
     'horizontalalignment','right','EdgeColor','k','Margin',2);

    
% % 002 region RTMs
% ii = ypos-4;
% indplt = ind002;
% for i = 1: length(indplt)
%     hfN = numhf(indplt(i));
%     
%     lfN = numlf(indplt(i));
%     lf=resproplf(indplt(i),1:lfN);
%     lfm = mean(lf);
%     lfstd = std(lf);
%     lfsem = lfstd/sqrt(lfN);
%     ci95 = tinv([0.025 0.975], lfN-1);
%     offcor = bsxfun(@times, lfsem, ci95(:));
%     lfci95 = lfm+offcor;
%     lfCI95 = CIFcn(lf,95);
%     
%     ypos = ii-i+1;
%     if i == 2
%         ypos = -1.5;
%     elseif i == 3
%         ypos = -3.5;
%     elseif i == 4
%         ypos = -5;
%     end
%     if angbest(indplt(i))>=0 && angbest(indplt(i))<90       % NE
%         xpos = offset(indplt(i));   % sign unchanged
%         e1 = errorbar(ax,xpos,ypos,offcor(1),offcor(2),'horizontal','o','markersize',4.5,'color',...
%          'k','linewidth',0.8,'MarkerEdgeColor','k','MarkerFaceColor',[0 1 1],'CapSize',5);
%     elseif angbest(indplt(i))>=180 && angbest(indplt(i))<270    % SW
%         xpos = -offset(indplt(i));   % sign flipped
%         e3 = errorbar(ax,xpos,ypos,offcor(1),offcor(2),'horizontal','o','markersize',4.5,'color',...
%          'k','linewidth',0.8,'MarkerEdgeColor','k','MarkerFaceColor',[1 96/255 0],'CapSize',5);
%     elseif angbest(indplt(i))>=90 && angbest(indplt(i))<180     % SE
%         xpos = offset(indplt(i));   % sign unchanged
%         e2 = errorbar(ax,xpos,ypos,offcor(1),offcor(2),'horizontal','o','markersize',4.5,'color',...
%          'k','linewidth',0.8,'MarkerEdgeColor','k','MarkerFaceColor','k','CapSize',5);
%     end
%     text(ax,xpos+0.1,ypos+0.4,num2str(indplt(i)),'fontsize',7,...
%         'HorizontalAlignment','left');
%     text(ax,4.9,ypos,strcat(num2str(hfN),'; ',num2str(lfN)),'fontsize',7,...
%         'HorizontalAlignment','right');
% end
% scatter(ax,-0.66, -0.5, 12,'MarkerEdgeColor','k','linewidth',1);
% scatter(ax,-1.27, -0.5, 12,'MarkerEdgeColor','k','linewidth',1);
% plot(ax,[-0.66 -1.13], [-0.5 -0.5], 'k:','linewidth',1.5);
% scatter(ax, 0.22, -4, 12,'MarkerEdgeColor','k','linewidth',1);
% scatter(ax,-0.19, -4, 12,'MarkerEdgeColor','k','linewidth',1);
% plot(ax,[0.22 -0.19], [-4 -4], 'k:','linewidth',1.5);
% plot(ax,[-100 100],[ii-length(ind002ne)-0.5 ii-length(ind002ne)-0.5],'--','color',[0 0 0],'linew',0.5);
% text(ax,0.75,0.3,'Fam 002 region','FontSize',10,'unit','normalized',...
%      'horizontalalignment','center','EdgeColor','k','Margin',2);

% for excluding fam 017
yran = [-3,32];
% for excluding fam 068

% for excluding fam 099

% for excluding fam 043

xran = [-5,5];
% axis(ax,[xran yran]);
xvect = [-1.5 -4.5];
yvect = [-2 -2];
drawArrow(ax,xvect,yvect,xran,yran,'linewidth',1.5,'linestyle','-',...
          'color',[0.5 0.5 0.5]);
xvect = [1.5 4.5];
yvect = [-2 -2];
drawArrow(ax,xvect,yvect,xran,yran,'linewidth',1.5,'linestyle','-',...
          'color',[0.5 0.5 0.5]);
text(ax,0.2,0.05,strcat({'WSW'}),'fontsize',10,'fontw','bold','unit','normalized',...
     'horizontalalignment','center');
text(ax,0.8,0.05,strcat({'ENE'}),'fontsize',10,'fontw','bold','unit','normalized',...
     'horizontalalignment','center');
 
% drawArrow(ax,[-4.2 -1],[5 5],xran,yran,'linewidth',1,'linestyle','-','color',[0.65 0.65 0.65]);
% text(ax,-4.8,5,'I','FontSize',13,'VerticalAlignment','middle'); 
% drawArrow(ax,[-4.2 -1],[13 13],xran,yran,'linewidth',1,'linestyle','-','color',[0.65 0.65 0.65]);
% text(ax,-4.8,13,'II','FontSize',13,'VerticalAlignment','middle');
% drawArrow(ax,[-4.2 -1],[16 16],xran,yran,'linewidth',1,'linestyle','-','color',[0.65 0.65 0.65]);
% text(ax,-4.8,16,'III','FontSize',13,'VerticalAlignment','middle');
% drawArrow(ax,[-4.2 -1],[21 21],xran,yran,'linewidth',1,'linestyle','-','color',[0.65 0.65 0.65]);
% text(ax,-4.8,21,'IV','FontSize',13,'VerticalAlignment','middle');

 
xlabel(ax,'LF-HF (km)','fontsize',11);
% ylabel(ax,'RTM number in chronological order','fontsize',11);
legend(ax,[e1, e3],{'ENE propagation','WSW propagation'},'FontSize',7,...
       'Position',[0.156 0.63 0.065 0.05]);       % ,'edgecolor',[0.5 0.5 0.5]
xticks(ax,-5:1:5);
ax.YTick = [];
ax.Box='on';
grid(ax,'on');
set(ax,'GridLineStyle','--');
ax.FontSize = 9;
hold(ax,'off');


% Other RTMs
ax = f5.ax(2);
hold(ax,'on');
patarea = [0 -100;
           -100 -100;
           -100 100;
            0 100;
            0 -100];
patch(ax,patarea(:,1),patarea(:,2),'k','Facealpha',0.05,'edgecolor','none');
ii = length(indlzb)+ length(ind002);
indplt = indelse;
for i = 1: length(indplt)
    hfN = numhf(indplt(i));
    
    lfN = numlf(indplt(i));
    lf=resproplf(indplt(i),1:lfN);
    lfm = mean(lf);
    lfstd = std(lf);
    lfsem = lfstd/sqrt(lfN);
    ci95 = tinv([0.025 0.975], lfN-1);
    offcor = bsxfun(@times, lfsem, ci95(:));
    lfci95 = lfm+offcor;
    lfCI95 = CIFcn(lf,95);
    
    ypos = ii-i+1;
    if angbest(indplt(i))>0 && angbest(indplt(i))<=90       % NE
        xpos = offset(indplt(i));   % sign unchanged
        e1 = errorbar(ax,xpos,ypos,offcor(1),offcor(2),'horizontal','o','markersize',4.5,'color',...
         'k','linewidth',0.8,'MarkerEdgeColor','k','MarkerFaceColor',[0 1 1],'CapSize',5);
    elseif angbest(indplt(i))>90 && angbest(indplt(i))<=180     % SE
        xpos = offset(indplt(i));   % sign unchanged
        e2 = errorbar(ax,xpos,ypos,offcor(1),offcor(2),'horizontal','o','markersize',4.5,'color',...
         'k','linewidth',0.8,'MarkerEdgeColor','k','MarkerFaceColor','k','CapSize',5);
    elseif angbest(indplt(i))>180 && angbest(indplt(i))<=270    % SW
        xpos = offset(indplt(i));   % sign unchanged
        e3 = errorbar(ax,xpos,ypos,offcor(1),offcor(2),'horizontal','o','markersize',4.5,'color',...
         'k','linewidth',0.8,'MarkerEdgeColor','k','MarkerFaceColor',[1 96/255 0],'CapSize',5);
    else
        xpos = offset(indplt(i));   % sign unchanged
        e4 = errorbar(ax,xpos,ypos,offcor(1),offcor(2),'horizontal','o','markersize',4.5,'color',...
         'k','linewidth',0.8,'MarkerEdgeColor','k','MarkerFaceColor',[0.5 0.5 0.5],'CapSize',5);        
    end
    text(ax,xpos+0.1,ypos+0.4,num2str(indplt(i)),'fontsize',7,...
        'HorizontalAlignment','left');
    text(ax,4.9,ypos,strcat(num2str(hfN),'; ',num2str(lfN)),'fontsize',7,...
        'HorizontalAlignment','right');
end

text(ax,0.04,0.96,'b','FontSize',11,'unit','normalized','EdgeColor','k','Margin',2);
text(ax,0.96,0.96,'Other RTMs','FontSize',10,'unit','normalized',...
     'horizontalalignment','right','EdgeColor','k','Margin',2);
% text(ax,0.7,0.97,'LZB','FontSize',11,'unit','normalized','EdgeColor','k','Margin',2);
text(ax,0.2,0.07,'LF','FontSize',10,'unit','normalized','HorizontalAlignment','center','fontw',...
    'bold');
text(ax,0.2,0.04,'lags','FontSize',10,'unit','normalized','HorizontalAlignment','center','fontw',...
    'bold');
text(ax,0.8,0.07,'LF','FontSize',10,'unit','normalized','HorizontalAlignment','center','fontw',...
    'bold');
text(ax,0.8,0.04,'leads','FontSize',10,'unit','normalized','HorizontalAlignment','center','fontw',...
    'bold');
% for excluding fam 017
yran = [2 31];
% for excluding fam 068

% for excluding fam 099

% for excluding fam 043

xran = [-5 5];
% drawArrow(ax,[-4.2 -2.5],[5 5],xran,yran,'linewidth',1.5,'linestyle','-','color',[0.65 0.65 0.65]);
% text(ax,-4.8,5,'IV','FontSize',13,'VerticalAlignment','middle','backgroundcolor','w',...
%      'Margin',0.5);
axis(ax,[xran yran]);
xlabel(ax,'LF-HF (km)','fontsize',11);
% ylabel(ax,'RTM number in chronological order','fontsize',11);
legend(ax,[e1, e2, e3, e4],{'NE propagation','SE propagation','SW propagation', 'NW propagation'},...
       'FontSize',7,'Position',[0.58 0.65 0.065 0.05]);  
xticks(ax,-5:1:5);
ax.YTick = [];
% yticks(0:1:8);
ax.Box='on';
grid(ax,'on');
set(ax,'GridLineStyle','--');
ax.FontSize = 9;
hold(ax,'off');

%%% save figure
print(f5.fig,'-dpdf',strcat(rstpath,'/LZB.ex',famex,'.autortmv3.mig.lfit.sum.',num2str(winlenhf),'_',...
    num2str(winlenlf),'.',num2str(lofflf),'.',num2str(ccminlf),'.pdf'));















    
    
    
    