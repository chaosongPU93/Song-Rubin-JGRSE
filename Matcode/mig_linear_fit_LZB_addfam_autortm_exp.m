% function mig_linear_fit_LZB_addfam_autortm_v3_exp
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script is used to plot an example of migration in which LF starts earlier
% than HF, the corresponding serial number in 'mig_linear_fit_LZB_addfam_autortm_v3'
% is #13. In that script, the time range was selected automatically based on the 
% tremor bursts determination algorithm to HF detections. Here we plot the longer
% range inlcuding the earlier LF part.
%
% The reason behind this migration is NOT because the portion of fault that hosts
% the earlier LF detections is not producing HF. In fact it produces HF at other
% time. The actual reason is that HF is coming from a broader source region so that
% the coherency is worse at earlier time which is lower than our detecting threshold.
% BUT this point has to be confirmed by trying to lower the threshold.
%
% This is created for the migration example #13, since we want to lower the 
% the HF detection CC threshold to redo the detection for day 2004 198 (Jul 14),
% in order to see if during the earlier time we can see HF. This is to answer
% why we missed HF in the current detection regime (2021/01/26)
%
%
%    2020/09/10, i am adding a fam 001 back again for last try, 13 fam in total now
%                this is for adding new fams 006, 001
%
% Chao Song, chaosong@princeton.edu
% First created date:   2021/01/25
% Last modified date:   2021/01/25
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Initialization
format short e   % Set the format to 5-digit floating point
clear
close all
clc

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
            '017';
            '006';
            '001';
%             '158';      % 158, 20200916,testing purpose
            ];     % 2020/09/07, i am adding a new family 006 to the pool, final version
        
nfam = size(nfampool,1);
disp(nfam); 

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
fname1 = strcat(rstpath, '/evtloc.all',num2str(nfam),'fam.dcutnodou.',SUFFIXhf);
fname2 = strcat(rstpath, '/evtloc.all',num2str(nfam),'fam.nodcutnodou.',SUFFIXhf);
fname3 = strcat(rstpath, '/evtloc.all',num2str(nfam),'fam.eq8kmdcutnodou.',SUFFIXhf);
fname4 = strcat(rstpath, '/evtloc.all',num2str(nfam),'fam.eq12kmdcutnodou.',SUFFIXhf);

fname = fname1;
hftime = load(fname);
% 25 cols, format is:
%   E(043) N(043) E(own) N(own) dist(own) lon lat dep off12 off13 off12sec off13sec date
%   main_arrival_time  cent_of_win  avecc ccadd1 ccadd2 ccadd3 ccadd4 ampmax ampmin eng normeng famnum
%


SUFFIXlf = strcat('lf.time.',num2str(winlenhf),'_',num2str(winlenlf),'.',num2str(lofflf),...
                  '.',num2str(ccminlf));
fname1 = strcat(rstpath, '/evtloc.all',num2str(nfam),'fam.dcutnodou.',SUFFIXlf);
fname2 = strcat(rstpath, '/evtloc.all',num2str(nfam),'fam.nodcutnodou.',SUFFIXlf);
fname3 = strcat(rstpath, '/evtloc.all',num2str(nfam),'fam.eq8kmdcutnodou.',SUFFIXlf);
fname4 = strcat(rstpath, '/evtloc.all',num2str(nfam),'fam.eq12kmdcutnodou.',SUFFIXlf);

fname = fname1;
lftime = load(fname);

% this is inverted from (0,0) of all fams, same order, location of control points
loccont = [-123.492667 48.451500 38.1400; 
           -123.772167 48.493000 35.5900; 
           -123.863167 48.528167 35.2100;
           -123.603333 48.440167 36.7100;
           -123.800167 48.408833 34.5200;
           -123.893333 48.536500 35.0700;
           -123.864500 48.498667 34.8800;
           -123.753333 48.525667 36.2000;
           -123.703667 48.502667 36.4100;
           -123.814333 48.538667 35.7900;
           -123.838500 48.544833 35.6600;
           -123.908000 48.494167 34.5100;       % 006
           -123.879667 48.446167 34.2600;       % 001
%            -123.967000 48.458667 33.7800;       % 158, 20200916,testing purpose
           ];
       

relacont = loccont;
[dx, dy] = absloc2relaloc(loccont(:,1),loccont(:,2),loccont(2,1),loccont(2,2));
relacont(:,1) = dx;
relacont(:,2) = dy;

%%% load the location resolution points, +-1 & 2 samples
reslzb1 = load(fullfile(rstpath,'evtloc.13fam_locres_hf1spl'));
reslzb2 = load(fullfile(rstpath,'evtloc.13fam_locres_hf2spl'));

% sort according to day, sec
hftime = sortrows(hftime, [13, 15]);
lftime = sortrows(lftime, [13, 15]);

% this is new time ranges that come from 'identify_RTMs_v3.m'
trange = [
    2004198   8.3520e+04   8.6400e+04   % combine, y
    ]; 
    
% hardcopy of angbest from either min. rmse or max. slope of accepted ones ntranacc
angbest = [
   240
   ];

indspeed = [
    2
    27
    34
    37
    47
    ];



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

% choose the end of RTM 42 as the division point, get its location x y
x0 = -3.5474e+00;
y0 = 2.0923e+00;
% the division line passes through this point and orientates sse
xdiv = -12:0.01:0;
a1 = 1/tan(deg2rad(180-22.5));
b1 = y0-a1*x0;
ydiv = linefcn(xdiv,a1,b1);

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
       
         
for i = 1: size(trange,1)
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
%     % only use the detections ouside the bound
%     [is,ion] = inpolygon(mighfdum(:,1),mighfdum(:,2),reg002(:,1),reg002(:,2));
%     isinreg = is | ion;
%     mighfdum = mighfdum(isinreg ~= 1, :);
%     
%     [is,ion] = inpolygon(miglfdum(:,1),miglfdum(:,2),reg002(:,1),reg002(:,2));
%     isinreg = is | ion;
%     miglfdum = miglfdum(isinreg ~= 1, :);
    
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

    
    % assign the location resolution from one fam
%     ifam = find(famperchf(i,:) == max(famperchf(i,:))); % choose fam contributes most detections
    ifam = 2;   % choose fam 043 as the representative
    tmp1 = reslzb1((ifam-1)*4+1: ifam*4, :);
    [dx, dy] = absloc2relaloc(tmp1(:,1),tmp1(:,2),loccont(ifam,1),loccont(ifam,2));
    vecreslzb1 = [dx dy];
    tmp2 = reslzb2((ifam-1)*4+1: ifam*4, :);
    [dx, dy] = absloc2relaloc(tmp2(:,1),tmp2(:,2),loccont(ifam,1),loccont(ifam,2));
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
    text(f.ax(1),0.58,0.1,strcat('I', {', including earlier time'}),'FontSize',11,'unit','normalized',...
         'horizontalalignment','center');
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
    text(f.ax(2),0.58,0.1,strcat('I', {', including earlier time'}),'FontSize',11,'unit','normalized',...
         'horizontalalignment','center');
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
    
    
    %%% save figure
    if strcmp(fname, fname1)
        print(f.fig,'-dpdf',strcat(rstpath,'/LZB.',num2str(nfam),'autortmv3.sel1long.',...
            num2str(trange(i,1)),'_',num2str(trange(i,2)),'-',num2str(trange(i,3)),...
            '.',num2str(winlenhf),'_',num2str(winlenlf),'.',num2str(lofflf),'.',...
            num2str(ccminlf),'.pdf'));
    elseif strcmp(fname, fname2)
        print(f.fig,'-dpdf',strcat(rstpath,'/LZB.',num2str(nfam),'autortmv3.sel1long.nodcut.',...
            num2str(trange(i,1)),'_',num2str(trange(i,2)),'-',num2str(trange(i,3)),...
            '.',num2str(winlenhf),'_',num2str(winlenlf),'.',num2str(lofflf),'.',...
            num2str(ccminlf),'.pdf'));
    elseif strcmp(fname, fname3)
        print(f.fig,'-dpdf',strcat(rstpath,'/LZB.',num2str(nfam),'autortmv3.sel1long.eq8kmdcut.',...
            num2str(trange(i,1)),'_',num2str(trange(i,2)),'-',num2str(trange(i,3)),...
            '.',num2str(winlenhf),'_',num2str(winlenlf),'.',num2str(lofflf),'.',...
            num2str(ccminlf),'.pdf'));
    elseif strcmp(fname, fname4)
        print(f.fig,'-dpdf',strcat(rstpath,'/LZB.',num2str(nfam),'autortmv3.sel1long.eq12kmdcut.',...
            num2str(trange(i,1)),'_',num2str(trange(i,2)),'-',num2str(trange(i,3)),...
            '.',num2str(winlenhf),'_',num2str(winlenlf),'.',num2str(lofflf),'.',...
            num2str(ccminlf),'.pdf'));        
    end
        
%     print(f.fig,'-dpdf',strcat('/home/data2/chaosong/CurrentResearch/Song_Rubin_2020/Figures/LZB',...
%           num2str(i),'.pdf'));

end









    
    
    
    