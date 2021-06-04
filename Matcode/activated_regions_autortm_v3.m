% function activated_regions_autortm_v3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script is to use the density plot to lighten up the activated regions
% during the ETS periods in hf and lf
% Use the new automatic time ranges of migrations that come from 'identify_RTMs_v3.m'
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
fname = strcat(rstpath, '/evtloc.all',num2str(nfam),'fam.dcutnodou.',SUFFIXhf);
hftime = load(fname);
% 25 cols, format is:
%   E(043) N(043) E(own) N(own) dist(own) lon lat dep off12 off13 off12sec off13sec date
%   main_arrival_time  cent_of_win  avecc ccadd1 ccadd2 ccadd3 ccadd4 ampmax ampmin eng normeng famnum
%


SUFFIXlf = strcat('lf.time.',num2str(winlenhf),'_',num2str(winlenlf),'.',num2str(lofflf),...
                  '.',num2str(ccminlf));
fname = strcat(rstpath, '/evtloc.all',num2str(nfam),'fam.dcutnodou.',SUFFIXlf);
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


% sort according to day, sec
hftime = sortrows(hftime, [13, 15]);
lftime = sortrows(lftime, [13, 15]);

%% 
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
    
% hardcopy of angbest from either min. rmse or max. slope of accepted ones ntranacc
angbest = [
   250
    90
   220
   185
   225
   270
   235
   200
   230
    80
    65
   235
   245
   310
   250
   120
   195
    90
   140
   130
    70
   255
   250
   260
   120
   110
    80
   240
   220
    80
   115
   105
   260
    95
   260
   130
   230
   250
   195
   115
   250
   165
   245
   235
    10
    95
    85
   185
   255
   115
   250
   235
   ];

indspeed = [
    2
    27
    34
    37
    47
    ];

%%
% all detections in RTMs, lf and hf, that are located at main lzb region in each propagation group
mighfa = [];
miglfa = [];


for i = 1: length(trange)

%     disp(trange(i,:));
    indhf = find(hftime(:,13)==trange(i,1) & hftime(:,15)>=trange(i,2) & ...
                 hftime(:,15)<=trange(i,3));
    mighf = hftime(indhf,:);
    
    mighfdum = mighf;
    for j = 1: size(mighf,1)
        x0 = mighfdum(j,1);
        y0 = mighfdum(j,2);
        [newx,newy] = coordinate_rot(x0,y0,-(angbest(i)-90),0,0);
        mighfdum(j,1) = newx;
        mighfdum(j,2) = newy;
    end
    
    indlf = find(lftime(:,13)==trange(i,1) & lftime(:,15)>=trange(i,2) & ...
                 lftime(:,15)<=trange(i,3));
    miglf = lftime(indlf,:); 
    
    miglfdum = miglf;
    for j = 1: size(miglf,1)
        x0 = miglfdum(j,1);
        y0 = miglfdum(j,2);
        [newx,newy] = coordinate_rot(x0,y0,-(angbest(i)-90),0,0);
        miglfdum(j,1) = newx;
        miglfdum(j,2) = newy;
    end
    
    
    mintime = max(min(mighfdum(:,15)), min(miglfdum(:,15)));
    maxtime = min(max(mighfdum(:,15)), max(miglfdum(:,15)));
    mighfdum2 = mighfdum(mighfdum(:,15)>=mintime & mighfdum(:,15)<=maxtime, :);
    miglfdum2 = miglfdum(miglfdum(:,15)>=mintime & miglfdum(:,15)<=maxtime, :);

        
    % rotate back the projected locations to get the original coordinates
    tmp = mighfdum2;
    [newx,newy] = coordinate_rot(tmp(:,1),tmp(:,2),(angbest(i)-90),0,0);
    tmp(:,1) = newx;
    tmp(:,2) = newy;
    mighfa = [mighfa; tmp];
    tmp = miglfdum2;
    [newx,newy] = coordinate_rot(tmp(:,1),tmp(:,2),(angbest(i)-90),0,0);
    tmp(:,1) = newx;
    tmp(:,2) = newy;
    miglfa = [miglfa; tmp];
   
    
end

%% density in each ETS, 2004 and 2005
f1.fig=figure;
f1.fig.Renderer='Painters';
widin = 8;  % maximum width allowed is 8.5 inches
htin = 9;   % maximum height allowed is 11 inches
set(f1.fig,'Position',[scrsz(1)+1*scrsz(3)/10 scrsz(2)+scrsz(4)/20 widin*res htin*res]);
nrow = 2;
ncol = 2;
for isub = 1:nrow*ncol
    f1.ax(isub) = subplot(nrow,ncol,isub);
end

%%% reposition
set(f1.ax(1), 'position', [ 0.08, 0.6, 0.36, 0.32]);
set(f1.ax(2), 'position', [ 0.52, 0.6, 0.36, 0.32]);
set(f1.ax(3), 'position', [ 0.08, 0.18, 0.36, 0.32]);
set(f1.ax(4), 'position', [ 0.52, 0.18, 0.36, 0.32]);
    

xran = [-20 25];
yran = [-20 20];

msizehf = 4;
msizelf = 10;

% subplot 1 of figure i
hfplt = hftime(hftime(:,13) > 2004*1000 & hftime(:,13) < 2005*1000, :);
lfplt = lftime(lftime(:,13) > 2004*1000 & lftime(:,13) < 2005*1000, :);

ax = f1.ax(1);
hold(ax,'on');
plot(ax,[-100 100],[0 0],'k--');
plot(ax,[0 0],[-100 100],'k--');
ax.FontSize = 9;
% create a density matrix to store the number of detections in each small grid
dxhf = 0.2;
dyhf = 0.2;
[density1d,xloc2d,yloc2d,density2d] = density_matrix(hfplt(:,1),hfplt(:,2),...
    xran,yran,dxhf,dyhf);
dumhf = density1d(density1d(:,3)>0, :);
dumhf(dumhf(:,3)>1, :) = [];
scatter(ax,dumhf(:,1),dumhf(:,2), msizehf, log10(dumhf(:,3)),'s','linew',0.2);  %, 'MarkerEdgeColor', 'w')
dumhf = sortrows(density1d(density1d(:,3)>0, :), 3);
dumhf(dumhf(:,3)==1, :) = [];
scatter(ax,dumhf(:,1),dumhf(:,2), msizehf, log10(dumhf(:,3)), 'filled','s');  %, 'MarkerEdgeColor', 'w')
% imagesc(ax,[xran(1)+0.5*dx xran(end)-0.5*dx], [yran(1)+0.5*dy yran(end)-0.5*dy], density2d');
oldcmap = colormap(ax,'jet');
% colormap(ax, flipud(oldcmap) );
c=colorbar(ax,'SouthOutside');
pos = ax.Position;
c.Position = [pos(1), pos(2)-0.018, pos(3), 0.02];
c.Label.String = strcat('log_{10}(N) of detections');
c.Label.FontSize = 11;
text(ax, 0.85, 0.93, '2004','FontSize',12,'unit','normalized','horizontalalignment','center',...
     'EdgeColor','k','Margin',2);
text(ax,0.85,0.8,'HF','FontSize',12,'unit','normalized','horizontalalignment','center');
text(ax,0.04,0.93,'a','FontSize',11,'unit','normalized','EdgeColor','k','Margin',2);
text(ax,0.5,0.93,strcat(num2str(length(hfplt(:,1))),{' detections'}),'FontSize',10,'unit','normalized',...
     'horizontalalignment','center');
ax.Box = 'on';
grid(ax, 'on');
axis(ax, 'equal');
ax.GridLineStyle = '--';
ax.XAxisLocation = 'top';
xlim(ax,xran);
ylim(ax,yran);
xticks(ax,xran(1):5:xran(2));
yticks(ax,yran(1):5:yran(2));
xlabel(ax,'E (km)','fontsize',11);
ylabel(ax,'N (km)','fontsize',11);
hold(ax,'off');

% subplot 2 of figure i
ax = f1.ax(2);
hold(ax,'on');
plot(ax,[-100 100],[0 0],'k--');
plot(ax,[0 0],[-100 100],'k--');
ax.FontSize = 9;
% create a density matrix to store the number of detections in each small grid
dxlf = 0.5;
dylf = 0.5;
[density1d,xloc2d,yloc2d,density2d] = density_matrix(lfplt(:,1),lfplt(:,2),...
    xran,yran,dxlf,dylf);
dumhf = density1d(density1d(:,3)>0, :);
dumhf(dumhf(:,3)>1, :) = [];
scatter(ax,dumhf(:,1),dumhf(:,2), msizelf, log10(dumhf(:,3)),'s','linew',0.2);  %, 'MarkerEdgeColor', 'w')
dumhf = sortrows(density1d(density1d(:,3)>0, :), 3);
dumhf(dumhf(:,3)==1, :) = [];
scatter(ax,dumhf(:,1),dumhf(:,2), msizelf, log10(dumhf(:,3)), 'filled','s');  %, 'MarkerEdgeColor', 'w')
% imagesc(ax,[xran(1)+0.5*dx xran(end)-0.5*dx], [yran(1)+0.5*dy yran(end)-0.5*dy], density2d');
oldcmap = colormap(ax,'jet');
% colormap(ax, flipud(oldcmap) );
c=colorbar(ax,'SouthOutside');
pos = ax.Position;
c.Position = [pos(1), pos(2)-0.018, pos(3), 0.02];
c.Label.String = strcat('log_{10}(N) of detections');
c.Label.FontSize = 11;
text(ax, 0.85, 0.93, '2004','FontSize',12,'unit','normalized','horizontalalignment','center',...
     'EdgeColor','k','Margin',2);
text(ax,0.85,0.8,'LF','FontSize',12,'unit','normalized','horizontalalignment','center');
text(ax,0.04,0.93,'b','FontSize',11,'unit','normalized','EdgeColor','k','Margin',2);
text(ax,0.5,0.93,strcat(num2str(length(lfplt(:,1))),{' detections'}),'FontSize',10,'unit','normalized',...
     'horizontalalignment','center');
ax.Box = 'on';
grid(ax, 'on');
axis(ax, 'equal');
ax.GridLineStyle = '--';
ax.XAxisLocation = 'top';
xlim(ax,xran);
ylim(ax,yran);
xticks(ax,xran(1):5:xran(2));
yticks(ax,yran(1):5:yran(2));
xlabel(ax,'E (km)','fontsize',11);
ylabel(ax,'N (km)','fontsize',11);
hold(ax,'off');


% subplot 3 of figure i
hfplt = hftime(hftime(:,13) > 2005*1000, :);
lfplt = lftime(lftime(:,13) > 2005*1000, :);

ax = f1.ax(3);
hold(ax,'on');
plot(ax,[-100 100],[0 0],'k--');
plot(ax,[0 0],[-100 100],'k--');
ax.FontSize = 9;
% create a density matrix to store the number of detections in each small grid
dxhf = 0.2;
dyhf = 0.2;
[density1d,xloc2d,yloc2d,density2d] = density_matrix(hfplt(:,1),hfplt(:,2),...
    xran,yran,dxhf,dyhf);
dumhf = density1d(density1d(:,3)>0, :);
dumhf(dumhf(:,3)>1, :) = [];
scatter(ax,dumhf(:,1),dumhf(:,2), msizehf, log10(dumhf(:,3)),'s','linew',0.2);  %, 'MarkerEdgeColor', 'w')
dumhf = sortrows(density1d(density1d(:,3)>0, :), 3);
dumhf(dumhf(:,3)==1, :) = [];
scatter(ax,dumhf(:,1),dumhf(:,2), msizehf, log10(dumhf(:,3)), 'filled','s');  %, 'MarkerEdgeColor', 'w')
% imagesc(ax,[xran(1)+0.5*dx xran(end)-0.5*dx], [yran(1)+0.5*dy yran(end)-0.5*dy], density2d');
oldcmap = colormap(ax,'jet');
% colormap(ax, flipud(oldcmap) );
c=colorbar(ax,'SouthOutside');
pos = ax.Position;
c.Position = [pos(1), pos(2)-0.018, pos(3), 0.02];
c.Label.String = strcat('log_{10}(N) of detections');
c.Label.FontSize = 11;
text(ax, 0.85, 0.93, '2005','FontSize',12,'unit','normalized','horizontalalignment','center',...
     'EdgeColor','k','Margin',2);
text(ax,0.85,0.8,'HF','FontSize',12,'unit','normalized','horizontalalignment','center');
text(ax,0.04,0.93,'c','FontSize',11,'unit','normalized','EdgeColor','k','Margin',2);
text(ax,0.5,0.93,strcat(num2str(length(hfplt(:,1))),{' detections'}),'FontSize',10,'unit','normalized',...
     'horizontalalignment','center');
ax.Box = 'on';
grid(ax, 'on');
axis(ax, 'equal');
ax.GridLineStyle = '--';
ax.XAxisLocation = 'top';
xlim(ax,xran);
ylim(ax,yran);
xticks(ax,xran(1):5:xran(2));
yticks(ax,yran(1):5:yran(2));
xlabel(ax,'E (km)','fontsize',11);
ylabel(ax,'N (km)','fontsize',11);
hold(ax,'off');

% subplot 4 of figure i
ax = f1.ax(4);
hold(ax,'on');
plot(ax,[-100 100],[0 0],'k--');
plot(ax,[0 0],[-100 100],'k--');
ax.FontSize = 9;

% create a density matrix to store the number of detections in each small grid
dxlf = 0.5;
dylf = 0.5;
[density1d,xloc2d,yloc2d,density2d] = density_matrix(lfplt(:,1),lfplt(:,2),...
    xran,yran,dxlf,dylf);
dumhf = density1d(density1d(:,3)>0, :);
dumhf(dumhf(:,3)>1, :) = [];
scatter(ax,dumhf(:,1),dumhf(:,2), msizelf, log10(dumhf(:,3)),'s','linew',0.2);  %, 'MarkerEdgeColor', 'w')
dumhf = sortrows(density1d(density1d(:,3)>0, :), 3);
dumhf(dumhf(:,3)==1, :) = [];
scatter(ax,dumhf(:,1),dumhf(:,2), msizelf, log10(dumhf(:,3)), 'filled','s');  %, 'MarkerEdgeColor', 'w')
% imagesc(ax,[xran(1)+0.5*dx xran(end)-0.5*dx], [yran(1)+0.5*dy yran(end)-0.5*dy], density2d');
oldcmap = colormap(ax,'jet');
% colormap(ax, flipud(oldcmap) );
c=colorbar(ax,'SouthOutside');
pos = ax.Position;
c.Position = [pos(1), pos(2)-0.018, pos(3), 0.02];
c.Label.String = strcat('log_{10}(N) of detections');
c.Label.FontSize = 11;
text(ax, 0.85, 0.93, '2005','FontSize',12,'unit','normalized','horizontalalignment','center',...
     'EdgeColor','k','Margin',2);
text(ax,0.85,0.8,'LF','FontSize',12,'unit','normalized','horizontalalignment','center');
text(ax,0.04,0.93,'d','FontSize',11,'unit','normalized','EdgeColor','k','Margin',2);
text(ax,0.5,0.93,strcat(num2str(length(lfplt(:,1))),{' detections'}),'FontSize',10,'unit','normalized',...
     'horizontalalignment','center');
ax.Box = 'on';
grid(ax, 'on');
axis(ax, 'equal');
ax.GridLineStyle = '--';
ax.XAxisLocation = 'top';
xlim(ax,xran);
ylim(ax,yran);
xticks(ax,xran(1):5:xran(2));
yticks(ax,yran(1):5:yran(2));
xlabel(ax,'E (km)','fontsize',11);
ylabel(ax,'N (km)','fontsize',11);
hold(ax,'off');

print(f1.fig,'-dpdf',strcat(rstpath,'/autortmv3_density_ETS2004_2005.pdf'));


%%
f2.fig=figure;
f2.fig.Renderer='Painters';
widin = 8;  % maximum width allowed is 8.5 inches
htin = 9;   % maximum height allowed is 11 inches
set(f2.fig,'Position',[scrsz(1)+1*scrsz(3)/10 scrsz(2)+scrsz(4)/20 widin*res htin*res]);
nrow = 2;
ncol = 2;
for isub = 1:nrow*ncol
    f2.ax(isub) = subplot(nrow,ncol,isub);
end

%%% reposition
set(f2.ax(1), 'position', [ 0.08, 0.6, 0.36, 0.32]);
set(f2.ax(2), 'position', [ 0.52, 0.6, 0.36, 0.32]);
set(f2.ax(3), 'position', [ 0.08, 0.18, 0.36, 0.32]);
set(f2.ax(4), 'position', [ 0.52, 0.18, 0.36, 0.32]);

% xran = [-14 9];
% yran = [-9 13];

msizehf = 4;
msizelf = 10;

% subplot 1 of figure i
hfplt = hftime(hftime(:,13) > 2004*1000, :);
lfplt = lftime(lftime(:,13) > 2004*1000, :);

ax = f2.ax(1);
hold(ax,'on');
plot(ax,[-100 100],[0 0],'k--');
plot(ax,[0 0],[-100 100],'k--');
ax.FontSize = 9;
%     ind = find(density1d>=0);
%     scatter(ax,xloc(ind),yloc(ind), 4, log(density1d(ind)), 'filled','s');  %, 'MarkerEdgeColor', 'w')
%     aaa = density2d';
%     bbb = aaa(aaa>=1);
%     ccc = log(bbb);
%     imagesc(ax,xran+dx, yran+dy, ccc);
% create a density matrix to store the number of detections in each small grid
dxhf = 0.2;
dyhf = 0.2;
[density1d,xloc2d,yloc2d,density2d] = density_matrix(hfplt(:,1),hfplt(:,2),...
    xran,yran,dxhf,dyhf);
dumhf = density1d(density1d(:,3)>0, :);
dumhf(dumhf(:,3)>1, :) = [];
scatter(ax,dumhf(:,1),dumhf(:,2), msizehf, log10(dumhf(:,3)),'s','linew',0.2);  %, 'MarkerEdgeColor', 'w')
dumhf = sortrows(density1d(density1d(:,3)>0, :), 3);
dumhf(dumhf(:,3)==1, :) = [];
scatter(ax,dumhf(:,1),dumhf(:,2), msizehf, log10(dumhf(:,3)), 'filled','s');  %, 'MarkerEdgeColor', 'w')
% imagesc(ax,[xran(1)+0.5*dx xran(end)-0.5*dx], [yran(1)+0.5*dy yran(end)-0.5*dy], density2d');
oldcmap = colormap(ax,'jet');
% colormap(ax, flipud(oldcmap) );
c=colorbar(ax,'SouthOutside');
pos = ax.Position;
c.Position = [pos(1), pos(2)-0.018, pos(3), 0.02];
c.Label.String = strcat('log_{10}(N) of detections');
c.Label.FontSize = 11;
text(ax, 0.95, 0.93, '2004+2005','FontSize',12,'unit','normalized','horizontalalignment','right',...
     'EdgeColor','k','Margin',2);
text(ax,0.85,0.8,'HF','FontSize',12,'unit','normalized','horizontalalignment','center');
text(ax,0.04,0.93,'a','FontSize',11,'unit','normalized','EdgeColor','k','Margin',2);
text(ax,0.4,0.93,strcat(num2str(length(hfplt(:,1))),{' detections'}),'FontSize',10,'unit','normalized',...
     'horizontalalignment','center');
ax.Box = 'on';
grid(ax, 'on');
axis(ax, 'equal');
ax.GridLineStyle = '--';
ax.XAxisLocation = 'top';
xlim(ax,xran);
ylim(ax,yran);
xticks(ax,xran(1):5:xran(2));
yticks(ax,yran(1):5:yran(2));
xlabel(ax,'E (km)','fontsize',11);
ylabel(ax,'N (km)','fontsize',11);
hold(ax,'off');

% subplot 2 of figure i
ax = f2.ax(2);
hold(ax,'on');
plot(ax,[-100 100],[0 0],'k--');
plot(ax,[0 0],[-100 100],'k--');
ax.FontSize = 9;
% create a density matrix to store the number of detections in each small grid
dxlf = 0.5;
dylf = 0.5;
[density1d,xloc2d,yloc2d,density2d] = density_matrix(lfplt(:,1),lfplt(:,2),...
    xran,yran,dxlf,dylf);
dumhf = density1d(density1d(:,3)>0, :);
dumhf(dumhf(:,3)>1, :) = [];
scatter(ax,dumhf(:,1),dumhf(:,2), msizelf, log10(dumhf(:,3)),'s','linew',0.2);  %, 'MarkerEdgeColor', 'w')
dumhf = sortrows(density1d(density1d(:,3)>0, :), 3);
dumhf(dumhf(:,3)==1, :) = [];
scatter(ax,dumhf(:,1),dumhf(:,2), msizelf, log10(dumhf(:,3)), 'filled','s');  %, 'MarkerEdgeColor', 'w')
% imagesc(ax,[xran(1)+0.5*dx xran(end)-0.5*dx], [yran(1)+0.5*dy yran(end)-0.5*dy], density2d');
oldcmap = colormap(ax,'jet');
% colormap(ax, flipud(oldcmap) );
c=colorbar(ax,'SouthOutside');
pos = ax.Position;
c.Position = [pos(1), pos(2)-0.018, pos(3), 0.02];
c.Label.String = strcat('log_{10}(N) of detections');
c.Label.FontSize = 11;
text(ax, 0.95, 0.93, '2004+2005','FontSize',12,'unit','normalized','horizontalalignment','right',...
     'EdgeColor','k','Margin',2);
text(ax,0.85,0.8,'LF','FontSize',12,'unit','normalized','horizontalalignment','center');
text(ax,0.04,0.93,'b','FontSize',11,'unit','normalized','EdgeColor','k','Margin',2);
text(ax,0.4,0.93,strcat(num2str(length(lfplt(:,1))),{' detections'}),'FontSize',10,'unit','normalized',...
     'horizontalalignment','center');
ax.Box = 'on';
grid(ax, 'on');
axis(ax, 'equal');
ax.GridLineStyle = '--';
ax.XAxisLocation = 'top';
xlim(ax,xran);
ylim(ax,yran);
xticks(ax,xran(1):5:xran(2));
yticks(ax,yran(1):5:yran(2));
xlabel(ax,'E (km)','fontsize',11);
ylabel(ax,'N (km)','fontsize',11);
hold(ax,'off');


% subplot 3 of figure i
hfplt = mighfa;
lfplt = miglfa;

ax = f2.ax(3);
hold(ax,'on');
plot(ax,[-100 100],[0 0],'k--');
plot(ax,[0 0],[-100 100],'k--');
ax.FontSize = 9;
% create a density matrix to store the number of detections in each small grid
dxhf = 0.2;
dyhf = 0.2;
[density1d,xloc2d,yloc2d,density2d] = density_matrix(hfplt(:,1),hfplt(:,2),...
    xran,yran,dxhf,dyhf);
dumhf = density1d(density1d(:,3)>0, :);
dumhf(dumhf(:,3)>1, :) = [];
scatter(ax,dumhf(:,1),dumhf(:,2), msizehf, log10(dumhf(:,3)),'s','linew',0.2);  %, 'MarkerEdgeColor', 'w')
dumhf = sortrows(density1d(density1d(:,3)>0, :), 3);
dumhf(dumhf(:,3)==1, :) = [];
scatter(ax,dumhf(:,1),dumhf(:,2), msizehf, log10(dumhf(:,3)), 'filled','s');  %, 'MarkerEdgeColor', 'w')
% imagesc(ax,[xran(1)+0.5*dx xran(end)-0.5*dx], [yran(1)+0.5*dy yran(end)-0.5*dy], density2d');
oldcmap = colormap(ax,'jet');
% colormap(ax, flipud(oldcmap) );
c=colorbar(ax,'SouthOutside');
pos = ax.Position;
c.Position = [pos(1), pos(2)-0.018, pos(3), 0.02];
c.Label.String = strcat('log_{10}(N) of detections');
c.Label.FontSize = 11;
text(ax, 0.95, 0.93, 'All RTMs','FontSize',12,'unit','normalized','horizontalalignment','right',...
     'EdgeColor','k','Margin',2);
text(ax,0.85,0.8,'HF','FontSize',12,'unit','normalized','horizontalalignment','center');
text(ax,0.04,0.93,'c','FontSize',11,'unit','normalized','EdgeColor','k','Margin',2);
text(ax,0.4,0.93,strcat(num2str(length(hfplt(:,1))),{' detections'}),'FontSize',10,'unit','normalized',...
     'horizontalalignment','center');
ax.Box = 'on';
grid(ax, 'on');
axis(ax, 'equal');
ax.GridLineStyle = '--';
ax.XAxisLocation = 'top';
xlim(ax,xran);
ylim(ax,yran);
xticks(ax,xran(1):5:xran(2));
yticks(ax,yran(1):5:yran(2));
xlabel(ax,'E (km)','fontsize',11);
ylabel(ax,'N (km)','fontsize',11);
hold(ax,'off');

% subplot 4 of figure i
ax = f2.ax(4);
hold(ax,'on');
plot(ax,[-100 100],[0 0],'k--');
plot(ax,[0 0],[-100 100],'k--');
ax.FontSize = 9;
% create a density matrix to store the number of detections in each small grid
dxlf = 0.5;
dylf = 0.5;
[density1d,xloc2d,yloc2d,density2d] = density_matrix(lfplt(:,1),lfplt(:,2),...
    xran,yran,dxlf,dylf);
dumhf = density1d(density1d(:,3)>0, :);
dumhf(dumhf(:,3)>1, :) = [];
scatter(ax,dumhf(:,1),dumhf(:,2), msizelf, log10(dumhf(:,3)),'s','linew',0.2);  %, 'MarkerEdgeColor', 'w')
dumhf = sortrows(density1d(density1d(:,3)>0, :), 3);
dumhf(dumhf(:,3)==1, :) = [];
scatter(ax,dumhf(:,1),dumhf(:,2), msizelf, log10(dumhf(:,3)), 'filled','s');  %, 'MarkerEdgeColor', 'w')
% imagesc(ax,[xran(1)+0.5*dx xran(end)-0.5*dx], [yran(1)+0.5*dy yran(end)-0.5*dy], density2d');
oldcmap = colormap(ax,'jet');
% colormap(ax, flipud(oldcmap) );
c=colorbar(ax,'SouthOutside');
pos = ax.Position;
c.Position = [pos(1), pos(2)-0.018, pos(3), 0.02];
c.Label.String = strcat('log_{10}(N) of detections');
c.Label.FontSize = 11;
text(ax, 0.95, 0.93, 'All RTMs','FontSize',12,'unit','normalized','horizontalalignment','right',...
     'EdgeColor','k','Margin',2);
text(ax,0.85,0.8,'LF','FontSize',12,'unit','normalized','horizontalalignment','center');
text(ax,0.04,0.93,'d','FontSize',11,'unit','normalized','EdgeColor','k','Margin',2);
text(ax,0.4,0.93,strcat(num2str(length(lfplt(:,1))),{' detections'}),'FontSize',10,'unit','normalized',...
     'horizontalalignment','center');
ax.Box = 'on';
grid(ax, 'on');
axis(ax, 'equal');
ax.GridLineStyle = '--';
ax.XAxisLocation = 'top';
xlim(ax,xran);
ylim(ax,yran);
xticks(ax,xran(1):5:xran(2));
yticks(ax,yran(1):5:yran(2));
xlabel(ax,'E (km)','fontsize',11);
ylabel(ax,'N (km)','fontsize',11);
hold(ax,'off');

print(f2.fig,'-dpdf',strcat(rstpath,'/autortmv3_density_rtm.pdf'));


%%
f3.fig=figure;
f3.fig.Renderer='Painters';
widin = 8;  % maximum width allowed is 8.5 inches
htin = 9;   % maximum height allowed is 11 inches
set(f3.fig,'Position',[scrsz(1)+1*scrsz(3)/10 scrsz(2)+scrsz(4)/20 widin*res htin*res]);
nrow = 1;
ncol = 2;
for isub = 1:nrow*ncol
    f3.ax(isub) = subplot(nrow,ncol,isub);
end

%%% reposition
set(f3.ax(1), 'position', [ 0.08, 0.6, 0.36, 0.32]);
set(f3.ax(2), 'position', [ 0.52, 0.6, 0.36, 0.32]);
% set(f3.ax(3), 'position', [ 0.08, 0.18, 0.36, 0.32]);
% set(f3.ax(4), 'position', [ 0.52, 0.18, 0.36, 0.32]);

% xran = [-14 9];
% yran = [-9 13];

msizehf = 4;
msizelf = 10;

% subplot 1 of figure i
hfplt = hftime(hftime(:,13) > 2004*1000, :);
lfplt = lftime(lftime(:,13) > 2004*1000, :);

ax = f3.ax(1);
hold(ax,'on');
plot(ax,[-100 100],[0 0],'k--');
plot(ax,[0 0],[-100 100],'k--');
ax.FontSize = 9;
%     ind = find(density1d>=0);
%     scatter(ax,xloc(ind),yloc(ind), 4, log(density1d(ind)), 'filled','s');  %, 'MarkerEdgeColor', 'w')
%     aaa = density2d';
%     bbb = aaa(aaa>=1);
%     ccc = log(bbb);
%     imagesc(ax,xran+dx, yran+dy, ccc);
% create a density matrix to store the number of detections in each small grid
dxhf = 0.2;
dyhf = 0.2;
[density1d,xloc2d,yloc2d,density2d] = density_matrix(hfplt(:,1),hfplt(:,2),...
    xran,yran,dxhf,dyhf);
dumhf = density1d(density1d(:,3)>0, :);
dumhf(dumhf(:,3)>1, :) = [];
scatter(ax,dumhf(:,1),dumhf(:,2), msizehf, log10(dumhf(:,3)),'s','linew',0.2);  %, 'MarkerEdgeColor', 'w')
dumhf = sortrows(density1d(density1d(:,3)>0, :), 3);
dumhf(dumhf(:,3)==1, :) = [];
scatter(ax,dumhf(:,1),dumhf(:,2), msizehf, log10(dumhf(:,3)), 'filled','s');  %, 'MarkerEdgeColor', 'w')
% imagesc(ax,[xran(1)+0.5*dx xran(end)-0.5*dx], [yran(1)+0.5*dy yran(end)-0.5*dy], density2d');
oldcmap = colormap(ax,'jet');
% colormap(ax, flipud(oldcmap) );
c=colorbar(ax,'SouthOutside');
pos = ax.Position;
c.Position = [pos(1), pos(2)-0.018, pos(3), 0.02];
c.Label.String = strcat('log_{10}(N) of detections');
c.Label.FontSize = 11;
text(ax, 0.95, 0.93, '2004+2005','FontSize',12,'unit','normalized','horizontalalignment','right',...
     'EdgeColor','k','Margin',2);
text(ax,0.85,0.8,'HF','FontSize',12,'unit','normalized','horizontalalignment','center');
text(ax,0.04,0.93,'a','FontSize',11,'unit','normalized','EdgeColor','k','Margin',2);
text(ax,0.4,0.93,strcat(num2str(length(hfplt(:,1))),{' detections'}),'FontSize',10,'unit','normalized',...
     'horizontalalignment','center');
ax.Box = 'on';
grid(ax, 'on');
axis(ax, 'equal');
ax.GridLineStyle = '--';
ax.XAxisLocation = 'top';
xlim(ax,xran);
ylim(ax,yran);
xticks(ax,xran(1):5:xran(2));
yticks(ax,yran(1):5:yran(2));
xlabel(ax,'E (km)','fontsize',11);
ylabel(ax,'N (km)','fontsize',11);
hold(ax,'off');

% subplot 2 of figure i
ax = f3.ax(2);
hold(ax,'on');
plot(ax,[-100 100],[0 0],'k--');
plot(ax,[0 0],[-100 100],'k--');
ax.FontSize = 9;
% create a density matrix to store the number of detections in each small grid
dxlf = 0.5;
dylf = 0.5;
[density1d,xloc2d,yloc2d,density2d] = density_matrix(lfplt(:,1),lfplt(:,2),...
    xran,yran,dxlf,dylf);
dumhf = density1d(density1d(:,3)>0, :);
dumhf(dumhf(:,3)>1, :) = [];
scatter(ax,dumhf(:,1),dumhf(:,2), msizelf, log10(dumhf(:,3)),'s','linew',0.2);  %, 'MarkerEdgeColor', 'w')
dumhf = sortrows(density1d(density1d(:,3)>0, :), 3);
dumhf(dumhf(:,3)==1, :) = [];
scatter(ax,dumhf(:,1),dumhf(:,2), msizelf, log10(dumhf(:,3)), 'filled','s');  %, 'MarkerEdgeColor', 'w')
% imagesc(ax,[xran(1)+0.5*dx xran(end)-0.5*dx], [yran(1)+0.5*dy yran(end)-0.5*dy], density2d');
oldcmap = colormap(ax,'jet');
% colormap(ax, flipud(oldcmap) );
c=colorbar(ax,'SouthOutside');
pos = ax.Position;
c.Position = [pos(1), pos(2)-0.018, pos(3), 0.02];
c.Label.String = strcat('log_{10}(N) of detections');
c.Label.FontSize = 11;
text(ax, 0.95, 0.93, '2004+2005','FontSize',12,'unit','normalized','horizontalalignment','right',...
     'EdgeColor','k','Margin',2);
text(ax,0.85,0.8,'LF','FontSize',12,'unit','normalized','horizontalalignment','center');
text(ax,0.04,0.93,'b','FontSize',11,'unit','normalized','EdgeColor','k','Margin',2);
text(ax,0.4,0.93,strcat(num2str(length(lfplt(:,1))),{' detections'}),'FontSize',10,'unit','normalized',...
     'horizontalalignment','center');
ax.Box = 'on';
grid(ax, 'on');
axis(ax, 'equal');
ax.GridLineStyle = '--';
ax.XAxisLocation = 'top';
xlim(ax,xran);
ylim(ax,yran);
xticks(ax,xran(1):5:xran(2));
yticks(ax,yran(1):5:yran(2));
xlabel(ax,'E (km)','fontsize',11);
ylabel(ax,'N (km)','fontsize',11);
hold(ax,'off');

print(f3.fig,'-dpdf',strcat(rstpath,'/autortmv3_density_ETS2004+2005.pdf'));


%%

mighfne = [];
miglfne = [];
mighfse = [];
miglfse = [];
mighfsw = [];
miglfsw = [];
mighfnw = [];
miglfnw = [];

% indne = [];
% indse = [];
% indsw = [];
% indnw = [];

% used migrations with option 4 in 'mig_linear_fit_LZB_addfam_autortm_v3'
indne = [10,18,21,30,34,46,47];
indsw = [5,6,7,9,12,13,15,22,23,28,33,35,37,38,41,43,44,49,51,52];

indlargediff = [2,3,11,14,24,29,32,34,36,40,43,44,49,52];  %
indpoorlf = [40,4,5,33,3,15,11,23,35,21,45,17,2,31,30,24,32,39,6,9,27,37,51,19,48,26,46,22,8,29,18,47]; %

% indne = setdiff(indne, largediffind); 
% indsw = setdiff(indsw, largediffind);
% indne = setdiff(indne, poorlfind); 
% indsw = setdiff(indsw, poorlfind);

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

%     disp(trange(i,:));
    indhf = find(hftime(:,13)==trange(i,1) & hftime(:,15)>=trange(i,2) & ...
                 hftime(:,15)<=trange(i,3));
    mighf = hftime(indhf,:);
    
    indlf = find(lftime(:,13)==trange(i,1) & lftime(:,15)>=trange(i,2) & ...
                 lftime(:,15)<=trange(i,3));
    miglf = lftime(indlf,:); 
    
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
    
    % rotate back the projected locations to get the original coordinates
    tmphf = mighfdum2;
    [newx,newy] = coordinate_rot(tmphf(:,1),tmphf(:,2),(angbest(i)-90),0,0);
    tmphf(:,1) = newx;
    tmphf(:,2) = newy;
    tmplf = miglfdum2;
    [newx,newy] = coordinate_rot(tmplf(:,1),tmplf(:,2),(angbest(i)-90),0,0);
    tmplf(:,1) = newx;
    tmplf(:,2) = newy;
        
%     if angbest(i)>0 && angbest(i)<=90
%         indne = [indne; i];
%         mighfne = [mighfne; tmphf];
%         miglfne = [miglfne; tmplf];    
%     elseif angbest(i)>90 && angbest(i)<=180
%         indse = [indse; i];
%         mighfse = [mighfse; tmphf];
%         miglfse = [miglfse; tmplf];
%     elseif angbest(i)>180 && angbest(i)<=270
%         indsw = [indsw; i];
%         mighfsw = [mighfsw; tmphf];
%         miglfsw = [miglfsw; tmplf];
%     else
%         indnw = [indnw; i];
%         mighfnw = [mighfnw; tmphf];
%         miglfnw = [miglfnw; tmplf];
%     end

    if ~isempty(intersect(indne,i))
        mighfne = [mighfne; tmphf];
        miglfne = [miglfne; tmplf];    
    elseif ~isempty(intersect(indsw,i))
        mighfsw = [mighfsw; tmphf];
        miglfsw = [miglfsw; tmplf];
    end
    
end

%% RTM density, NE vs. SW
f.fig=figure;
f.fig.Renderer='Painters';
widin = 8;  % maximum width allowed is 8.5 inches
htin = 9;   % maximum height allowed is 11 inches
set(f.fig,'Position',[scrsz(1)+1*scrsz(3)/10 scrsz(2)+scrsz(4)/20 widin*res htin*res]);
nrow = 2;
ncol = 2;
for isub = 1:nrow*ncol
    f.ax(isub) = subplot(nrow,ncol,isub);
end

%%% reposition
% set(f.ax(1), 'position', [ 0.1, 0.55, 0.4, 0.4]);
% set(f.ax(2), 'position', [ 0.55, 0.55, 0.4, 0.4]);
% set(f.ax(3), 'position', [ 0.1, 0.1, 0.4, 0.4]);
% set(f.ax(4), 'position', [ 0.55, 0.1, 0.4, 0.4]);
set(f.ax(1), 'position', [ 0.08, 0.6, 0.36, 0.32]);
set(f.ax(2), 'position', [ 0.52, 0.6, 0.36, 0.32]);
set(f.ax(3), 'position', [ 0.08, 0.18, 0.36, 0.32]);
set(f.ax(4), 'position', [ 0.52, 0.18, 0.36, 0.32]);

msizehf = 4;
msizelf = 10;

% subplot 1 of figure i
hfplt = mighfne;
lfplt = miglfne;

ax = f.ax(1);
hold(ax,'on');
plot(ax,[-100 100],[0 0],'k--');
plot(ax,[0 0],[-100 100],'k--');
ax.FontSize = 9;
%     ind = find(density1d>=0);
%     scatter(ax,xloc(ind),yloc(ind), 4, log(density1d(ind)), 'filled','s');  %, 'MarkerEdgeColor', 'w')
%     aaa = density2d';
%     bbb = aaa(aaa>=1);
%     ccc = log(bbb);
%     imagesc(ax,xran+dx, yran+dy, ccc);
% create a density matrix to store the number of detections in each small grid
dxhf = 0.2;
dyhf = 0.2;
[density1d,xloc2d,yloc2d,density2d] = density_matrix(hfplt(:,1),hfplt(:,2),...
    xran,yran,dxhf,dyhf);
dumhf = density1d(density1d(:,3)>0, :);
dumhf(dumhf(:,3)>1, :) = [];
scatter(ax,dumhf(:,1),dumhf(:,2), msizehf, log10(dumhf(:,3)),'s','linew',0.2);  %, 'MarkerEdgeColor', 'w')
dumhf = sortrows(density1d(density1d(:,3)>0, :), 3);
dumhf(dumhf(:,3)==1, :) = [];
scatter(ax,dumhf(:,1),dumhf(:,2), msizehf, log10(dumhf(:,3)), 'filled','s');  %, 'MarkerEdgeColor', 'w')
% imagesc(ax,[xran(1)+0.5*dx xran(end)-0.5*dx], [yran(1)+0.5*dy yran(end)-0.5*dy], density2d');
oldcmap = colormap(ax,'jet');
% colormap(ax, flipud(oldcmap) );
c=colorbar(ax,'SouthOutside');
pos = ax.Position;
c.Position = [pos(1), pos(2)-0.018, pos(3), 0.02];
c.Label.String = strcat('log_{10}(N) of detections');
c.Label.FontSize = 11;
text(ax, 0.95, 0.93, 'ENE-propagating RTMs','FontSize',12,'unit','normalized','horizontalalignment',...
     'right','EdgeColor','k','Margin',2);
text(ax,0.04,0.93,'a','FontSize',11,'unit','normalized','EdgeColor','k','Margin',2);
text(ax,0.8,0.1,'HF','FontSize',12,'unit','normalized','horizontalalignment','center');
text(ax,0.99,0.8,strcat(num2str(size(hfplt,1)),{' detections'}),'FontSize',10,'unit','normalized',...
     'horizontalalignment','right');
ax.Box = 'on';
grid(ax, 'on');
axis(ax, 'equal');
ax.GridLineStyle = '--';
ax.XAxisLocation = 'top';
xlim(ax,xran);
ylim(ax,yran);
xticks(ax,xran(1):5:xran(2));
yticks(ax,yran(1):5:yran(2));
xlabel(ax,'E (km)','fontsize',11);
ylabel(ax,'N (km)','fontsize',11);
hold(ax,'off');

% subplot 2 of figure i
ax = f.ax(2);
hold(ax,'on');
plot(ax,[-100 100],[0 0],'k--');
plot(ax,[0 0],[-100 100],'k--');
ax.FontSize = 9;
% create a density matrix to store the number of detections in each small grid
dxlf = 0.5;
dylf = 0.5;
[density1d,xloc2d,yloc2d,density2d] = density_matrix(lfplt(:,1),lfplt(:,2),...
    xran,yran,dxlf,dylf);
dumlf = density1d(density1d(:,3)>0, :);
dumlf(dumlf(:,3)>1, :) = [];
scatter(ax,dumlf(:,1),dumlf(:,2), msizelf, log10(dumlf(:,3)),'s','linew',0.2);  %, 'MarkerEdgeColor', 'w')
dumlf = sortrows(density1d(density1d(:,3)>0, :), 3);
dumlf(dumlf(:,3)==1, :) = [];
scatter(ax,dumlf(:,1),dumlf(:,2), msizelf, log10(dumlf(:,3)), 'filled','s');  %, 'MarkerEdgeColor', 'w')
% imagesc(ax,[xran(1)+0.5*dx xran(end)-0.5*dx], [yran(1)+0.5*dy yran(end)-0.5*dy], density2d');
oldcmap = colormap(ax,'jet');
% colormap(ax, flipud(oldcmap) );
c=colorbar(ax,'SouthOutside');
pos = ax.Position;
c.Position = [pos(1), pos(2)-0.018, pos(3), 0.02];
c.Label.String = strcat('log_{10}(N) of detections');
c.Label.FontSize = 11;
text(ax, 0.95, 0.93, 'ENE-propagating RTMs','FontSize',12,'unit','normalized','horizontalalignment',...
     'right','EdgeColor','k','Margin',2);
text(ax,0.04,0.93,'b','FontSize',11,'unit','normalized','EdgeColor','k','Margin',2);
text(ax,0.8,0.1,'LF','FontSize',12,'unit','normalized','horizontalalignment','center');
text(ax,0.99,0.8,strcat(num2str(size(lfplt,1)),{' detections'}),'FontSize',10,'unit','normalized',...
     'horizontalalignment','right');
ax.Box = 'on';
grid(ax, 'on');
axis(ax, 'equal');
ax.GridLineStyle = '--';
ax.XAxisLocation = 'top';
xlim(ax,xran);
ylim(ax,yran);
xticks(ax,xran(1):5:xran(2));
yticks(ax,yran(1):5:yran(2));
xlabel(ax,'E (km)','fontsize',11);
ylabel(ax,'N (km)','fontsize',11);
hold(ax,'off');

% subplot 3 of figure i
hfplt = mighfsw;
lfplt = miglfsw;

ax = f.ax(3);
hold(ax,'on');
plot(ax,[-100 100],[0 0],'k--');
plot(ax,[0 0],[-100 100],'k--');
ax.FontSize = 9;
%     ind = find(density1d>=0);
%     scatter(ax,xloc(ind),yloc(ind), 4, log(density1d(ind)), 'filled','s');  %, 'MarkerEdgeColor', 'w')
%     aaa = density2d';
%     bbb = aaa(aaa>=1);
%     ccc = log(bbb);
%     imagesc(ax,xran+dx, yran+dy, ccc);
% create a density matrix to store the number of detections in each small grid
dxhf = 0.2;
dyhf = 0.2;
[density1d,xloc2d,yloc2d,density2d] = density_matrix(hfplt(:,1),hfplt(:,2),...
    xran,yran,dxhf,dyhf);
dumhf = density1d(density1d(:,3)>0, :);
dumhf(dumhf(:,3)>1, :) = [];
scatter(ax,dumhf(:,1),dumhf(:,2), msizehf, log10(dumhf(:,3)),'s','linew',0.2);  %, 'MarkerEdgeColor', 'w')
dumhf = sortrows(density1d(density1d(:,3)>0, :), 3);
dumhf(dumhf(:,3)==1, :) = [];
scatter(ax,dumhf(:,1),dumhf(:,2), msizehf, log10(dumhf(:,3)), 'filled','s');  %, 'MarkerEdgeColor', 'w')
% imagesc(ax,[xran(1)+0.5*dx xran(end)-0.5*dx], [yran(1)+0.5*dy yran(end)-0.5*dy], density2d');
oldcmap = colormap(ax,'jet');
% colormap(ax, flipud(oldcmap) );
c=colorbar(ax,'SouthOutside');
pos = ax.Position;
c.Position = [pos(1), pos(2)-0.018, pos(3), 0.02];
c.Label.String = strcat('log_{10}(N) of detections');
c.Label.FontSize = 11;
text(ax, 0.95, 0.93, 'WSW-propagating RTMs','FontSize',12,'unit','normalized','horizontalalignment',...
     'right','EdgeColor','k','Margin',2);
text(ax,0.04,0.93,'c','FontSize',11,'unit','normalized','EdgeColor','k','Margin',2);
text(ax,0.8,0.1,'HF','FontSize',12,'unit','normalized','horizontalalignment','center');
text(ax,0.99,0.8,strcat(num2str(size(hfplt,1)),{' detections'}),'FontSize',10,'unit','normalized',...
     'horizontalalignment','right');
ax.Box = 'on';
grid(ax, 'on');
axis(ax, 'equal');
ax.GridLineStyle = '--';
ax.XAxisLocation = 'top';
xlim(ax,xran);
ylim(ax,yran);
xticks(ax,xran(1):5:xran(2));
yticks(ax,yran(1):5:yran(2));
xlabel(ax,'E (km)','fontsize',11);
ylabel(ax,'N (km)','fontsize',11);
hold(ax,'off');

% subplot 4 of figure i
ax = f.ax(4);
hold(ax,'on');
plot(ax,[-100 100],[0 0],'k--');
plot(ax,[0 0],[-100 100],'k--');
ax.FontSize = 9;
% create a density matrix to store the number of detections in each small grid
dxlf = 0.5;
dylf = 0.5;
[density1d,xloc2d,yloc2d,density2d] = density_matrix(lfplt(:,1),lfplt(:,2),...
    xran,yran,dxlf,dylf);
dumlf = density1d(density1d(:,3)>0, :);
dumlf(dumlf(:,3)>1, :) = [];
scatter(ax,dumlf(:,1),dumlf(:,2), msizelf, log10(dumlf(:,3)),'s','linew',0.2);  %, 'MarkerEdgeColor', 'w')
dumlf = sortrows(density1d(density1d(:,3)>0, :), 3);
dumlf(dumlf(:,3)==1, :) = [];
scatter(ax,dumlf(:,1),dumlf(:,2), msizelf, log10(dumlf(:,3)), 'filled','s');  %, 'MarkerEdgeColor', 'w')
% imagesc(ax,[xran(1)+0.5*dx xran(end)-0.5*dx], [yran(1)+0.5*dy yran(end)-0.5*dy], density2d');
oldcmap = colormap(ax,'jet');
% colormap(ax, flipud(oldcmap) );
c=colorbar(ax,'SouthOutside');
pos = ax.Position;
c.Position = [pos(1), pos(2)-0.018, pos(3), 0.02];
c.Label.String = strcat('log_{10}(N) of detections');
c.Label.FontSize = 11;
text(ax, 0.95, 0.93, 'WSW-propagating RTMs','FontSize',12,'unit','normalized','horizontalalignment',...
     'right','EdgeColor','k','Margin',2);
text(ax,0.04,0.93,'d','FontSize',11,'unit','normalized','EdgeColor','k','Margin',2);
text(ax,0.8,0.1,'LF','FontSize',12,'unit','normalized','horizontalalignment','center');
text(ax,0.99,0.8,strcat(num2str(size(lfplt,1)),{' detections'}),'FontSize',10,'unit','normalized',...
     'horizontalalignment','right');
ax.Box = 'on';
grid(ax, 'on');
axis(ax, 'equal');
ax.GridLineStyle = '--';
ax.XAxisLocation = 'top';
xlim(ax,xran);
ylim(ax,yran);
xticks(ax,xran(1):5:xran(2));
yticks(ax,yran(1):5:yran(2));
xlabel(ax,'E (km)','fontsize',11);
ylabel(ax,'N (km)','fontsize',11);
hold(ax,'off');

print(f.fig,'-dpdf',strcat(rstpath,'/autortmv3_density_rtm_NE+SW.pdf'));















