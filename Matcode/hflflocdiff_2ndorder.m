% function mig_linear_fit_LZB_addfam_autortm_v3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% We thought the migration #13 in which LF starts earlier than HF is caused
% by the bad correlation between HF sources at earlier time.
% Then we tried to lower the CC threhold 0.4-->0.1-->0.01, even increase the
% loopoffmax 1.5-->2.5, however, it seems that HF signals are not increased
% noticably in that time window.
% Assuming our detecting range is enough and the missing is not due to the  
% distance cutoff, then it is possible that that region produces LF energy
% only.
% To test it, the first thing to do is to see in that polygon, when did the 
% sparse HF and more numerous LF occur? Are they random in time? How many 
% are within the migration, under the same distance cutoff (eg, 12 km)
%
%
%
% Chao Song, chaosong@princeton.edu
% First created date:   2021/02/15
% Last modified date:   2021/02/15
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Initialization
format short e   % Set the format to 5-digit floating point
clear
close all
clc

set(0,'DefaultFigureVisible','on');
% set(0,'DefaultFigureVisible','off');   % switch to show the plots or not

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
fname1 = strcat(rstpath, '/evtloc.all',num2str(nfam),'fam.nodcutnodou.',SUFFIXhf);
fname2 = strcat(rstpath, '/evtloc.all',num2str(nfam),'fam.eq8kmdcutnodou.',SUFFIXhf);
fname3 = strcat(rstpath, '/evtloc.all',num2str(nfam),'fam.eq12kmdcutnodou.',SUFFIXhf);

fname = fname3;
hftime = load(fname);
% 25 cols, format is:
%   E(043) N(043) E(own) N(own) dist(own) lon lat dep off12 off13 off12sec off13sec date
%   main_arrival_time  cent_of_win  avecc ccadd1 ccadd2 ccadd3 ccadd4 ampmax ampmin eng normeng famnum
%


SUFFIXlf = strcat('lf.time.',num2str(winlenhf),'_',num2str(winlenlf),'.',num2str(lofflf),...
                  '.',num2str(ccminlf));
fname1 = strcat(rstpath, '/evtloc.all',num2str(nfam),'fam.nodcutnodou.',SUFFIXlf);
fname2 = strcat(rstpath, '/evtloc.all',num2str(nfam),'fam.eq8kmdcutnodou.',SUFFIXlf);
fname3 = strcat(rstpath, '/evtloc.all',num2str(nfam),'fam.eq12kmdcutnodou.',SUFFIXlf);

fname = fname3;
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

hftime = hftime(hftime(:,13) > 2004*1000, :);
lftime = lftime(lftime(:,13) > 2004*1000, :);

%% 
% this is finally accepted new time ranges that come from 'identify_RTMs_v3.m'
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

trange1(:,1) = trange(:,1);
trange1(:,2) = trange(:,2)-0.2*3600;
trange1(:,3) = trange(:,3)+0.2*3600;

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

% a box around fam 002 region
% reg1 = [ 11 -3;
%            20 -1.5;
%            20 5;
%            11 5;
%            11 -3];
       
reg1 = [ 10 0;
           20 0;
           20 5;
           10 5;
           10 0];
       
%%
f2.fig=figure;
f2.fig.Renderer='Painters';
widin = 8;  % maximum width allowed is 8.5 inches
htin = 9;   % maximum height allowed is 11 inches
set(f2.fig,'Position',[scrsz(1)+1*scrsz(3)/10 scrsz(2)+scrsz(4)/20 widin*res htin*res]);
nrow = 1;
ncol = 2;
for isub = 1:nrow*ncol
    f2.ax(isub) = subplot(nrow,ncol,isub);
end

%%% reposition
set(f2.ax(1), 'position', [ 0.08, 0.6, 0.36, 0.32]);
set(f2.ax(2), 'position', [ 0.52, 0.6, 0.36, 0.32]);
% set(f2.ax(3), 'position', [ 0.08, 0.18, 0.36, 0.32]);
% set(f2.ax(4), 'position', [ 0.52, 0.18, 0.36, 0.32]);

xran = [-25 30];
yran = [-20 30];

msizehf = 4;
msizelf = 10;

% subplot 1 of figure i
% hfplt = hftime(hftime(:,13) > 2004*1000, :);
% lfplt = lftime(lftime(:,13) > 2004*1000, :);
hfplt = hftime;
lfplt = lftime;

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
plot(ax,reg1(:,1),reg1(:,2),'k-','linew',2);
% imagesc(ax,[xran(1)+0.5*dx xran(end)-0.5*dx], [yran(1)+0.5*dy yran(end)-0.5*dy], density2d');
oldcmap = colormap(ax,'jet');
% colormap(ax, flipud(oldcmap) );
c=colorbar(ax,'SouthOutside');
pos = ax.Position;
c.Position = [pos(1), pos(2)-0.018, pos(3), 0.02];
c.Label.String = strcat('log_{10}(N) of detections');
c.Label.FontSize = 11;
text(ax, 0.95, 0.93, 'Entire catalog','FontSize',12,'unit','normalized','horizontalalignment','right',...
     'EdgeColor','k','Margin',2);
text(ax,0.85,0.8,'HF','FontSize',12,'unit','normalized','horizontalalignment','center');
text(ax,0.04,0.93,'a','FontSize',11,'unit','normalized','EdgeColor','k','Margin',2);
text(ax,0.4,0.93,strcat(num2str(length(hfplt)),{' detections'}),'FontSize',10,'unit','normalized',...
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
plot(ax,reg1(:,1),reg1(:,2),'k-','linew',2);
% imagesc(ax,[xran(1)+0.5*dx xran(end)-0.5*dx], [yran(1)+0.5*dy yran(end)-0.5*dy], density2d');
oldcmap = colormap(ax,'jet');
% colormap(ax, flipud(oldcmap) );
c=colorbar(ax,'SouthOutside');
pos = ax.Position;
c.Position = [pos(1), pos(2)-0.018, pos(3), 0.02];
c.Label.String = strcat('log_{10}(N) of detections');
c.Label.FontSize = 11;
text(ax, 0.95, 0.93, 'Entire catalog','FontSize',12,'unit','normalized','horizontalalignment','right',...
     'EdgeColor','k','Margin',2);
text(ax,0.85,0.8,'LF','FontSize',12,'unit','normalized','horizontalalignment','center');
text(ax,0.04,0.93,'b','FontSize',11,'unit','normalized','EdgeColor','k','Margin',2);
text(ax,0.4,0.93,strcat(num2str(length(lfplt)),{' detections'}),'FontSize',10,'unit','normalized',...
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


%%
% find detections inside the bound
[is,ion] = inpolygon(hftime(:,1),hftime(:,2),reg1(:,1),reg1(:,2));
isinreg = is | ion;
hfreg1 = hftime(isinreg == 1, :);

[is,ion] = inpolygon(lftime(:,1),lftime(:,2),reg1(:,1),reg1(:,2));
isinreg = is | ion;
lfreg1 = lftime(isinreg == 1, :);


% how many detections are inside the migration?        
hfreg1mig = cell(size(trange1, 1), 1);
lfreg1mig = cell(size(trange1, 1), 1);
hfreg1in = [];
lfreg1in = [];
for i = 1: size(trange1, 1)
%     disp(trange1(i,:));
    indhf = find(hfreg1(:,13)==trange1(i,1) & hfreg1(:,15)>=trange1(i,2) & ...
                 hfreg1(:,15)<=trange1(i,3));
    mighf = hfreg1(indhf,:);
    hfreg1in = [hfreg1in; mighf];
    hfreg1mig{i} = mighf;
    
    indlf = find(lfreg1(:,13)==trange1(i,1) & lfreg1(:,15)>=trange1(i,2) & ...
                 lfreg1(:,15)<=trange1(i,3));
    miglf = lfreg1(indlf,:); 
    lfreg1in = [lfreg1in; miglf];
    lfreg1mig{i} = miglf;
    
end

% it seems that the rate of detections inside the bound relative to catalog is comparable
% but the rate of detections in migrations relative to that inside the bound for LF is clearly
% higher
perchfreg1 = size(hfreg1,1)/size(hftime,1)*100;
perclfreg1 = size(lfreg1,1)/size(lftime,1)*100;
ratreg1 = perclfreg1/perchfreg1

perchfin = size(hfreg1in,1)/size(hfreg1,1)*100;
perclfin = size(lfreg1in,1)/size(lfreg1,1)*100;
ratin = perclfin/perchfin


% next we check which migrations are of interest
indmighf = [];
indmiglf = [];
for i = 1: size(trange1, 1)
    if ~isempty(hfreg1mig{i})
        indmighf = [indmighf i];
    end
    if ~isempty(lfreg1mig{i})
        indmiglf = [indmiglf i];
    end
end

indmig = union(indmighf,indmiglf)'

indmig = [13,38,49];

%% plot the potentially useful migrations
xran = [-20 25];
yran = [-20 20];

for i = 1: length(indmig)
% for i = 1: 211
%     i=46;
    indhf = find(hftime(:,13)==trange1(indmig(i),1) & hftime(:,15)>=trange1(indmig(i),2) & ...
                 hftime(:,15)<=trange1(indmig(i),3));
    mighf = hftime(indhf,:);
    
    indlf = find(lftime(:,13)==trange1(indmig(i),1) & lftime(:,15)>=trange1(indmig(i),2) & ...
                 lftime(:,15)<=trange1(indmig(i),3));
    miglf = lftime(indlf,:); 
    
    %%% Actually used is the pre-determined best prop direc to do the fitting
    mighfdum = mighf;
    for j = 1: size(mighf,1)
        x0 = mighfdum(j,1);
        y0 = mighfdum(j,2);
        [newx,newy] = coordinate_rot(x0,y0,-(angbest(indmig(i))-90),0,0);
        mighfdum(j,1) = newx;
        mighfdum(j,2) = newy;
    end

    miglfdum = miglf;
    for j = 1: size(miglf,1)
        x0 = miglfdum(j,1);
        y0 = miglfdum(j,2);
        [newx,newy] = coordinate_rot(x0,y0,-(angbest(indmig(i))-90),0,0);
        miglfdum(j,1) = newx;
        miglfdum(j,2) = newy;
    end

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
    
    % subplot 1 of figure indmig(i)
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
    juldate = num2str(trange1(indmig(i),1));
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
%     c.Label.String = strcat(num2str(trange1(indmig(i),1)),' of HF',' (hr)');
    c.Label.FontSize = 11;
    caxis(f.ax(1),[trange1(indmig(i),2)/3600 trange1(indmig(i),3)/3600])
    text(f.ax(1),0.85,0.1,'HF','FontSize',12,'unit','normalized');
    text(f.ax(1),0.04,0.93,'a','FontSize',11,'unit','normalized','EdgeColor','k','Margin',2);
    text(f.ax(1),0.04,0.1,'LZB','FontSize',11,'unit','normalized','EdgeColor','k','Margin',2);
    f.ax(1).Box = 'on';
    grid(f.ax(1), 'on');
    axis(f.ax(1), 'equal');
    f.ax(1).GridLineStyle = '--';
    f.ax(1).XAxisLocation = 'top';
    medxhf = median(mighf(:,1));
    medyhf = median(mighf(:,2));
    [rotx, roty] = complex_rot(0,5,-angbest(indmig(i)));
    xvect = [medxhf-rotx medxhf+rotx];
    yvect = [medyhf-roty medyhf+roty];
    drawArrow(f.ax(1),xvect,yvect,xran,yran,'linewidth',1.5);
    [rotx, roty] = complex_rot(-5,0,-angbest(indmig(i)));
    xvect = [medxhf-rotx medxhf+rotx];
    yvect = [medyhf-roty medyhf+roty];
    drawArrow(f.ax(1),xvect,yvect,xran,yran,'linewidth',1.5,'linestyle','--','color',[0.4 0.4 0.4]);
    xticks(f.ax(1),xran(1):5:xran(2));
    yticks(f.ax(1),yran(1):5:yran(2));    
    xlabel(f.ax(1),'E (km)','fontsize',11);
    ylabel(f.ax(1),'N (km)','fontsize',11);
    text(f.ax(1),0.5,0.1,num2str(indmig(i)),'FontSize',12,'unit','normalized');
    text(f.ax(1),0.84,0.95,strcat(num2str(size(mighf,1)),{' detections'}),'FontSize',8,...
         'unit','normalized','horizontalalignment','center');
    text(f.ax(1),0.84,0.90,strcat({'in '},num2str(trange1(indmig(i),3)-trange1(indmig(i),2)),{' s'}),'FontSize',8,...
         'unit','normalized','horizontalalignment','center');
    rate = sprintf('%.3f',size(mighf,1)/(trange1(indmig(i),3)-trange1(indmig(i),2)));
    text(f.ax(1),0.84,0.85,strcat({'rate: '},rate),'FontSize',8,...
         'unit','normalized','horizontalalignment','center');
    text(f.ax(1),0.84,0.80,strcat(num2str(angbest(indmig(i))),{'{\circ}'}),'FontSize',10,...
         'unit','normalized','horizontalalignment','center'); 
    hold(f.ax(1),'off');

    % subplot 2 of figure indmig(i)
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
%     c.Label.String = strcat(num2str(trange1(indmig(i),1)),' of LF',' (hr)');
    c.Label.FontSize = 11;
    caxis(f.ax(2),[trange1(indmig(i),2)/3600 trange1(indmig(i),3)/3600])
    text(f.ax(2),0.85,0.1,'LF','FontSize',12,'unit','normalized');
    text(f.ax(2),0.04,0.93,'b','FontSize',11,'unit','normalized','EdgeColor','k','Margin',2);
    text(f.ax(2),0.04,0.1,'LZB','FontSize',11,'unit','normalized','EdgeColor','k','Margin',2);
    f.ax(2).Box = 'on';
    grid(f.ax(2), 'on');
    axis(f.ax(2), 'equal');
    f.ax(2).GridLineStyle = '--';
    f.ax(2).XAxisLocation = 'top';
    medxlf = median(miglf(:,1));
    medylf = median(miglf(:,2));
    [rotx, roty] = complex_rot(0,5,-angbest(indmig(i)));
    xvect = [medxlf-rotx medxlf+rotx];
    yvect = [medylf-roty medylf+roty];
    drawArrow(f.ax(2),xvect,yvect,xran,yran,'linewidth',1.5);
    [rotx, roty] = complex_rot(-5,0,-angbest(indmig(i)));
    xvect = [medxlf-rotx medxlf+rotx];
    yvect = [medylf-roty medylf+roty];
    drawArrow(f.ax(2),xvect,yvect,xran,yran,'linewidth',1.5,'linestyle','--','color',[0.4 0.4 0.4]);
    xticks(f.ax(2),xran(1):5:xran(2));
    yticks(f.ax(2),yran(1):5:yran(2));    
    xlabel(f.ax(2),'E (km)','fontsize',11);
    ylabel(f.ax(2),'N (km)','fontsize',11);
    text(f.ax(2),0.5,0.1,num2str(indmig(i)),'FontSize',12,'unit','normalized');
    text(f.ax(2),0.84,0.95,strcat(num2str(size(miglf,1)),{' detections'}),'FontSize',8,...
         'unit','normalized','horizontalalignment','center');
    text(f.ax(2),0.84,0.9,strcat({'in '},num2str(trange1(indmig(i),3)-trange1(indmig(i),2)),{' s'}),'FontSize',8,...
         'unit','normalized','horizontalalignment','center');
    rate = sprintf('%.3f',size(miglf,1)/(trange1(indmig(i),3)-trange1(indmig(i),2)));
    text(f.ax(2),0.84,0.85,strcat({'rate: '},rate),'FontSize',8,...
         'unit','normalized','horizontalalignment','center');
    text(f.ax(2),0.84,0.80,strcat(num2str(angbest(indmig(i))),{'{\circ}'}),'FontSize',10,...
         'unit','normalized','horizontalalignment','center');
    hold(f.ax(2),'off');

    % subplot 3 of figure i
    hold(f.ax(3),'on');
    f.ax(3).FontSize = 9;
    scatter(f.ax(3),mighfdum(:,15)/3600,mighfdum(:,1),20,[0.6 0.6 0.6],'filled','o',...
            'MarkerEdgeColor','k');    
    scatter(f.ax(3),miglfdum(:,15)/3600,miglfdum(:,1),20,[0.6 1 1],'filled','o',...
            'MarkerEdgeColor','k');
    f.ax(3).Box = 'on';
    grid(f.ax(3), 'on');
    f.ax(3).GridLineStyle = '--';
    xran1 = [trange1(indmig(i),2)/3600 trange1(indmig(i),3)/3600];
%     yran1 = [round(min([mighfdum(:,1);miglfdum(:,1)]))-1 ...
%              round(max([mighfdum(:,1);miglfdum(:,1)]))+1];
    aa = round(prctile([mighfdum(:,1);miglfdum(:,1)], 98));
    bb = round(prctile([mighfdum(:,1);miglfdum(:,1)], 2));
    aap = aa + ceil((aa-bb)/6);
    bbp = bb - ceil((aa-bb)/3);
    yran1 = [bbp aap];
    xlim(f.ax(3),xran1);
    ylim(f.ax(3),yran1);
    plot(f.ax(3),[trange(indmig(i),2)/3600  trange(indmig(i),2)/3600], yran1, 'b--','linew',2);
    plot(f.ax(3),[trange(indmig(i),3)/3600  trange(indmig(i),3)/3600], yran1, 'b--','linew',2);
    
    text(f.ax(3),0.04,0.91,'c','FontSize',11,'unit','normalized','EdgeColor','k','Margin',2);
    xlabel(f.ax(3),'Time (hr)','fontsize',11);
    ylabel(f.ax(3),'Dist. along prop. (km)','fontsize',11); 
    hold(f.ax(3),'off');
    
    
end

%%
% note the indice of migrations that count
induseful = [49,38,13];

% ** #49 interesting case, though no HF, later LF themselves are coherent enough, seem real
% #42 unfortunately only 1 LF in region, likely to be false detections
% #39 LF likely to be false detections
% ** #38 highlighted one, earlier LF detections seem real and tied to the same mig, coherent
% #32 LF likely to be false detections
% #31 later LF likely to be real with larger errors 
% #30 LF likely to be false detections
% #28 LF likely to be false detections
% #27 is 002 migration, several HF (047, 002) and LF (002) likely deviated from cluster
% #25 LF likely to be false detections
% #18 hard to tell if earlier LF detections are not real 
% #17 later LF detections not tied to mig, likely have large errors
% ** #13 highlighted mig, LF detections from 002/125/047, no HF
% #7 similar as 6, in region detections come from a later 002 migration, feels like location error in LF (002 ,125)
% #6 seems like true detections (002/047/125) with larger error from near 002 migration
% #5 feels like due to LF location error;
% #4 has HF but no LF, likely to be real with larger errors  
% 


%%% it may be not surprising that during any migration, at one region, LF/HF show up but the other
%%% don't, but at other times, the other show up. It may be because the coherency at that time is
%%% bad for one type. It is hard to rule out this possibility
%%% It is more interesting that one region ONLY produces one type of energy but no
%%% the other at all, the following examples belong to the first scenario.

% #8 seems like one where HF starts earlier, but surely the region generates LF at other times 
% #41 seems like one where HF starts earlier, but surely the region generates LF at other times 
% #43 seems like one where HF starts earlier, but surely the region generates LF at other times 
% #44 seems like one where HF starts earlier, but surely the region generates LF at other times 


%% plot the LF CC coeff. during time when HF is missing VS. when HF is present
for i = 1: length(indmig)
% for i = 1: 211
%     i=46;
    indhf = find(hftime(:,13)==trange(indmig(i),1) & hftime(:,15)>=trange(indmig(i),2) & ...
                 hftime(:,15)<=trange(indmig(i),3));
    mighf = hftime(indhf,:);
    
    indlf = find(lftime(:,13)==trange(indmig(i),1) & lftime(:,15)>=trange(indmig(i),2) & ...
                 lftime(:,15)<=trange(indmig(i),3));
    miglf = lftime(indlf,:); 

    if i < 3
        indlf = find(lftime(:,13)==trange(indmig(i),1) & lftime(:,15)>=trange1(indmig(i),2) & ...
                    lftime(:,15)<=trange(indmig(i),2));
        earlylf = lftime(indlf,:);
    else
        indlf = find(lftime(:,13)==trange(indmig(i),1) & lftime(:,15)>=trange(indmig(i),3) & ...
                    lftime(:,15)<=trange1(indmig(i),3));
        earlylf = lftime(indlf,:);
    end
%     earlylf = lfreg1mig{indmig(i), 1};
    
    
    %%% define and position the figure frame and axes of each plot
    f.fig=figure;
    widin = 8;  % maximum width allowed is 8.5 inches
    htin = 8;   % maximum height allowed is 11 inches
    set(f.fig,'Position',[scrsz(1)+1*scrsz(3)/10 scrsz(2)+scrsz(4)/20 widin*res htin*res]);
    nrow = 5;
    ncol = 1;    
    for isub = 1:nrow*ncol
        f.ax(isub) = subplot(nrow,ncol,isub);
    end
    
    % LZB trio 
    ax = f.ax(1);
    hold(ax,'on');
    ax.Box = 'on';
    grid(ax, 'on');
    plot(ax,earlylf(:,15)/3600, earlylf(:, 16), '.-','color','b'); 
    plot(ax,miglf(:,15)/3600, miglf(:, 16),  '.-','color','r');
    text(ax,0.9,0.8,'Trio','unit','normalized');
    
    % 4th station PGC
    ax = f.ax(2);
    hold(ax,'on');
    ax.Box = 'on';
    grid(ax, 'on');
    plot(ax,earlylf(:,15)/3600, earlylf(:, 17), '.-','color','b'); 
    plot(ax,miglf(:,15)/3600, miglf(:, 17),  '.-','color','r');
    text(ax,0.9,0.8,'PGC','unit','normalized');
    
    % 4th station SSIB
    ax = f.ax(3);
    hold(ax,'on');
    ax.Box = 'on';
    grid(ax, 'on');
    plot(ax,earlylf(:,15)/3600, earlylf(:, 18), '.-','color','b'); 
    plot(ax,miglf(:,15)/3600, miglf(:, 18),  '.-','color','r');
    text(ax,0.9,0.8,'SSIB','unit','normalized');
    
    % 4th station SILB
    ax = f.ax(4);
    hold(ax,'on');
    ax.Box = 'on';
    grid(ax, 'on');
    plot(ax,earlylf(:,15)/3600, earlylf(:, 19), '.-','color','b'); 
    plot(ax,miglf(:,15)/3600, miglf(:, 19),  '.-','color','r'); 
    text(ax,0.9,0.8,'SILB','unit','normalized');
    
    % 4th station KLNB
    ax = f.ax(5);
    hold(ax,'on');
    ax.Box = 'on';
    grid(ax, 'on');   
    plot(ax,earlylf(:,15)/3600, earlylf(:, 20), '.-','color','b'); 
    plot(ax,miglf(:,15)/3600, miglf(:, 20),  '.-','color','r'); 
    text(ax,0.9,0.8,'KLNB','unit','normalized');
    
end














