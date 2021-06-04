% function plot_all_rela_LZB_zoom
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function is used to plot the locations of all detections inverted by
% hypoinverse in km relative to offset (0,0), and w/ time
%
%   REMEMBER to check the input format !!!
%   Different from plot_dcut_rela_LZB.m which only plot the ones that pass
%   the distance cutoff and duplicates removal
%
% Looks like figure 4 in Peng et al. 2015
%
% Modified by Chao Song, chaosong@princeton.edu
% First created date:   2020/05/20
% Last modified date:   2020/06/04
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Initialization
format short e   % Set the format to 5-digit floating point
clear
close all
clc

set(0,'DefaultFigureVisible','on');
% set(0,'DefaultFigureVisible','off');   % switch to show the plots or not

% get the scrsz in pixels and number of pixels per inch of monitor 1
[scrsz, res] = pixelperinch(1);

workpath = getenv('MHOME');
rstpath = strcat(workpath, '/Seisbasics/hypoinverse/lzbrst');

winlenhf = 4;

winlenlf = 16;
lofflf = 4;
ccminlf = 0.35;

nfam = 12;

if nfam == 11
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
               -123.838500 48.544833 35.6600];
    
    % load the time results after dist cutoff and duplicates removal
    SUFFIXhf = strcat('up.hf.time.',num2str(winlenhf),'_',num2str(winlenlf),'.',num2str(lofflf),...
        '.',num2str(ccminlf));
    fname = strcat(rstpath, '/evtloc.allfam.nodcutnodou.',SUFFIXhf);
    hfrelatime = load(fname);
    % 25 cols, format is:
    %   E(043) N(043) E(own) N(own) dist(own) lon lat dep off12 off13 off12sec off13sec date
    %   main_arrival_time  cent_of_win  avecc ccadd1 ccadd2 ccadd3 ccadd4 ampmax ampmin eng normeng famnum
    
    SUFFIXlf = strcat('lf.time.',num2str(winlenhf),'_',num2str(winlenlf),'.',num2str(lofflf),...
        '.',num2str(ccminlf));
    fname = strcat(rstpath, '/evtloc.allfam.nodcutnodou.',SUFFIXlf);
    lfrelatime = load(fname);
    
elseif nfam == 12
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
               -123.879667 48.446167 34.2600];
    
    % load the time results after dist cutoff and duplicates removal
    SUFFIXhf = strcat('up.hf.time.',num2str(winlenhf),'_',num2str(winlenlf),'.',num2str(lofflf),...
        '.',num2str(ccminlf));
    fname = strcat(rstpath, '/evtloc.all12fam.nodcutnodou.',SUFFIXhf);
    hfrelatime = load(fname);
    % 25 cols, format is:
    %   E(043) N(043) E(own) N(own) dist(own) lon lat dep off12 off13 off12sec off13sec date
    %   main_arrival_time  cent_of_win  avecc ccadd1 ccadd2 ccadd3 ccadd4 ampmax ampmin eng normeng famnum
    
    SUFFIXlf = strcat('lf.time.',num2str(winlenhf),'_',num2str(winlenlf),'.',num2str(lofflf),...
        '.',num2str(ccminlf));
    fname = strcat(rstpath, '/evtloc.all12fam.nodcutnodou.',SUFFIXlf);
    lfrelatime = load(fname);
    
end

relacont = loccont;
[dx, dy] = absloc2relaloc(loccont(:,1),loccont(:,2),loccont(2,1),loccont(2,2));
relacont(:,1) = dx;
relacont(:,2) = dy;

% sort the results according to their absolute locations
hftemp = sortrows(hfrelatime, [6, 7]);
lftemp = sortrows(lfrelatime, [6, 7]);

% hf, obtain the counts at same position and combine them from diff fams
[hftempuniq,iuniq,~] = unique(hftemp(:,6:7),'rows','stable');
ctothf = zeros(size(hftempuniq,1),1);
medcchf = zeros(size(hftempuniq,1),5);
for i = 1:size(hftempuniq,1)
    [idup,~] = find(hftemp(:,6)==hftempuniq(i,1) & hftemp(:,7)==hftempuniq(i,2));
    ctothf(i,1) = length(idup);   % total count of hf
    medcchf(i,1:5) = median(hftemp(idup,16:20));   % this is median CC from trio and additional sta
end
if sum(ctothf) == size(hftemp,1)
    hfrelacount = [hftemp(iuniq,1:8) ctothf medcchf];
else
    disp('Inconsistency of hf total counts number')
end

% lf, obtain the counts at same position and combine them from diff fams
[lftempuniq,iuniq,~] = unique(lftemp(:,6:7),'rows','stable');
ctotlf = zeros(size(lftempuniq,1),1);
medcclf = zeros(size(lftempuniq,1),5);
for i = 1:size(lftempuniq,1)
    [idup,~] = find(lftemp(:,6)==lftempuniq(i,1) & lftemp(:,7)==lftempuniq(i,2));
    ctotlf(i,1) = length(idup);   % total count of hf
    medcclf(i,1:5) = median(lftemp(idup,16:20));   % this is median CC from trio and additional sta
end
if sum(ctotlf) == size(lftemp,1)
    lfrelacount = [lftemp(iuniq,1:8) ctotlf medcclf];   % this would give 8+1+5=14 cols
else
    disp('Inconsistency of lf total counts number')
end


%% read the detections checked by additional fam, choose KLNB only
if nfam == 11
    SUFFIXhf = strcat('up.hf.time.',num2str(winlenhf),'_',num2str(winlenlf),'.',num2str(lofflf),...
        '.',num2str(ccminlf));
    fname = strcat(rstpath, '/evtloc.KLNBcheck.allfam.nodcutnodou.',SUFFIXhf);
    hfaddcheck = load(fname);
    % 25 cols, format is:
    %   E(043) N(043) E(own) N(own) dist(own) lon lat dep off12 off13 off12sec off13sec date
    %   main_arrival_time  cent_of_win  avecc ccadd1 ccadd2 ccadd3 ccadd4 ampmax ampmin eng normeng famnum
    
    SUFFIXlf = strcat('lf.time.',num2str(winlenhf),'_',num2str(winlenlf),'.',num2str(lofflf),...
        '.',num2str(ccminlf));
    fname = strcat(rstpath, '/evtloc.KLNBcheck.allfam.nodcutnodou.',SUFFIXlf);
    lfaddcheck = load(fname);
    
elseif nfam == 12
    SUFFIXhf = strcat('up.hf.time.',num2str(winlenhf),'_',num2str(winlenlf),'.',num2str(lofflf),...
        '.',num2str(ccminlf));
    fname = strcat(rstpath, '/evtloc.KLNBcheck.all12fam.nodcutnodou.',SUFFIXhf);
    hfaddcheck = load(fname);
    % 25 cols, format is:
    %   E(043) N(043) E(own) N(own) dist(own) lon lat dep off12 off13 off12sec off13sec date
    %   main_arrival_time  cent_of_win  avecc ccadd1 ccadd2 ccadd3 ccadd4 ampmax ampmin eng normeng famnum
    
    SUFFIXlf = strcat('lf.time.',num2str(winlenhf),'_',num2str(winlenlf),'.',num2str(lofflf),...
        '.',num2str(ccminlf));
    fname = strcat(rstpath, '/evtloc.KLNBcheck.all12fam.nodcutnodou.',SUFFIXlf);
    lfaddcheck = load(fname);
end


%% plot the accumulative hit count
%%% figure 1, plot the overall hf and lf detections
% define and position the figure frame and axes of each plot
f.fig = figure;
widin = 8;  % maximum width allowed is 8.5 inches
htin = 10;   % maximum height allowed is 11 inches
set(f.fig,'Position',[1*scrsz(3)/20 scrsz(4)/10 widin*res htin*res]);
nrow = 2;
ncol = 2;
for isub = 1:nrow*ncol
    f.ax(isub) = subplot(nrow,ncol,isub);
    f.ax(isub).Box = 'on';
    grid(f.ax(isub), 'on');
    f.ax(isub).GridLineStyle = '--';
    axis(f.ax(isub), 'equal');
    axis(f.ax(isub),[-30 50 -30 40]);
    
end

% reposition
set(f.ax(3),'Position', [0.1, 0.1, 0.4, 0.4]);
set(f.ax(4),'Position', [0.55, 0.1, 0.4, 0.4]);
set(f.ax(1),'Position', [0.1, 0.55, 0.4, 0.4]);
set(f.ax(2),'Position', [0.55, 0.55, 0.4, 0.4]);

% marker size
msizehf = 1;
msizelf = 6;

% subplot 1 of figure 1
hold(f.ax(1),'on');
dumhf = hfrelacount;
dumhf(dumhf(:,9)>1, :) = [];
scatter(f.ax(1),dumhf(:,1),dumhf(:,2), msizehf, log10(dumhf(:,9)),'o','linew',0.2);  %, 'MarkerEdgeColor', 'w')
dumhf = sortrows(hfrelacount,9);
dumhf(dumhf(:,9)==1, :) = [];
scatter(f.ax(1),dumhf(:,1),dumhf(:,2), msizehf, log10(dumhf(:,9)), 'filled','o');  %, 'MarkerEdgeColor', 'w')
scatter(f.ax(1),relacont(:,1),relacont(:,2),20,'k','o','LineWidth',1);
plot(f.ax(1),[-100 100],[0 0],'k--');
plot(f.ax(1),[0 0],[-100 100],'k--');
text(f.ax(1),0.05,0.93,'a','FontSize',12,'unit','normalized','horizontalalignment','center',...
    'EdgeColor','k','Margin',2);
text(f.ax(1),0.93,0.93,'LZB','FontSize',12,'unit','normalized','horizontalalignment','center',...
    'EdgeColor','k','Margin',2);
text(f.ax(1),0.93,0.1,'HF','FontSize',15,'unit','normalized','horizontalalignment','center');
colormap(f.ax(1),'jet');
c=colorbar(f.ax(1),'SouthOutside');
c.Label.String = 'log_{10}(hits)';
c.Label.FontSize = 12;
% caxis(f.ax(1),[0 1.5]);
f.ax(1).YLabel.String = 'N (km)';
hold(f.ax(1),'off');

% subplot 2 of figure 1
hold(f.ax(2),'on');
dumlf = lfrelacount;
dumlf(dumlf(:,9)>1, :) = [];
scatter(f.ax(2),dumlf(:,1),dumlf(:,2), msizelf, log10(dumlf(:,9)),'o','linew',0.2);   %, 'MarkerEdgeColor', 'w')
dumlf = sortrows(lfrelacount,9);
dumlf(dumlf(:,9)==1, :) = [];
scatter(f.ax(2),dumlf(:,1),dumlf(:,2), msizelf, log10(dumlf(:,9)), 'filled','o');  %, 'MarkerEdgeColor', 'w')
scatter(f.ax(2),relacont(:,1),relacont(:,2),20,'k','o','LineWidth',1);
plot(f.ax(2),[-100 100],[0 0],'k--');
plot(f.ax(2),[0 0],[-100 100],'k--');
text(f.ax(2),0.05,0.93,'b','FontSize',12,'unit','normalized','horizontalalignment','center',...
    'EdgeColor','k','Margin',2);
text(f.ax(2),0.93,0.93,'LZB','FontSize',12,'unit','normalized','horizontalalignment','center',...
    'EdgeColor','k','Margin',2);
text(f.ax(2),0.93,0.1,'LF','FontSize',15,'unit','normalized','horizontalalignment','center');
colormap(f.ax(2),'jet');
c=colorbar(f.ax(2),'SouthOutside');
c.Label.String = 'log_{10}(hits)';
c.Label.FontSize = 12;
caxis(f.ax(2),[0 1.7]);
hold(f.ax(2),'off');

% subplot 3 of figure 1
hold(f.ax(3),'on');
dumhf = sortrows(hfrelacount,10);
scatter(f.ax(3),dumhf(:,1),dumhf(:,2), msizehf, dumhf(:,10), 'filled','o');  %, 'MarkerEdgeColor', 'w')
scatter(f.ax(3),relacont(:,1),relacont(:,2),20,'k','o','LineWidth',1);
plot(f.ax(3),[-100 100],[0 0],'k--');
plot(f.ax(3),[0 0],[-100 100],'k--');
text(f.ax(3),0.05,0.93,'c','FontSize',12,'unit','normalized','horizontalalignment','center',...
    'EdgeColor','k','Margin',2);
text(f.ax(3),0.93,0.93,'LZB','FontSize',12,'unit','normalized','horizontalalignment','center',...
    'EdgeColor','k','Margin',2);
text(f.ax(3),0.93,0.1,'HF','FontSize',15,'unit','normalized','horizontalalignment','center');
colormap(f.ax(3),'jet');
c=colorbar(f.ax(3),'NorthOutside');
c.Label.String = 'Median CC';
c.Label.FontSize = 12;
caxis(f.ax(3),[0.35 0.65]);
f.ax(3).XLabel.String = 'E (km)';
f.ax(3).YLabel.String = 'N (km)';
hold(f.ax(3),'off');

% subplot 4 of figure 1
hold(f.ax(4),'on');
dumlf = sortrows(lfrelacount,10);
scatter(f.ax(4),dumlf(:,1),dumlf(:,2), msizelf, dumlf(:,10), 'filled','o');  %, 'MarkerEdgeColor', 'w')
scatter(f.ax(4),relacont(:,1),relacont(:,2),20,'k','o','LineWidth',1);
plot(f.ax(4),[-100 100],[0 0],'k--');
plot(f.ax(4),[0 0],[-100 100],'k--');
text(f.ax(4),0.05,0.93,'d','FontSize',12,'unit','normalized','horizontalalignment','center',...
    'EdgeColor','k','Margin',2);
text(f.ax(4),0.93,0.93,'LZB','FontSize',12,'unit','normalized','horizontalalignment','center',...
    'EdgeColor','k','Margin',2);
text(f.ax(4),0.93,0.1,'LF','FontSize',15,'unit','normalized','horizontalalignment','center');
colormap(f.ax(4),'jet');
c=colorbar(f.ax(4),'NorthOutside');
c.Label.String = 'Median CC';
c.Label.FontSize = 12;
caxis(f.ax(4),[0.35 0.65]);
f.ax(4).XLabel.String = 'E (km)';
hold(f.ax(4),'off');

supertit(f.ax,strcat('absloc--no dist cutoff'));

% print(f.fig,'-dpdf',strcat(rstpath,'/asd.rela.pdf'));

%%% save figure to file
if nfam == 11
    print(f.fig,'-dpdf',strcat(rstpath,'/allfam.hflf.alld.nodcut.accuhit.',num2str(winlenhf),'_',...
          num2str(winlenlf),'.',num2str(lofflf),'.',num2str(ccminlf),'.rela.pdf'));
elseif nfam == 12
    print(f.fig,'-dpdf',strcat(rstpath,'/all12fam.hflf.alld.nodcut.accuhit.',num2str(winlenhf),'_',...
          num2str(winlenlf),'.',num2str(lofflf),'.',num2str(ccminlf),'.rela.pdf'));
end
% keyboard



%% it seems that LF & HF have different density patterns, so zoom in and check timing there
reg1 = [-20 -20;
    0 -20;
    0   0;
    -20   0;
    -20 -20];

reg2 = [ 5 -20;
    25 -20;
    25   0;
    5   0;
    5 -20];

% events & tremor inside bound 1
[is,ion] = inpolygon(hfrelatime(:,1),hfrelatime(:,2),reg1(:,1),reg1(:,2));
isinreg1 = is | ion;
hfreg1 = hfrelatime(isinreg1 == 1, :);

[is,ion] = inpolygon(lfrelatime(:,1),lfrelatime(:,2),reg1(:,1),reg1(:,2));
isinreg1 = is | ion;
lfreg1 = lfrelatime(isinreg1 == 1, :);

% events & tremor inside bound 2
[is,ion] = inpolygon(hfrelatime(:,1),hfrelatime(:,2),reg2(:,1),reg2(:,2));
isinreg2 = is | ion;
hfreg2 = hfrelatime(isinreg2 == 1, :);

[is,ion] = inpolygon(lfrelatime(:,1),lfrelatime(:,2),reg2(:,1),reg2(:,2));
isinreg2 = is | ion;
lfreg2 = lfrelatime(isinreg2 == 1, :);

% % rotate to the strike N30W and dip N60E direction
% rotang = 30;  % counter-clockwise from y to strike
% % now the 1,2 col is changing from E(043) and N(043), to down-dip and strike
% [hfreg1(:,1),hfreg1(:,2)] = coordinate_rot(hfreg1(:,1),hfreg1(:,2),rotang,0,0);
% [lfreg1(:,1),lfreg1(:,2)] = coordinate_rot(lfreg1(:,1),lfreg1(:,2),rotang,0,0);
% [hfreg2(:,1),hfreg2(:,2)] = coordinate_rot(hfreg2(:,1),hfreg2(:,2),rotang,0,0);
% [lfreg2(:,1),lfreg2(:,2)] = coordinate_rot(lfreg2(:,1),lfreg2(:,2),rotang,0,0);

% sort them according to their occurring time
hfreg1 = sortrows(hfreg1,[13,15]);
hfreg2 = sortrows(hfreg2,[13,15]);
lfreg1 = sortrows(lfreg1,[13,15]);
lfreg2 = sortrows(lfreg2,[13,15]);


dateall = unique(hfrelatime(:,13));

% for 2005 only, rotate to the strike N30W and dip N60E direction, 2003 and 2004 have the main front
% migration direction to N
rotang = 30;  % counter-clockwise from y to strike, 2005

% obtain the relative time in hr
ih103 = find(hfreg1(:,13) < 2004*1000);
hfreg1(ih103,26) = (hfreg1(ih103,13)-2003060).*24+hfreg1(ih103,15)./3600;
ih104 = find(hfreg1(:,13) < 2005*1000 & hfreg1(:,13) > 2004*1000);
hfreg1(ih104,26) = (hfreg1(ih104,13)-2004194).*24+hfreg1(ih104,15)./3600;
ih105 = find(hfreg1(:,13) > 2005*1000);
hfreg1(ih105,26) = (hfreg1(ih105,13)-2005254).*24+hfreg1(ih105,15)./3600;
% rotate to the strike N30W and dip N60E direction,now the 1,2 col is changing from E(043) and N(043),
% to down-dip and strike
[hfreg1(ih105,1),hfreg1(ih105,2)] = coordinate_rot(hfreg1(ih105,1),hfreg1(ih105,2),rotang,0,0);

il103 = find(lfreg1(:,13) < 2004*1000);
lfreg1(il103,26) = (lfreg1(il103,13)-2003060).*24+lfreg1(il103,15)./3600;
il104 = find(lfreg1(:,13) < 2005*1000 & lfreg1(:,13) > 2004*1000);
lfreg1(il104,26) = (lfreg1(il104,13)-2004194).*24+lfreg1(il104,15)./3600;
il105 = find(lfreg1(:,13) > 2005*1000);
lfreg1(il105,26) = (lfreg1(il105,13)-2005254).*24+lfreg1(il105,15)./3600;
% rotate to the strike N30W and dip N60E direction,now the 1,2 col is changing from E(043) and N(043),
% to down-dip and strike
[lfreg1(il105,1),lfreg1(il105,2)] = coordinate_rot(lfreg1(il105,1),lfreg1(il105,2),rotang,0,0);

ih203 = find(hfreg2(:,13) < 2004*1000);
hfreg2(ih203,26) = (hfreg2(ih203,13)-2003060).*24+hfreg2(ih203,15)./3600;
ih204 = find(hfreg2(:,13) < 2005*1000 & hfreg2(:,13) > 2004*1000);
hfreg2(ih204,26) = (hfreg2(ih204,13)-2004194).*24+hfreg2(ih204,15)./3600;
ih205 = find(hfreg2(:,13) > 2005*1000);
hfreg2(ih205,26) = (hfreg2(ih205,13)-2005254).*24+hfreg2(ih205,15)./3600;
% rotate to the strike N30W and dip N60E direction,now the 1,2 col is changing from E(043) and N(043),
% to down-dip and strike
[hfreg2(ih205,1),hfreg2(ih205,2)] = coordinate_rot(hfreg2(ih205,1),hfreg2(ih205,2),rotang,0,0);

il203 = find(lfreg2(:,13) < 2004*1000);
lfreg2(il203,26) = (lfreg2(il203,13)-2003060).*24+lfreg2(il203,15)./3600;
il204 = find(lfreg2(:,13) < 2005*1000 & lfreg2(:,13) > 2004*1000);
lfreg2(il204,26) = (lfreg2(il204,13)-2004194).*24+lfreg2(il204,15)./3600;
il205 = find(lfreg2(:,13) > 2005*1000);
lfreg2(il205,26) = (lfreg2(il205,13)-2005254).*24+lfreg2(il205,15)./3600;
% rotate to the strike N30W and dip N60E direction,now the 1,2 col is changing from E(043) and N(043),
% to down-dip and strike
[lfreg2(il205,1),lfreg2(il205,2)] = coordinate_rot(lfreg2(il205,1),lfreg2(il205,2),rotang,0,0);


%%% for detections passed the check by KLNB
% tremors that passed the additional check by KLNB
[is,ion] = inpolygon(hfaddcheck(:,1),hfaddcheck(:,2),reg1(:,1),reg1(:,2));
isinreg1 = is | ion;
hfareg1 = hfaddcheck(isinreg1 == 1, :);

[is,ion] = inpolygon(lfaddcheck(:,1),lfaddcheck(:,2),reg1(:,1),reg1(:,2));
isinreg1 = is | ion;
lfareg1 = lfaddcheck(isinreg1 == 1, :);

% events & tremor inside bound 2
[is,ion] = inpolygon(hfaddcheck(:,1),hfaddcheck(:,2),reg2(:,1),reg2(:,2));
isinreg2 = is | ion;
hfareg2 = hfaddcheck(isinreg2 == 1, :);

[is,ion] = inpolygon(lfaddcheck(:,1),lfaddcheck(:,2),reg2(:,1),reg2(:,2));
isinreg2 = is | ion;
lfareg2 = lfaddcheck(isinreg2 == 1, :);

% sort based on time 
hfareg1 = sortrows(hfareg1,[13,15]);
hfareg2 = sortrows(hfareg2,[13,15]);
lfareg1 = sortrows(lfareg1,[13,15]);
lfareg2 = sortrows(lfareg2,[13,15]);

% obtain the relative time in hr
iha103 = find(hfareg1(:,13) < 2004*1000);
hfareg1(iha103,23) = (hfareg1(iha103,13)-2003060).*24+hfareg1(iha103,15)./3600;
iha104 = find(hfareg1(:,13) < 2005*1000 & hfareg1(:,13) > 2004*1000);
hfareg1(iha104,23) = (hfareg1(iha104,13)-2004194).*24+hfareg1(iha104,15)./3600;
iha105 = find(hfareg1(:,13) > 2005*1000);
hfareg1(iha105,23) = (hfareg1(iha105,13)-2005254).*24+hfareg1(iha105,15)./3600;
% rotate to the strike N30W and dip N60E direction,now the 1,2 col is changing from E(043) and N(043),
% to down-dip and strike
[hfareg1(iha105,1),hfareg1(iha105,2)] = coordinate_rot(hfareg1(iha105,1),hfareg1(iha105,2),rotang,0,0);

ila103 = find(lfareg1(:,13) < 2004*1000);
lfareg1(ila103,23) = (lfareg1(ila103,13)-2003060).*24+lfareg1(ila103,15)./3600;
ila104 = find(lfareg1(:,13) < 2005*1000 & lfareg1(:,13) > 2004*1000);
lfareg1(ila104,23) = (lfareg1(ila104,13)-2004194).*24+lfareg1(ila104,15)./3600;
ila105 = find(lfareg1(:,13) > 2005*1000);
lfareg1(ila105,23) = (lfareg1(ila105,13)-2005254).*24+lfareg1(ila105,15)./3600;
% rotate to the strike N30W and dip N60E direction,now the 1,2 col is changing from E(043) and N(043),
% to down-dip and strike
[lfareg1(ila105,1),lfareg1(ila105,2)] = coordinate_rot(lfareg1(ila105,1),lfareg1(ila105,2),rotang,0,0);

iha203 = find(hfareg2(:,13) < 2004*1000);
hfareg2(iha203,23) = (hfareg2(iha203,13)-2003060).*24+hfareg2(iha203,15)./3600;
iha204 = find(hfareg2(:,13) < 2005*1000 & hfareg2(:,13) > 2004*1000);
hfareg2(iha204,23) = (hfareg2(iha204,13)-2004194).*24+hfareg2(iha204,15)./3600;
iha205 = find(hfareg2(:,13) > 2005*1000);
hfareg2(iha205,23) = (hfareg2(iha205,13)-2005254).*24+hfareg2(iha205,15)./3600;
% rotate to the strike N30W and dip N60E direction,now the 1,2 col is changing from E(043) and N(043),
% to down-dip and strike
[hfareg2(iha205,1),hfareg2(iha205,2)] = coordinate_rot(hfareg2(iha205,1),hfareg2(iha205,2),rotang,0,0);

ila203 = find(lfareg2(:,13) < 2004*1000);
lfareg2(ila203,23) = (lfareg2(ila203,13)-2003060).*24+lfareg2(ila203,15)./3600;
ila204 = find(lfareg2(:,13) < 2005*1000 & lfareg2(:,13) > 2004*1000);
lfareg2(ila204,23) = (lfareg2(ila204,13)-2004194).*24+lfareg2(ila204,15)./3600;
ila205 = find(lfareg2(:,13) > 2005*1000);
lfareg2(ila205,23) = (lfareg2(ila205,13)-2005254).*24+lfareg2(ila205,15)./3600;
% rotate to the strike N30W and dip N60E direction,now the 1,2 col is changing from E(043) and N(043),
% to down-dip and strike
[lfareg2(ila205,1),lfareg2(ila205,2)] = coordinate_rot(lfareg2(ila205,1),lfareg2(ila205,2),rotang,0,0);



%% region 1, along strike
f4.fig = figure;
f4.fig.Renderer = 'painters';
widin = 12;  % maximum width allowed is 8.5 inches
htin = 9.5;   % maximum height allowed is 11 inches
set(f4.fig,'Position',[1*scrsz(3)/20 scrsz(4)/10 widin*res htin*res]);
nrow = 6;
ncol = 1;
for isub = 1:nrow*ncol
    f4.ax(isub) = subplot(nrow,ncol,isub);
end

% reposition
for isub = 1:nrow*ncol
    f4.ax(isub).Box = 'on';
    grid(f4.ax(isub), 'on');
    f4.ax(isub).GridLineStyle = '--';
end

% marker size
msizehf = 2;

% ax = f4.ax(1);
% hold(ax,'on');
% scatter(ax,hfreg1(ih103,24),hfreg1(ih103,2), msizehf, 'filled','o');  %, 'MarkerEdgeColor', 'w')
% scatter(ax,lfreg1(il103,24),lfreg1(il103,2), msizehf, 'filled','ro');  %, 'MarkerEdgeColor', 'w'))
% text(ax,0.05,0.93,'2003','FontSize',10,'unit','normalized','horizontalalignment','left',...
%     'EdgeColor','k','Margin',2);
% axis(ax,[0 96 -15 5]);
% ax.YLabel.String = 'Along strike (km)';
% hold(ax,'off');
% 
% ax = f4.ax(2);
% hold(ax,'on');
% scatter(ax,hfreg1(ih103,24),hfreg1(ih103,2),msizehf,hfreg1(ih103,21)-hfreg1(ih103,22),'filled','o');  %, 'MarkerEdgeColor', 'w')
% text(ax,0.05,0.93,{'2003-hf diff amp'},'FontSize',10,'unit','normalized','horizontalalignment','left',...
%     'EdgeColor','k','Margin',2);
% axis(ax,[0 96 -15 5]);
% colormap(ax,'jet');
% colorbar(ax);
% hold(ax,'off');
% 
% ax = f4.ax(3);
% hold(ax,'on');
% scatter(ax,lfreg1(il103,24),lfreg1(il103,2),msizehf,median(lfreg1(il103,17:20),2),'filled','o');  %, 'MarkerEdgeColor', 'w'))
% text(ax,0.05,0.93,{'2003-lf med add CC'},'FontSize',10,'unit','normalized','horizontalalignment','left',...
%     'EdgeColor','k','Margin',2);
% axis(ax,[0 96 -15 5]);
% colormap(ax,'jet');
% colorbar(ax);
% hold(ax,'off');
% 
% ax = f4.ax(4);
% hold(ax,'on');
% scatter(ax,hfreg1(ih104,24),hfreg1(ih104,2), msizehf, 'filled','o');  %, 'MarkerEdgeColor', 'w')
% scatter(ax,lfreg1(il104,24),lfreg1(il104,2), msizehf, 'filled','ro');  %, 'MarkerEdgeColor', 'w'))
% text(ax,0.05,0.93,'2004','FontSize',10,'unit','normalized','horizontalalignment','left',...
%     'EdgeColor','k','Margin',2);
% axis(ax,[0 96 -15 5]);
% hold(ax,'off');
% 
% ax = f4.ax(5);
% hold(ax,'on');
% scatter(ax,hfreg1(ih104,24),hfreg1(ih104,2),msizehf,hfreg1(ih104,21)-hfreg1(ih104,22),'filled','o');  %, 'MarkerEdgeColor', 'w')
% text(ax,0.05,0.93,{'2004-hf diff amp'},'FontSize',10,'unit','normalized','horizontalalignment','left',...
%     'EdgeColor','k','Margin',2);
% axis(ax,[0 96 -15 5]);
% colormap(ax,'jet');
% colorbar(ax);
% hold(ax,'off');
% 
% ax = f4.ax(6);
% hold(ax,'on');
% scatter(ax,lfreg1(il104,24),lfreg1(il104,2),msizehf,median(lfreg1(il104,17:20),2),'filled','o');  %, 'MarkerEdgeColor', 'w'))
% text(ax,0.05,0.93,{'2004-lf med add CC'},'FontSize',10,'unit','normalized','horizontalalignment','left',...
%     'EdgeColor','k','Margin',2);
% axis(ax,[0 96 -15 5]);
% colormap(ax,'jet');
% colorbar(ax);
% hold(ax,'off');[0.6 0.6 0.6]

ax = f4.ax(1);
hold(ax,'on');
scatter(ax,hfreg1(ih105,26),hfreg1(ih105,2),msizehf,'filled','ko');  %, 'MarkerEdgeColor', 'w')
scatter(ax,lfreg1(il105,26),lfreg1(il105,2),msizehf,'filled','ro');  %, 'MarkerEdgeColor', 'w'))
text(ax,0.05,0.93,'2005-all detections','FontSize',10,'unit','normalized','horizontalalignment','left',...
    'EdgeColor','k','Margin',2);
axis(ax,[0 60 -15 5]);
ax.YLabel.String = 'Along strike (km)';
hold(ax,'off');

ax = f4.ax(2);
hold(ax,'on');
scatter(ax,hfreg1(ih105,26),hfreg1(ih105,2),msizehf,log(hfreg1(ih105,23)),'filled','o');  %, 'MarkerEdgeColor', 'w')
text(ax,0.05,0.93,{'2005-hf energy'},'FontSize',10,'unit','normalized','horizontalalignment','left',...
    'EdgeColor','k','Margin',2);
axis(ax,[0 60 -15 5]);
colormap(ax,'jet');
colorbar(ax);
cran = [prctile(log(hfreg1(ih105,23)),2), prctile(log(hfreg1(ih105,23)),98)];
caxis(ax,cran);
% caxis(ax,[0.1 0.7]);
hold(ax,'off');

ax = f4.ax(3);
hold(ax,'on');
scatter(ax,lfreg1(il105,26),lfreg1(il105,2),msizehf,log(lfreg1(il105,23)),'filled','o');  %, 'MarkerEdgeColor', 'w')
text(ax,0.05,0.93,{'2005-lf energy'},'FontSize',10,'unit','normalized','horizontalalignment','left',...
    'EdgeColor','k','Margin',2);
axis(ax,[0 60 -15 5]);
colormap(ax,'jet');
colorbar(ax);
cran = [prctile(log(lfreg1(il105,23)),2), prctile(log(lfreg1(il105,23)),98)];
caxis(ax,cran);
% caxis(ax,[0.1 0.5]);
hold(ax,'off');

ax = f4.ax(4);
hold(ax,'on');
scatter(ax,hfareg1(iha105,23),hfareg1(iha105,2),msizehf,hfareg1(iha105,21),'filled','o');  %, 'MarkerEdgeColor', 'w'))
text(ax,0.05,0.93,{'2005-hf KLNB add CC'},'FontSize',10,'unit','normalized','horizontalalignment','left',...
    'EdgeColor','k','Margin',2);
axis(ax,[0 60 -15 5]);
colormap(ax,'jet');
colorbar(ax);
cran = [prctile(hfareg1(iha105,21),2), prctile(hfareg1(iha105,21),98)];
caxis(ax,cran);
% caxis(ax,[0.35 0.55]);
hold(ax,'off');

ax = f4.ax(5);
hold(ax,'on');
scatter(ax,lfareg1(ila105,23),lfareg1(ila105,2),msizehf,lfareg1(ila105,21),'filled','o');  %, 'MarkerEdgeColor', 'w'))
text(ax,0.05,0.93,{'2005-lf KLNB add CC'},'FontSize',10,'unit','normalized','horizontalalignment','left',...
    'EdgeColor','k','Margin',2);
axis(ax,[0 60 -15 5]);
colormap(ax,'jet');
colorbar(ax);
cran = [prctile(lfareg1(ila105,21),2), prctile(lfareg1(ila105,21),98)];
caxis(ax,cran);
% caxis(ax,[0.35 0.55]);
hold(ax,'off');

ax = f4.ax(6);
hold(ax,'on');
scatter(ax,hfareg1(iha105,23),hfareg1(iha105,2),msizehf,'filled','ko');  %, 'MarkerEdgeColor', 'w')
scatter(ax,lfareg1(ila105,23),lfareg1(ila105,2),msizehf,'filled','ro');  %, 'MarkerEdgeColor', 'w'))
text(ax,0.05,0.93,'2005-KLNB checked','FontSize',10,'unit','normalized','horizontalalignment','left',...
    'EdgeColor','k','Margin',2);
axis(ax,[0 60 -15 5]);
ax.XLabel.String = 'Time (hr)';
hold(ax,'off');

print(f4.fig,'-depsc2',strcat(rstpath,'/reg1strike2005.eps'));



%% region 2, along strike
f6.fig = figure;
f6.fig.Renderer = 'painters';
widin = 12;  % maximum width allowed is 8.5 inches
htin = 9.5;   % maximum height allowed is 11 inches
set(f6.fig,'Position',[1*scrsz(3)/20 scrsz(4)/10 widin*res htin*res]);
nrow = 6;
ncol = 1;
for isub = 1:nrow*ncol
    f6.ax(isub) = subplot(nrow,ncol,isub);
end

% reposition
for isub = 1:nrow*ncol
    f6.ax(isub).Box = 'on';
    grid(f6.ax(isub), 'on');
    f6.ax(isub).GridLineStyle = '--';
end

% marker size
msizehf = 2;

% ax = f6.ax(1);
% hold(ax,'on');
% scatter(ax,hfreg2(ih203,24),hfreg2(ih203,2), msizehf, 'filled','o');  %, 'MarkerEdgeColor', 'w')
% scatter(ax,lfreg2(il203,24),lfreg2(il203,2), msizehf, 'filled','ro');  %, 'MarkerEdgeColor', 'w'))
% text(ax,0.05,0.93,{'2003-detections'},'FontSize',10,'unit','normalized','horizontalalignment','center',...
%     'EdgeColor','k','Margin',2);
% axis(ax,[0 96 -29 -2]);
% ax.YLabel.String = 'Along strike (km)';
% hold(ax,'off');
% 
% ax = f6.ax(2);
% hold(ax,'on');
% scatter(ax,hfreg2(ih203,24),hfreg2(ih203,2),msizehf,hfreg2(ih203,21)-hfreg2(ih203,22),'filled','o');  %, 'MarkerEdgeColor', 'w')
% text(ax,0.05,0.93,{'2003-hf diff amp'},'FontSize',10,'unit','normalized','horizontalalignment','center',...
%     'EdgeColor','k','Margin',2);
% axis(ax,[0 96 -29 -2]);
% colormap jet
% colorbar
% hold(ax,'off');
% 
% ax = f6.ax(3);
% hold(ax,'on');
% scatter(ax,lfreg2(il203,24),lfreg2(il203,2),msizehf,median(lfreg2(il203,17:20),2),'filled','o');  %, 'MarkerEdgeColor', 'w'))
% text(ax,0.05,0.93,{'2003-lf med add CC'},'FontSize',10,'unit','normalized','horizontalalignment','center',...
%     'EdgeColor','k','Margin',2);
% axis(ax,[0 96 -29 -2]);
% colormap jet
% colorbar
% hold(ax,'off');
% 
% ax = f6.ax(4);
% hold(ax,'on');
% scatter(ax,hfreg2(ih204,24),hfreg2(ih204,2), msizehf, 'filled','o');  %, 'MarkerEdgeColor', 'w')
% scatter(ax,lfreg2(il204,24),lfreg2(il204,2), msizehf, 'filled','ro');  %, 'MarkerEdgeColor', 'w'))
% text(ax,0.05,0.93,'2004','FontSize',13,'unit','normalized','horizontalalignment','center',...
%     'EdgeColor','k','Margin',2);
% axis(ax,[0 96 -29 -2]);
% hold(ax,'off');
% 
% ax = f6.ax(5);
% hold(ax,'on');
% scatter(ax,hfreg2(ih204,24),hfreg2(ih204,2),msizehf,hfreg2(ih204,21)-hfreg2(ih204,22),'filled','o');  %, 'MarkerEdgeColor', 'w')
% text(ax,0.05,0.93,{'2004-hf diff amp'},'FontSize',10,'unit','normalized','horizontalalignment','center',...
%     'EdgeColor','k','Margin',2);
% axis(ax,[0 96 -29 -2]);
% colormap jet
% colorbar
% hold(ax,'off');
% 
% ax = f6.ax(6);
% hold(ax,'on');
% scatter(ax,lfreg2(il204,24),lfreg2(il204,2),msizehf,median(lfreg2(il204,17:20),2),'filled','o');  %, 'MarkerEdgeColor', 'w'))
% text(ax,0.05,0.93,{'2004-lf med add CC'},'FontSize',10,'unit','normalized','horizontalalignment','center',...
%     'EdgeColor','k','Margin',2);
% axis(ax,[0 96 -29 -2]);
% colormap jet
% colorbar
% hold(ax,'off');

ax = f6.ax(1);
hold(ax,'on');
scatter(ax,hfreg2(ih205,26),hfreg2(ih205,2), msizehf, 'filled','ko');  %, 'MarkerEdgeColor', 'w')
scatter(ax,lfreg2(il205,26),lfreg2(il205,2), msizehf, 'filled','ro');  %, 'MarkerEdgeColor', 'w'))
text(ax,0.05,0.93,'2005-all detections','FontSize',10,'unit','normalized','horizontalalignment','left',...
    'EdgeColor','k','Margin',2);
axis(ax,[0 60 -29 -2]);
ax.YLabel.String = 'Along strike (km)';
hold(ax,'off');

ax = f6.ax(2);
hold(ax,'on');
scatter(ax,hfreg2(ih205,26),hfreg2(ih205,2),msizehf,log(hfreg2(ih205,23)),'filled','o');  %, 'MarkerEdgeColor', 'w')
text(ax,0.05,0.93,{'2005-hf energy'},'FontSize',10,'unit','normalized','horizontalalignment','left',...
    'EdgeColor','k','Margin',2);
axis(ax,[0 60 -29 -2]);
colormap(ax,'jet');
colorbar(ax);
cran = [prctile(log(hfreg2(ih205,23)),2), prctile(log(hfreg2(ih205,23)),98)];
caxis(ax,cran);
% caxis(ax,[0.1 0.7]);
hold(ax,'off');

ax = f6.ax(3);
hold(ax,'on');
scatter(ax,lfreg2(il205,26),lfreg2(il205,2),msizehf,log(lfreg2(il205,23)),'filled','o');  %, 'MarkerEdgeColor', 'w')
text(ax,0.05,0.93,{'2005-lf energy'},'FontSize',10,'unit','normalized','horizontalalignment','left',...
    'EdgeColor','k','Margin',2);
axis(ax,[0 60 -29 -2]);
colormap(ax,'jet');
colorbar(ax);
cran = [prctile(log(lfreg2(il205,23)),2), prctile(log(lfreg2(il205,23)),98)];
caxis(ax,cran);
% caxis(ax,[0.1 0.5]);
hold(ax,'off');

ax = f6.ax(4);
hold(ax,'on');
scatter(ax,hfareg2(iha205,23),hfareg2(iha205,2),msizehf,hfareg2(iha205,21),'filled','o');  %, 'MarkerEdgeColor', 'w'))
text(ax,0.05,0.93,{'2005-hf KLNB add CC'},'FontSize',10,'unit','normalized','horizontalalignment','left',...
    'EdgeColor','k','Margin',2);
axis(ax,[0 60 -29 -2]);
colormap(ax,'jet');
colorbar(ax);
cran = [prctile(hfareg2(iha205,21),2), prctile(hfareg2(iha205,21),98)];
caxis(ax,cran);
% caxis(ax,[0.35 0.55]);
hold(ax,'off');

ax = f6.ax(5);
hold(ax,'on');
scatter(ax,lfareg2(ila205,23),lfareg2(ila205,2),msizehf,lfareg2(ila205,21),'filled','o');  %, 'MarkerEdgeColor', 'w'))
text(ax,0.05,0.93,{'2005-lf KLNB add CC'},'FontSize',10,'unit','normalized','horizontalalignment','left',...
    'EdgeColor','k','Margin',2);
axis(ax,[0 60 -29 -2]);
colormap(ax,'jet');
colorbar(ax);
cran = [prctile(lfareg2(ila205,21),2), prctile(lfareg2(ila205,21),98)];
caxis(ax,cran);
% caxis(ax,[0.35 0.55]);
hold(ax,'off');

ax = f6.ax(6);
hold(ax,'on');
scatter(ax,hfareg2(iha205,23),hfareg2(iha205,2), msizehf, 'filled','ko');  %, 'MarkerEdgeColor', 'w')
scatter(ax,lfareg2(ila205,23),lfareg2(ila205,2), msizehf, 'filled','ro');  %, 'MarkerEdgeColor', 'w'))
text(ax,0.05,0.93,'2005-KLNB checked','FontSize',10,'unit','normalized','horizontalalignment','left',...
    'EdgeColor','k','Margin',2);
axis(ax,[0 60 -29 -2]);
ax.XLabel.String = 'Time (hr)';
hold(ax,'off');

print(f6.fig,'-depsc2',strcat(rstpath,'/reg2strike2005.eps'));



%% try the entire catalog
hfall = hfrelatime;
lfall = lfrelatime;

% sort them according to their occurring time
hfall = sortrows(hfall,[13,15]);
lfall = sortrows(lfall,[13,15]);

% obtain the relative time in hr
ih03 = find(hfall(:,13) < 2004*1000);
hfall(ih03,26) = (hfall(ih03,13)-2003060).*24+hfall(ih03,15)./3600;
ih04 = find(hfall(:,13) < 2005*1000 & hfall(:,13) > 2004*1000);
hfall(ih04,26) = (hfall(ih04,13)-2004194).*24+hfall(ih04,15)./3600;
ih05 = find(hfall(:,13) > 2005*1000);
hfall(ih05,26) = (hfall(ih05,13)-2005254).*24+hfall(ih05,15)./3600;
% rotate to the strike N30W and dip N60E direction,now the 1,2 col is changing from E(043) and N(043),
% to down-dip and strike
[hfall(ih05,1),hfall(ih05,2)] = coordinate_rot(hfall(ih05,1),hfall(ih05,2),rotang,0,0);

il03 = find(lfall(:,13) < 2004*1000);
lfall(il03,26) = (lfall(il03,13)-2003060).*24+lfall(il03,15)./3600;
il04 = find(lfall(:,13) < 2005*1000 & lfall(:,13) > 2004*1000);
lfall(il04,26) = (lfall(il04,13)-2004194).*24+lfall(il04,15)./3600;
il05 = find(lfall(:,13) > 2005*1000);
lfall(il05,26) = (lfall(il05,13)-2005254).*24+lfall(il05,15)./3600;
% rotate to the strike N30W and dip N60E direction,now the 1,2 col is changing from E(043) and N(043),
% to down-dip and strike
[lfall(il05,1),lfall(il05,2)] = coordinate_rot(lfall(il05,1),lfall(il05,2),rotang,0,0);



%%% for the detections that pass check by KLNB
hfaall = hfaddcheck;
lfaall = lfaddcheck;

% sort them according to their occurring time
hfaall = sortrows(hfaall,[13,15]);
lfaall = sortrows(lfaall,[13,15]);

% obtain the relative time in hr
iha03 = find(hfaall(:,13) < 2004*1000);
hfaall(iha03,23) = (hfaall(iha03,13)-2003060).*24+hfaall(iha03,15)./3600;
iha04 = find(hfaall(:,13) < 2005*1000 & hfaall(:,13) > 2004*1000);
hfaall(iha04,23) = (hfaall(iha04,13)-2004194).*24+hfaall(iha04,15)./3600;
iha05 = find(hfaall(:,13) > 2005*1000);
hfaall(iha05,23) = (hfaall(iha05,13)-2005254).*24+hfaall(iha05,15)./3600;
% rotate to the strike N30W and dip N60E direction,now the 1,2 col is changing from E(043) and N(043),
% to down-dip and strike
[hfaall(iha05,1),hfaall(iha05,2)] = coordinate_rot(hfaall(iha05,1),hfaall(iha05,2),rotang,0,0);

ila03 = find(lfaall(:,13) < 2004*1000);
lfaall(ila03,23) = (lfaall(ila03,13)-2003060).*24+lfaall(ila03,15)./3600;
ila04 = find(lfaall(:,13) < 2005*1000 & lfaall(:,13) > 2004*1000);
lfaall(ila04,23) = (lfaall(ila04,13)-2004194).*24+lfaall(ila04,15)./3600;
ila05 = find(lfaall(:,13) > 2005*1000);
lfaall(ila05,23) = (lfaall(ila05,13)-2005254).*24+lfaall(ila05,15)./3600;
% rotate to the strike N30W and dip N60E direction,now the 1,2 col is changing from E(043) and N(043),
% to down-dip and strike
[lfaall(ila05,1),lfaall(ila05,2)] = coordinate_rot(lfaall(ila05,1),lfaall(ila05,2),rotang,0,0);


%% all catalog, along strike, 2004
f8.fig = figure;
f8.fig.Renderer = 'painters';
widin = 12;  % maximum width allowed is 8.5 inches
htin = 9.5;   % maximum height allowed is 11 inches
set(f8.fig,'Position',[1*scrsz(3)/20 scrsz(4)/10 widin*res htin*res]);
nrow = 6;
ncol = 1;
for isub = 1:nrow*ncol
    f8.ax(isub) = subplot(nrow,ncol,isub);
end

% reposition
for isub = 1:nrow*ncol
    f8.ax(isub).Box = 'on';
    grid(f8.ax(isub), 'on');
    f8.ax(isub).GridLineStyle = '--';
end

% marker size
msizehf = 2;

% ax = f8.ax(1);
% hold(ax,'on');
% scatter(ax,hfall(ih03,24),hfall(ih03,2), msizehf, 'filled','ko');  %, 'MarkerEdgeColor', 'w')
% scatter(ax,lfall(il03,24),lfall(il03,2), msizehf, 'filled','ro');  %, 'MarkerEdgeColor', 'w'))
% text(ax,0.05,0.93,'2003-detections','FontSize',10,'unit','normalized','horizontalalignment','left',...
%     'EdgeColor','k','Margin',2);
% axis(ax,[0 96 -26 18]);
% ax.YLabel.String = 'Along strike (km)';
% hold(ax,'off');
% 
% ax = f8.ax(2);
% hold(ax,'on');
% scatter(ax,hfall(ih03,24),hfall(ih03,2),msizehf,hfall(ih03,21)-hfall(ih03,22),'filled','o');  %, 'MarkerEdgeColor', 'w')
% text(ax,0.05,0.93,{'2003-hf diff amp'},'FontSize',10,'unit','normalized','horizontalalignment','left',...
%     'EdgeColor','k','Margin',2);
% axis(ax,[0 96 -26 18]);
% colormap(ax,'jet');
% colorbar(ax);
% caxis(ax,[0.1 0.7]);
% hold(ax,'off');
% 
% ax = f8.ax(3);
% hold(ax,'on');
% scatter(ax,lfall(il03,24),lfall(il03,2),msizehf,median(lfall(il03,17:20),2),'filled','o');  %, 'MarkerEdgeColor', 'w'))
% text(ax,0.05,0.93,{'2003-lf med add CC'},'FontSize',10,'unit','normalized','horizontalalignment','left',...
%     'EdgeColor','k','Margin',2);
% axis(ax,[0 96 -26 18]);
% colormap(ax,'jet');
% colorbar(ax);
% caxis(ax,[0.35 0.55]);
% hold(ax,'off');

ax = f8.ax(1);
hold(ax,'on');
scatter(ax,hfall(ih04,26),hfall(ih04,2), msizehf, 'filled','ko');  %, 'MarkerEdgeColor', 'w')
scatter(ax,lfall(il04,26),lfall(il04,2), msizehf, 'filled','ro');  %, 'MarkerEdgeColor', 'w'))
text(ax,0.05,0.93,'2004-all detections','FontSize',10,'unit','normalized','horizontalalignment','left',...
    'EdgeColor','k','Margin',2);
axis(ax,[0 96 -25 15]);
ax.YLabel.String = 'Along strike (km)';
hold(ax,'off');

ax = f8.ax(2);
hold(ax,'on');
scatter(ax,hfall(ih04,26),hfall(ih04,2),msizehf,log(hfall(ih04,23)),'filled','o');  %, 'MarkerEdgeColor', 'w')
text(ax,0.05,0.93,{'2004-hf energy'},'FontSize',10,'unit','normalized','horizontalalignment','left',...
    'EdgeColor','k','Margin',2);
axis(ax,[0 96 -25 15]);
colormap(ax,'jet');
colorbar(ax);
cran = [prctile(log(hfall(ih04,23)),2), prctile(log(hfall(ih04,23)),98)];
% caxis(ax,cran);
caxis(ax,[-3 2]);
hold(ax,'off');

ax = f8.ax(3);
hold(ax,'on');
scatter(ax,lfall(il04,26),lfall(il04,2),msizehf,log(lfall(il04,23)),'filled','o');  %, 'MarkerEdgeColor', 'w')
text(ax,0.05,0.93,{'2004-lf energy'},'FontSize',10,'unit','normalized','horizontalalignment','left',...
    'EdgeColor','k','Margin',2);
axis(ax,[0 96 -25 15]);
colormap(ax,'jet');
colorbar(ax);
cran = [prctile(log(lfall(il04,23)),2), prctile(log(lfall(il04,23)),98)];
% caxis(ax,cran);
caxis(ax,[-3 2]);
hold(ax,'off');

ax = f8.ax(4);
hold(ax,'on');
scatter(ax,hfaall(iha04,23),hfaall(iha04,2),msizehf,hfaall(iha04,21),'filled','o');  %, 'MarkerEdgeColor', 'w'))
text(ax,0.05,0.93,{'2004-hf KLNB add CC'},'FontSize',10,'unit','normalized','horizontalalignment','left',...
    'EdgeColor','k','Margin',2);
axis(ax,[0 96 -25 15]);
colormap(ax,'jet');
colorbar(ax);
cran = [prctile(hfaall(iha04,21),2), prctile(hfaall(iha04,21),98)];
% caxis(ax,cran);
caxis(ax,[0.35 0.75]);
hold(ax,'off');

ax = f8.ax(5);
hold(ax,'on');
scatter(ax,lfaall(ila04,23),lfaall(ila04,2),msizehf,lfaall(ila04,21),'filled','o');  %, 'MarkerEdgeColor', 'w'))
text(ax,0.05,0.93,{'2004-lf KLNB add CC'},'FontSize',10,'unit','normalized','horizontalalignment','left',...
    'EdgeColor','k','Margin',2);
axis(ax,[0 96 -25 15]);
colormap(ax,'jet');
colorbar(ax);
cran = [prctile(lfaall(ila04,21),2), prctile(lfaall(ila04,21),98)];
% caxis(ax,cran);
caxis(ax,[0.35 0.75]);
hold(ax,'off');

ax = f8.ax(6);
hold(ax,'on');
scatter(ax,hfaall(iha04,23),hfaall(iha04,2), msizehf, 'filled','ko');  %, 'MarkerEdgeColor', 'w')
scatter(ax,lfaall(ila04,23),lfaall(ila04,2), msizehf, 'filled','ro');  %, 'MarkerEdgeColor', 'w'))
text(ax,0.05,0.93,'2004-KLNB checked','FontSize',10,'unit','normalized','horizontalalignment','left',...
    'EdgeColor','k','Margin',2);
axis(ax,[0 96 -25 15]);
ax.XLabel.String = 'Time (hr)';
hold(ax,'off');

% ax = f8.ax(7);
% hold(ax,'on');
% scatter(ax,hfall(ih05,24),hfall(ih05,2), msizehf, 'filled','ko');  %, 'MarkerEdgeColor', 'w')
% scatter(ax,lfall(il05,24),lfall(il05,2), msizehf, 'filled','ro');  %, 'MarkerEdgeColor', 'w'))
% text(ax,0.05,0.93,'2005','FontSize',10,'unit','normalized','horizontalalignment','left',...
%     'EdgeColor','k','Margin',2);
% axis(ax,[0 96 -26 18]);
% hold(ax,'off');
% 
% ax = f8.ax(8);
% hold(ax,'on');
% scatter(ax,hfall(ih05,24),hfall(ih05,2),msizehf,hfall(ih05,21)-hfall(ih05,22),'filled','o');  %, 'MarkerEdgeColor', 'w')
% text(ax,0.05,0.93,{'2005-hf diff amp'},'FontSize',10,'unit','normalized','horizontalalignment','left',...
%     'EdgeColor','k','Margin',2);
% axis(ax,[0 96 -26 18]);
% colormap(ax,'jet');
% colorbar(ax);
% caxis(ax,[0.1 0.7]);
% hold(ax,'off');
% 
% ax = f8.ax(9);
% hold(ax,'on');
% scatter(ax,lfall(il05,24),lfall(il05,2),msizehf,median(lfall(il05,17:20),2),'filled','o');  %, 'MarkerEdgeColor', 'w'))
% text(ax,0.05,0.93,{'2005-lf med add CC'},'FontSize',10,'unit','normalized','horizontalalignment','left',...
%     'EdgeColor','k','Margin',2);
% axis(ax,[0 96 -26 18]);
% colormap(ax,'jet');
% colorbar(ax);
% caxis(ax,[0.35 0.55]);
% ax.XLabel.String = 'Time (hr)';
% hold(ax,'off');

print(f8.fig,'-depsc2',strcat(rstpath,'/nocutallstrike2004.eps'));


%% all catalog, along strike, 2005
f9.fig = figure;
f9.fig.Renderer = 'painters';
widin = 12;  % maximum width allowed is 8.5 inches
htin = 9.5;   % maximum height allowed is 11 inches
set(f9.fig,'Position',[1*scrsz(3)/20 scrsz(4)/10 widin*res htin*res]);
nrow = 6;
ncol = 1;
for isub = 1:nrow*ncol
    f9.ax(isub) = subplot(nrow,ncol,isub);
end

% reposition
for isub = 1:nrow*ncol
    f9.ax(isub).Box = 'on';
    grid(f9.ax(isub), 'on');
    f9.ax(isub).GridLineStyle = '--';
end

% marker size
msizehf = 2;

ax = f9.ax(1);
hold(ax,'on');
scatter(ax,hfall(ih05,26),hfall(ih05,2), msizehf, 'filled','ko');  %, 'MarkerEdgeColor', 'w')
scatter(ax,lfall(il05,26),lfall(il05,2), msizehf, 'filled','ro');  %, 'MarkerEdgeColor', 'w'))
text(ax,0.05,0.93,'2005-all detections','FontSize',10,'unit','normalized','horizontalalignment','left',...
    'EdgeColor','k','Margin',2);
axis(ax,[0 96 -28 20]);
ax.YLabel.String = 'Along strike (km)';
hold(ax,'off');

ax = f9.ax(2);
hold(ax,'on');
scatter(ax,hfall(ih05,26),hfall(ih05,2),msizehf,log(hfall(ih05,23)),'filled','o');  %, 'MarkerEdgeColor', 'w')
text(ax,0.05,0.93,{'2005-hf energy'},'FontSize',10,'unit','normalized','horizontalalignment','left',...
    'EdgeColor','k','Margin',2);
axis(ax,[0 96 -28 20]);
colormap(ax,'jet');
colorbar(ax);
cran = [prctile(log(hfall(ih05,23)),2), prctile(log(hfall(ih05,23)),98)];
% caxis(ax,cran);
caxis(ax,[-3 2]);
hold(ax,'off');

ax = f9.ax(3);
hold(ax,'on');
scatter(ax,lfall(il05,26),lfall(il05,2),msizehf,log(lfall(il05,23)),'filled','o');  %, 'MarkerEdgeColor', 'w')
text(ax,0.05,0.93,{'2005-lf energy'},'FontSize',10,'unit','normalized','horizontalalignment','left',...
    'EdgeColor','k','Margin',2);
axis(ax,[0 96 -28 20]);
colormap(ax,'jet');
colorbar(ax);
cran = [prctile(log(lfall(il05,23)),2), prctile(log(lfall(il05,23)),98)];
% caxis(ax,cran);
caxis(ax,[-3 2]);
hold(ax,'off');

ax = f9.ax(4);
hold(ax,'on');
scatter(ax,hfaall(iha05,23),hfaall(iha05,2),msizehf,hfaall(iha05,21),'filled','o');  %, 'MarkerEdgeColor', 'w'))
text(ax,0.05,0.93,{'2005-hf KLNB add CC'},'FontSize',10,'unit','normalized','horizontalalignment','left',...
    'EdgeColor','k','Margin',2);
axis(ax,[0 96 -28 20]);
colormap(ax,'jet');
colorbar(ax);
cran = [prctile(hfaall(iha05,21),2), prctile(hfaall(iha05,21),98)];
% caxis(ax,cran);
caxis(ax,[0.35 0.75]);
hold(ax,'off');

ax = f9.ax(5);
hold(ax,'on');
scatter(ax,lfaall(ila05,23),lfaall(ila05,2),msizehf,lfaall(ila05,21),'filled','o');  %, 'MarkerEdgeColor', 'w'))
text(ax,0.05,0.93,{'2005-lf KLNB add CC'},'FontSize',10,'unit','normalized','horizontalalignment','left',...
    'EdgeColor','k','Margin',2);
axis(ax,[0 96 -28 20]);
colormap(ax,'jet');
colorbar(ax);
cran = [prctile(lfaall(ila05,21),2), prctile(lfaall(ila05,21),98)];
% caxis(ax,cran);
caxis(ax,[0.35 0.75]);
hold(ax,'off');

ax = f9.ax(6);
hold(ax,'on');
scatter(ax,hfaall(iha05,23),hfaall(iha05,2), msizehf, 'filled','ko');  %, 'MarkerEdgeColor', 'w')
scatter(ax,lfaall(ila05,23),lfaall(ila05,2), msizehf, 'filled','ro');  %, 'MarkerEdgeColor', 'w'))
text(ax,0.05,0.93,'2005-KLNB checked','FontSize',10,'unit','normalized','horizontalalignment','left',...
    'EdgeColor','k','Margin',2);
axis(ax,[0 96 -28 20]);
ax.XLabel.String = 'Time (hr)';
hold(ax,'off');

print(f9.fig,'-depsc2',strcat(rstpath,'/nocutallstrike2005.eps'));


%% region 2, along dip
f7.fig = figure;
f7.fig.Renderer = 'painters';
widin = 12;  % maximum width allowed is 8.5 inches
htin = 8;   % maximum height allowed is 11 inches
set(f7.fig,'Position',[1*scrsz(3)/20 scrsz(4)/10 widin*res htin*res]);
nrow = 4;
ncol = 1;
for isub = 1:nrow*ncol
    f7.ax(isub) = subplot(nrow,ncol,isub);
end

% reposition
for isub = 1:nrow*ncol
    f7.ax(isub).Box = 'on';
    grid(f7.ax(isub), 'on');
    f7.ax(isub).GridLineStyle = '--';
end

% marker size
msizehf = 4;

ax = f7.ax(1);
hold(ax,'on');
scatter(ax,hfall(ih05,26),hfall(ih05,1), msizehf, 'filled','ko');  %, 'MarkerEdgeColor', 'w')
scatter(ax,lfall(il05,26),lfall(il05,1), msizehf, 'filled','ro');  %, 'MarkerEdgeColor', 'w'))
text(ax,0.05,0.93,'2005-detections','FontSize',10,'unit','normalized','horizontalalignment','left',...
    'EdgeColor','k','Margin',2);
axis(ax,[70 100 -20 35]);
ax.YLabel.String = 'Along down-dip (km)';
hold(ax,'off');

ax = f7.ax(2);
hold(ax,'on');
reg3 = [75 -27;
        100 -27;
        100 -13;
        75  -13 ;
        75 -27];
scatter(ax,hfall(ih05,26),hfall(ih05,2), msizehf, 'filled','ko');  %, 'MarkerEdgeColor', 'w')
scatter(ax,lfall(il05,26),lfall(il05,2), msizehf, 'filled','ro');  %, 'MarkerEdgeColor', 'w'))
plot(ax,reg3(:,1),reg3(:,2),'b-','linew',1);
text(ax,0.05,0.93,'2005-detections','FontSize',10,'unit','normalized','horizontalalignment','left',...
    'EdgeColor','k','Margin',2);
axis(ax,[70 100 -28 20]);
ax.YLabel.String = 'Along strike (km)';
hold(ax,'off');

ax = f7.ax(3);
hold(ax,'on');
scatter(ax,hfall(ih05,26),hfall(ih05,1), msizehf, 'filled','ko');  %, 'MarkerEdgeColor', 'w')
scatter(ax,lfall(il05,26),lfall(il05,1), msizehf, 'filled','ro');  %, 'MarkerEdgeColor', 'w'))
text(ax,0.05,0.93,'2005-detections','FontSize',10,'unit','normalized','horizontalalignment','left',...
    'EdgeColor','k','Margin',2);
axis(ax,[160 190 -20 35]);
ax.YLabel.String = 'Along down-dip (km)';
hold(ax,'off');

ax = f7.ax(4);
hold(ax,'on');
scatter(ax,hfall(ih05,26),hfall(ih05,2), msizehf, 'filled','ko');  %, 'MarkerEdgeColor', 'w')
scatter(ax,lfall(il05,26),lfall(il05,2), msizehf, 'filled','ro');  %, 'MarkerEdgeColor', 'w'))
text(ax,0.05,0.93,'2005-detections','FontSize',10,'unit','normalized','horizontalalignment','left',...
    'EdgeColor','k','Margin',2);
axis(ax,[160 190 -28 20]);
ax.YLabel.String = 'Along strike (km)';
ax.XLabel.String = 'Time (hr)';
hold(ax,'off');

print(f9.fig,'-depsc2',strcat(rstpath,'/nocutalldip2005.eps'));


%% check the weird belt in lf in 2005, see the time and space plot
reg3 = [75 -27;
        100 -27;
        100 -13;
        75  -13 ;
        75 -27];

% events & tremor inside bound 1
[is,ion] = inpolygon(lfall(il05,26),lfall(il05,2),reg3(:,1),reg3(:,2));
isinreg3 = is | ion;
ind = il05(isinreg3 == 1);
lfreg3 = lfrelatime(ind, :);

f11.fig = figure;
widin = 8;  % maximum width allowed is 8.5 inches
htin = 10;   % maximum height allowed is 11 inches
set(f11.fig,'Position',[1*scrsz(3)/20 scrsz(4)/10 widin*res htin*res]);
nrow = 2;
ncol = 2;
for isub = 1:nrow*ncol
    f11.ax(isub) = subplot(nrow,ncol,isub);
    f11.ax(isub).Box = 'on';
    grid(f11.ax(isub), 'on');
    f11.ax(isub).GridLineStyle = '--';
    axis(f11.ax(isub), 'equal');
    axis(f11.ax(isub),[-30 50 -30 40]);
    
end

% reposition
set(f11.ax(3),'Position', [0.1, 0.1, 0.4, 0.4]);
set(f11.ax(4),'Position', [0.55, 0.1, 0.4, 0.4]);
set(f11.ax(1),'Position', [0.1, 0.55, 0.4, 0.4]);
set(f11.ax(2),'Position', [0.55, 0.55, 0.4, 0.4]);

% marker size
msizelf = 6;


% subplot 2 of figure 1
hold(f11.ax(2),'on');
dumlf = lfreg3;
scatter(f11.ax(2),dumlf(:,1),dumlf(:,2), msizelf,'r','filled','o');  %, 'MarkerEdgeColor', 'w')
scatter(f11.ax(2),relacont(:,1),relacont(:,2),20,'k','o','LineWidth',1);
plot(f11.ax(2),[-100 100],[0 0],'k--');
plot(f11.ax(2),[0 0],[-100 100],'k--');
text(f11.ax(2),0.05,0.93,'b','FontSize',12,'unit','normalized','horizontalalignment','center',...
    'EdgeColor','k','Margin',2);
text(f11.ax(2),0.93,0.93,'LZB','FontSize',12,'unit','normalized','horizontalalignment','center',...
    'EdgeColor','k','Margin',2);
text(f11.ax(2),0.93,0.1,'LF','FontSize',15,'unit','normalized','horizontalalignment','center');

hold(f11.ax(2),'off');

print(f11.fig,'-depsc2',strcat(rstpath,'/weirdbelt2005.eps'));






