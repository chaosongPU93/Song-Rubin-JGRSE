% function plot_all_rela_LZB
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
% First created date:   2019/11/21
% Last modified date:   2019/11/21
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

print(f.fig,'-dpdf',strcat(rstpath,'/asd.rela.pdf'));

% %%% save figure to file
% print(f.fig,'-dpdf',strcat(rstpath,'/allfam.hflf.alld.nodcut.accuhit.',num2str(winlenhf),'_',num2str(winlenlf),...
%       '.',num2str(lofflf),'.',num2str(ccminlf),'.rela.pdf'));
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


%% region 1, along down-dip
f3.fig = figure;
f3.fig.Renderer = 'painters';
widin = 8;  % maximum width allowed is 8.5 inches
htin = 9;   % maximum height allowed is 11 inches
set(f3.fig,'Position',[1*scrsz(3)/20 scrsz(4)/10 widin*res htin*res]);
nrow = 10;
ncol = 1;
for isub = 1:nrow*ncol
    f3.ax(isub) = subplot(nrow,ncol,isub);
end

% reposition
for isub = 1:nrow*ncol
    %     if mod(isub,ncol)==0
    %         icol = ncol;
    %         irow = floor(isub/ncol);
    %     else
    %         icol = mod(isub,ncol);
    %         irow = floor(isub/ncol)+1;
    %     end
    %     set(f3.ax(isub),'Position', [(icol-0.85)/ncol, 1-irow/nrow*0.9, 1/ncol*0.8 1/nrow*0.8]);
    f3.ax(isub).Box = 'on';
    grid(f3.ax(isub), 'on');
    f3.ax(isub).GridLineStyle = '--';
end

% marker size
msizehf = 2;

ax = f3.ax(1);
hold(ax,'on');
scatter(ax,hfreg1(ih103,24),hfreg1(ih103,1), msizehf, 'filled','o');  %, 'MarkerEdgeColor', 'w')
scatter(ax,lfreg1(il103,24),lfreg1(il103,1), msizehf, 'filled','ro');  %, 'MarkerEdgeColor', 'w'))
text(ax,0.05,0.93,'2003','FontSize',13,'unit','normalized','horizontalalignment','center',...
    'EdgeColor','k','Margin',2);
axis(ax,[0 60 -18 0]);
% axis(ax,[0 60 5 25]);
ax.YLabel.String = 'Along down-dip (km)';
% ax.YLabel.String = 'Along strike (km)';
hold(ax,'off');

ax = f3.ax(2);
hold(ax,'on');
scatter(ax,hfreg1(ih103,24),hfreg1(ih103,1), msizehf, 'filled','o');  %, 'MarkerEdgeColor', 'w')
scatter(ax,lfreg1(il103,24),lfreg1(il103,1), msizehf, 'filled','ro');  %, 'MarkerEdgeColor', 'w'))
axis(ax,[60 120 -18 0]);
% axis(ax,[60 120 5 25]);
% ax.YLabel.String = 'Along down-dip (km)';
% ax.YLabel.String = 'Along strike (km)';
hold(ax,'off');

ax = f3.ax(3);
hold(ax,'on');
scatter(ax,hfreg1(ih104,24),hfreg1(ih104,1), msizehf, 'filled','o');  %, 'MarkerEdgeColor', 'w')
scatter(ax,lfreg1(il104,24),lfreg1(il104,1), msizehf, 'filled','ro');  %, 'MarkerEdgeColor', 'w'))
text(ax,0.05,0.93,'2004','FontSize',13,'unit','normalized','horizontalalignment','center',...
    'EdgeColor','k','Margin',2);
axis(ax,[0 60 -18 0]);
% axis(ax,[0 60 5 25]);
% ax.YLabel.String = 'Along down-dip (km)';
% ax.YLabel.String = 'Along strike (km)';
hold(ax,'off');

ax = f3.ax(4);
hold(ax,'on');
scatter(ax,hfreg1(ih104,24),hfreg1(ih104,1), msizehf, 'filled','o');  %, 'MarkerEdgeColor', 'w')
scatter(ax,lfreg1(il104,24),lfreg1(il104,1), msizehf, 'filled','ro');  %, 'MarkerEdgeColor', 'w'))
axis(ax,[60 120 -18 0]);
% axis(ax,[60 120 5 25]);
% ax.YLabel.String = 'Along down-dip (km)';
% ax.YLabel.String = 'Along strike (km)';
hold(ax,'off');

ax = f3.ax(5);
hold(ax,'on');
scatter(ax,hfreg1(ih104,24),hfreg1(ih104,1), msizehf, 'filled','o');  %, 'MarkerEdgeColor', 'w')
scatter(ax,lfreg1(il104,24),lfreg1(il104,1), msizehf, 'filled','ro');  %, 'MarkerEdgeColor', 'w'))
axis(ax,[120 180 -18 0]);
% axis(ax,[120 180 5 25]);
% ax.YLabel.String = 'Along down-dip (km)';
% ax.YLabel.String = 'Along strike (km)';
hold(ax,'off');

ax = f3.ax(6);
hold(ax,'on');
scatter(ax,hfreg1(ih104,24),hfreg1(ih104,1), msizehf, 'filled','o');  %, 'MarkerEdgeColor', 'w')
scatter(ax,lfreg1(il104,24),lfreg1(il104,1), msizehf, 'filled','ro');  %, 'MarkerEdgeColor', 'w'))
axis(ax,[180 240 -18 0]);
% axis(ax,[180 240 5 25]);
% ax.YLabel.String = 'Along down-dip (km)';
% ax.YLabel.String = 'Along strike (km)';
hold(ax,'off');

ax = f3.ax(7);
hold(ax,'on');
scatter(ax,hfreg1(ih105,24),hfreg1(ih105,1), msizehf, 'filled','o');  %, 'MarkerEdgeColor', 'w')
scatter(ax,lfreg1(il105,24),lfreg1(il105,1), msizehf, 'filled','ro');  %, 'MarkerEdgeColor', 'w'))
text(ax,0.05,0.93,'2005','FontSize',13,'unit','normalized','horizontalalignment','center',...
    'EdgeColor','k','Margin',2);
axis(ax,[0 60 -18 0]);
% axis(ax,[0 60 5 25]);
% ax.YLabel.String = 'Along down-dip (km)';
% ax.YLabel.String = 'Along strike (km)';
hold(ax,'off');

ax = f3.ax(8);
hold(ax,'on');
scatter(ax,hfreg1(ih105,24),hfreg1(ih105,1), msizehf, 'filled','o');  %, 'MarkerEdgeColor', 'w')
scatter(ax,lfreg1(il105,24),lfreg1(il105,1), msizehf, 'filled','ro');  %, 'MarkerEdgeColor', 'w'))
axis(ax,[60 120 -18 0]);
% axis(ax,[60 120 5 25]);
% ax.YLabel.String = 'Along down-dip (km)';
% ax.YLabel.String = 'Along strike (km)';
hold(ax,'off');

ax = f3.ax(9);
hold(ax,'on');
scatter(ax,hfreg1(ih105,24),hfreg1(ih105,1), msizehf, 'filled','o');  %, 'MarkerEdgeColor', 'w')
scatter(ax,lfreg1(il105,24),lfreg1(il105,1), msizehf, 'filled','ro');  %, 'MarkerEdgeColor', 'w'))
axis(ax,[120 180 -18 0]);
% axis(ax,[120 180 5 25]);
% ax.YLabel.String = 'Along down-dip (km)';
% ax.YLabel.String = 'Along strike (km)';
hold(ax,'off');

ax = f3.ax(10);
hold(ax,'on');
scatter(ax,hfreg1(ih105,24),hfreg1(ih105,1), msizehf, 'filled','o');  %, 'MarkerEdgeColor', 'w')
scatter(ax,lfreg1(il105,24),lfreg1(il105,1), msizehf, 'filled','ro');  %, 'MarkerEdgeColor', 'w'))
axis(ax,[180 240 -18 0]);
% axis(ax,[180 240 5 25]);
% ax.YLabel.String = 'Along down-dip (km)';
% ax.YLabel.String = 'Along strike (km)';
ax.XLabel.String = 'Time (hr)';
hold(ax,'off');

print(f3.fig,'-dpdf',strcat(rstpath,'/reg1downdip.pdf'));


%% region 1, along strike
f4.fig = figure;
f4.fig.Renderer = 'painters';
widin = 8;  % maximum width allowed is 8.5 inches
htin = 9;   % maximum height allowed is 11 inches
set(f4.fig,'Position',[1*scrsz(3)/20 scrsz(4)/10 widin*res htin*res]);
nrow = 10;
ncol = 1;
for isub = 1:nrow*ncol
    f4.ax(isub) = subplot(nrow,ncol,isub);
end

% reposition
for isub = 1:nrow*ncol
    %     if mod(isub,ncol)==0
    %         icol = ncol;
    %         irow = floor(isub/ncol);
    %     else
    %         icol = mod(isub,ncol);
    %         irow = floor(isub/ncol)+1;
    %     end
    %     set(f4.ax(isub),'Position', [(icol-0.85)/ncol, 1-irow/nrow*0.9, 1/ncol*0.8 1/nrow*0.8]);
    f4.ax(isub).Box = 'on';
    grid(f4.ax(isub), 'on');
    f4.ax(isub).GridLineStyle = '--';
end

% marker size
msizehf = 2;

ax = f4.ax(1);
hold(ax,'on');
scatter(ax,hfreg1(ih103,26),hfreg1(ih103,2), msizehf, 'filled','o');  %, 'MarkerEdgeColor', 'w')
scatter(ax,lfreg1(il103,26),lfreg1(il103,2), msizehf, 'filled','ro');  %, 'MarkerEdgeColor', 'w'))
text(ax,0.05,0.93,'2003','FontSize',13,'unit','normalized','horizontalalignment','center',...
    'EdgeColor','k','Margin',2);
axis(ax,[0 60 -15 5]);
% ax.YLabel.String = 'Along down-dip (km)';
ax.YLabel.String = 'Along strike (km)';
hold(ax,'off');

ax = f4.ax(2);
hold(ax,'on');
scatter(ax,hfreg1(ih103,26),hfreg1(ih103,2), msizehf, 'filled','o');  %, 'MarkerEdgeColor', 'w')
scatter(ax,lfreg1(il103,26),lfreg1(il103,2), msizehf, 'filled','ro');  %, 'MarkerEdgeColor', 'w'))
axis(ax,[60 120 -15 5]);
% ax.YLabel.String = 'Along down-dip (km)';
% ax.YLabel.String = 'Along strike (km)';
hold(ax,'off');

ax = f4.ax(3);
hold(ax,'on');
scatter(ax,hfreg1(ih104,26),hfreg1(ih104,2), msizehf, 'filled','o');  %, 'MarkerEdgeColor', 'w')
scatter(ax,lfreg1(il104,26),lfreg1(il104,2), msizehf, 'filled','ro');  %, 'MarkerEdgeColor', 'w'))
text(ax,0.05,0.93,'2004','FontSize',13,'unit','normalized','horizontalalignment','center',...
    'EdgeColor','k','Margin',2);
axis(ax,[0 60 -15 5]);
% ax.YLabel.String = 'Along down-dip (km)';
% ax.YLabel.String = 'Along strike (km)';
hold(ax,'off');

ax = f4.ax(4);
hold(ax,'on');
scatter(ax,hfreg1(ih104,26),hfreg1(ih104,2), msizehf, 'filled','o');  %, 'MarkerEdgeColor', 'w')
scatter(ax,lfreg1(il104,26),lfreg1(il104,2), msizehf, 'filled','ro');  %, 'MarkerEdgeColor', 'w'))
axis(ax,[60 120 -15 5]);
% ax.YLabel.String = 'Along down-dip (km)';
% ax.YLabel.String = 'Along strike (km)';
hold(ax,'off');

ax = f4.ax(5);
hold(ax,'on');
scatter(ax,hfreg1(ih104,26),hfreg1(ih104,2), msizehf, 'filled','o');  %, 'MarkerEdgeColor', 'w')
scatter(ax,lfreg1(il104,26),lfreg1(il104,2), msizehf, 'filled','ro');  %, 'MarkerEdgeColor', 'w'))
axis(ax,[120 180 -15 5]);
% ax.YLabel.String = 'Along down-dip (km)';
% ax.YLabel.String = 'Along strike (km)';
hold(ax,'off');

ax = f4.ax(6);
hold(ax,'on');
scatter(ax,hfreg1(ih104,26),hfreg1(ih104,2), msizehf, 'filled','o');  %, 'MarkerEdgeColor', 'w')
scatter(ax,lfreg1(il104,26),lfreg1(il104,2), msizehf, 'filled','ro');  %, 'MarkerEdgeColor', 'w'))
axis(ax,[180 240 -15 5]);
% ax.YLabel.String = 'Along down-dip (km)';
% ax.YLabel.String = 'Along strike (km)';
hold(ax,'off');

ax = f4.ax(7);
hold(ax,'on');
scatter(ax,hfreg1(ih105,26),hfreg1(ih105,2), msizehf, 'filled','o');  %, 'MarkerEdgeColor', 'w')
scatter(ax,lfreg1(il105,26),lfreg1(il105,2), msizehf, 'filled','ro');  %, 'MarkerEdgeColor', 'w'))
text(ax,0.05,0.93,'2005','FontSize',13,'unit','normalized','horizontalalignment','center',...
    'EdgeColor','k','Margin',2);
axis(ax,[0 60 -15 5]);
% ax.YLabel.String = 'Along down-dip (km)';
% ax.YLabel.String = 'Along strike (km)';
hold(ax,'off');

ax = f4.ax(8);
hold(ax,'on');
scatter(ax,hfreg1(ih105,26),hfreg1(ih105,2), msizehf, 'filled','o');  %, 'MarkerEdgeColor', 'w')
scatter(ax,lfreg1(il105,26),lfreg1(il105,2), msizehf, 'filled','ro');  %, 'MarkerEdgeColor', 'w'))
axis(ax,[60 120 -15 5]);
% ax.YLabel.String = 'Along down-dip (km)';
% ax.YLabel.String = 'Along strike (km)';
hold(ax,'off');

ax = f4.ax(9);
hold(ax,'on');
scatter(ax,hfreg1(ih105,26),hfreg1(ih105,2), msizehf, 'filled','o');  %, 'MarkerEdgeColor', 'w')
scatter(ax,lfreg1(il105,26),lfreg1(il105,2), msizehf, 'filled','ro');  %, 'MarkerEdgeColor', 'w'))
axis(ax,[120 180 -15 5]);
% ax.YLabel.String = 'Along down-dip (km)';
% ax.YLabel.String = 'Along strike (km)';
hold(ax,'off');

ax = f4.ax(10);
hold(ax,'on');
scatter(ax,hfreg1(ih105,26),hfreg1(ih105,2), msizehf, 'filled','o');  %, 'MarkerEdgeColor', 'w')
scatter(ax,lfreg1(il105,26),lfreg1(il105,2), msizehf, 'filled','ro');  %, 'MarkerEdgeColor', 'w'))
axis(ax,[180 240 -15 5]);
% ax.YLabel.String = 'Along down-dip (km)';
% ax.YLabel.String = 'Along strike (km)';
ax.XLabel.String = 'Time (hr)';
hold(ax,'off');

print(f4.fig,'-dpdf',strcat(rstpath,'/reg1strike.pdf'));


%% region 2, along down-dip
f5.fig = figure;
f5.fig.Renderer = 'painters';
widin = 8;  % maximum width allowed is 8.5 inches
htin = 9;   % maximum height allowed is 11 inches
set(f5.fig,'Position',[1*scrsz(3)/20 scrsz(4)/10 widin*res htin*res]);
nrow = 10;
ncol = 1;
for isub = 1:nrow*ncol
    f5.ax(isub) = subplot(nrow,ncol,isub);
end

% reposition
for isub = 1:nrow*ncol
    %     if mod(isub,ncol)==0
    %         icol = ncol;
    %         irow = floor(isub/ncol);
    %     else
    %         icol = mod(isub,ncol);
    %         irow = floor(isub/ncol)+1;
    %     end
    %     set(f5.ax(isub),'Position', [(icol-0.85)/ncol, 1-irow/nrow*0.9, 1/ncol*0.8 1/nrow*0.8]);
    f5.ax(isub).Box = 'on';
    grid(f5.ax(isub), 'on');
    f5.ax(isub).GridLineStyle = '--';
end

% marker size
msizehf = 2;

ax = f5.ax(1);
hold(ax,'on');
scatter(ax,hfreg2(ih203,24),hfreg2(ih203,1), msizehf, 'filled','o');  %, 'MarkerEdgeColor', 'w')
scatter(ax,lfreg2(il203,24),lfreg2(il203,1), msizehf, 'filled','ro');  %, 'MarkerEdgeColor', 'w'))
text(ax,0.05,0.93,'2003','FontSize',13,'unit','normalized','horizontalalignment','center',...
    'EdgeColor','k','Margin',2);
axis(ax,[0 60 -4 22]);
% axis(ax,[0 60 5 25]);
ax.YLabel.String = 'Along down-dip (km)';
% ax.YLabel.String = 'Along strike (km)';
hold(ax,'off');

ax = f5.ax(2);
hold(ax,'on');
scatter(ax,hfreg2(ih203,24),hfreg2(ih203,1), msizehf, 'filled','o');  %, 'MarkerEdgeColor', 'w')
scatter(ax,lfreg2(il203,24),lfreg2(il203,1), msizehf, 'filled','ro');  %, 'MarkerEdgeColor', 'w'))
axis(ax,[60 120 -4 22]);
% axis(ax,[60 120 5 25]);
% ax.YLabel.String = 'Along down-dip (km)';
% ax.YLabel.String = 'Along strike (km)';
hold(ax,'off');

ax = f5.ax(3);
hold(ax,'on');
scatter(ax,hfreg2(ih204,24),hfreg2(ih204,1), msizehf, 'filled','o');  %, 'MarkerEdgeColor', 'w')
scatter(ax,lfreg2(il204,24),lfreg2(il204,1), msizehf, 'filled','ro');  %, 'MarkerEdgeColor', 'w'))
text(ax,0.05,0.93,'2004','FontSize',13,'unit','normalized','horizontalalignment','center',...
    'EdgeColor','k','Margin',2);
axis(ax,[0 60 -4 22]);
% axis(ax,[0 60 5 25]);
% ax.YLabel.String = 'Along down-dip (km)';
% ax.YLabel.String = 'Along strike (km)';
hold(ax,'off');

ax = f5.ax(4);
hold(ax,'on');
scatter(ax,hfreg2(ih204,24),hfreg2(ih204,1), msizehf, 'filled','o');  %, 'MarkerEdgeColor', 'w')
scatter(ax,lfreg2(il204,24),lfreg2(il204,1), msizehf, 'filled','ro');  %, 'MarkerEdgeColor', 'w'))
axis(ax,[60 120 -4 22]);
% axis(ax,[60 120 5 25]);
% ax.YLabel.String = 'Along down-dip (km)';
% ax.YLabel.String = 'Along strike (km)';
hold(ax,'off');

ax = f5.ax(5);
hold(ax,'on');
scatter(ax,hfreg2(ih204,24),hfreg2(ih204,1), msizehf, 'filled','o');  %, 'MarkerEdgeColor', 'w')
scatter(ax,lfreg2(il204,24),lfreg2(il204,1), msizehf, 'filled','ro');  %, 'MarkerEdgeColor', 'w'))
axis(ax,[120 180 -4 22]);
% axis(ax,[120 180 5 25]);
% ax.YLabel.String = 'Along down-dip (km)';
% ax.YLabel.String = 'Along strike (km)';
hold(ax,'off');

ax = f5.ax(6);
hold(ax,'on');
scatter(ax,hfreg2(ih204,24),hfreg2(ih204,1), msizehf, 'filled','o');  %, 'MarkerEdgeColor', 'w')
scatter(ax,lfreg2(il204,24),lfreg2(il204,1), msizehf, 'filled','ro');  %, 'MarkerEdgeColor', 'w'))
axis(ax,[180 240 -4 22]);
% axis(ax,[180 240 5 25]);
% ax.YLabel.String = 'Along down-dip (km)';
% ax.YLabel.String = 'Along strike (km)';
hold(ax,'off');

ax = f5.ax(7);
hold(ax,'on');
scatter(ax,hfreg2(ih205,24),hfreg2(ih205,1), msizehf, 'filled','o');  %, 'MarkerEdgeColor', 'w')
scatter(ax,lfreg2(il205,24),lfreg2(il205,1), msizehf, 'filled','ro');  %, 'MarkerEdgeColor', 'w'))
text(ax,0.05,0.93,'2005','FontSize',13,'unit','normalized','horizontalalignment','center',...
    'EdgeColor','k','Margin',2);
axis(ax,[0 60 -4 22]);
% axis(ax,[0 60 5 25]);
% ax.YLabel.String = 'Along down-dip (km)';
% ax.YLabel.String = 'Along strike (km)';
hold(ax,'off');

ax = f5.ax(8);
hold(ax,'on');
scatter(ax,hfreg2(ih205,24),hfreg2(ih205,1), msizehf, 'filled','o');  %, 'MarkerEdgeColor', 'w')
scatter(ax,lfreg2(il205,24),lfreg2(il205,1), msizehf, 'filled','ro');  %, 'MarkerEdgeColor', 'w'))
axis(ax,[60 120 -4 22]);
% axis(ax,[60 120 5 25]);
% ax.YLabel.String = 'Along down-dip (km)';
% ax.YLabel.String = 'Along strike (km)';
hold(ax,'off');

ax = f5.ax(9);
hold(ax,'on');
scatter(ax,hfreg2(ih205,24),hfreg2(ih205,1), msizehf, 'filled','o');  %, 'MarkerEdgeColor', 'w')
scatter(ax,lfreg2(il205,24),lfreg2(il205,1), msizehf, 'filled','ro');  %, 'MarkerEdgeColor', 'w'))
axis(ax,[120 180 -4 22]);
% axis(ax,[120 180 5 25]);
% ax.YLabel.String = 'Along down-dip (km)';
% ax.YLabel.String = 'Along strike (km)';
hold(ax,'off');

ax = f5.ax(10);
hold(ax,'on');
scatter(ax,hfreg2(ih205,24),hfreg2(ih205,1), msizehf, 'filled','o');  %, 'MarkerEdgeColor', 'w')
scatter(ax,lfreg2(il205,24),lfreg2(il205,1), msizehf, 'filled','ro');  %, 'MarkerEdgeColor', 'w'))
axis(ax,[180 240 -4 22]);
% axis(ax,[180 240 5 25]);
% ax.YLabel.String = 'Along down-dip (km)';
% ax.YLabel.String = 'Along strike (km)';
ax.XLabel.String = 'Time (hr)';
hold(ax,'off');

print(f5.fig,'-dpdf',strcat(rstpath,'/reg2downdip.pdf'));


%% region 2, along strike
f6.fig = figure;
f6.fig.Renderer = 'painters';
widin = 8;  % maximum width allowed is 8.5 inches
htin = 9;   % maximum height allowed is 11 inches
set(f6.fig,'Position',[1*scrsz(3)/20 scrsz(4)/10 widin*res htin*res]);
nrow = 10;
ncol = 1;
for isub = 1:nrow*ncol
    f6.ax(isub) = subplot(nrow,ncol,isub);
end

% reposition
for isub = 1:nrow*ncol
    %     if mod(isub,ncol)==0
    %         icol = ncol;
    %         irow = floor(isub/ncol);
    %     else
    %         icol = mod(isub,ncol);
    %         irow = floor(isub/ncol)+1;
    %     end
    %     set(f6.ax(isub),'Position', [(icol-0.85)/ncol, 1-irow/nrow*0.9, 1/ncol*0.8 1/nrow*0.8]);
    f6.ax(isub).Box = 'on';
    grid(f6.ax(isub), 'on');
    f6.ax(isub).GridLineStyle = '--';
end

% marker size
msizehf = 2;

ax = f6.ax(1);
hold(ax,'on');
scatter(ax,hfreg2(ih203,24),hfreg2(ih203,2), msizehf, 'filled','o');  %, 'MarkerEdgeColor', 'w')
scatter(ax,lfreg2(il203,24),lfreg2(il203,2), msizehf, 'filled','ro');  %, 'MarkerEdgeColor', 'w'))
text(ax,0.05,0.93,'2003','FontSize',13,'unit','normalized','horizontalalignment','center',...
    'EdgeColor','k','Margin',2);
axis(ax,[0 60 -29 -2]);
% ax.YLabel.String = 'Along down-dip (km)';
ax.YLabel.String = 'Along strike (km)';
hold(ax,'off');

ax = f6.ax(2);
hold(ax,'on');
scatter(ax,hfreg2(ih203,24),hfreg2(ih203,2), msizehf, 'filled','o');  %, 'MarkerEdgeColor', 'w')
scatter(ax,lfreg2(il203,24),lfreg2(il203,2), msizehf, 'filled','ro');  %, 'MarkerEdgeColor', 'w'))
axis(ax,[60 120 -29 -2]);
% ax.YLabel.String = 'Along down-dip (km)';
% ax.YLabel.String = 'Along strike (km)';
hold(ax,'off');

ax = f6.ax(3);
hold(ax,'on');
scatter(ax,hfreg2(ih204,24),hfreg2(ih204,2), msizehf, 'filled','o');  %, 'MarkerEdgeColor', 'w')
scatter(ax,lfreg2(il204,24),lfreg2(il204,2), msizehf, 'filled','ro');  %, 'MarkerEdgeColor', 'w'))
text(ax,0.05,0.93,'2004','FontSize',13,'unit','normalized','horizontalalignment','center',...
    'EdgeColor','k','Margin',2);
axis(ax,[0 60 -29 -2]);
% ax.YLabel.String = 'Along down-dip (km)';
% ax.YLabel.String = 'Along strike (km)';
hold(ax,'off');

ax = f6.ax(4);
hold(ax,'on');
scatter(ax,hfreg2(ih204,24),hfreg2(ih204,2), msizehf, 'filled','o');  %, 'MarkerEdgeColor', 'w')
scatter(ax,lfreg2(il204,24),lfreg2(il204,2), msizehf, 'filled','ro');  %, 'MarkerEdgeColor', 'w'))
axis(ax,[60 120 -29 -2]);
% ax.YLabel.String = 'Along down-dip (km)';
% ax.YLabel.String = 'Along strike (km)';
hold(ax,'off');

ax = f6.ax(5);
hold(ax,'on');
scatter(ax,hfreg2(ih204,24),hfreg2(ih204,2), msizehf, 'filled','o');  %, 'MarkerEdgeColor', 'w')
scatter(ax,lfreg2(il204,24),lfreg2(il204,2), msizehf, 'filled','ro');  %, 'MarkerEdgeColor', 'w'))
axis(ax,[120 180 -29 -2]);
% ax.YLabel.String = 'Along down-dip (km)';
% ax.YLabel.String = 'Along strike (km)';
hold(ax,'off');

ax = f6.ax(6);
hold(ax,'on');
scatter(ax,hfreg2(ih204,24),hfreg2(ih204,2), msizehf, 'filled','o');  %, 'MarkerEdgeColor', 'w')
scatter(ax,lfreg2(il204,24),lfreg2(il204,2), msizehf, 'filled','ro');  %, 'MarkerEdgeColor', 'w'))
axis(ax,[180 240 -29 -2]);
% ax.YLabel.String = 'Along down-dip (km)';
% ax.YLabel.String = 'Along strike (km)';
hold(ax,'off');

ax = f6.ax(7);
hold(ax,'on');
scatter(ax,hfreg2(ih205,24),hfreg2(ih205,2), msizehf, 'filled','o');  %, 'MarkerEdgeColor', 'w')
scatter(ax,lfreg2(il205,24),lfreg2(il205,2), msizehf, 'filled','ro');  %, 'MarkerEdgeColor', 'w'))
text(ax,0.05,0.93,'2005','FontSize',13,'unit','normalized','horizontalalignment','center',...
    'EdgeColor','k','Margin',2);
axis(ax,[0 60 -29 -2]);
% ax.YLabel.String = 'Along down-dip (km)';
% ax.YLabel.String = 'Along strike (km)';
hold(ax,'off');

ax = f6.ax(8);
hold(ax,'on');
scatter(ax,hfreg2(ih205,24),hfreg2(ih205,2), msizehf, 'filled','o');  %, 'MarkerEdgeColor', 'w')
scatter(ax,lfreg2(il205,24),lfreg2(il205,2), msizehf, 'filled','ro');  %, 'MarkerEdgeColor', 'w'))
axis(ax,[60 120 -29 -2]);
% ax.YLabel.String = 'Along down-dip (km)';
% ax.YLabel.String = 'Along strike (km)';
hold(ax,'off');

ax = f6.ax(9);
hold(ax,'on');
scatter(ax,hfreg2(ih205,24),hfreg2(ih205,2), msizehf, 'filled','o');  %, 'MarkerEdgeColor', 'w')
scatter(ax,lfreg2(il205,24),lfreg2(il205,2), msizehf, 'filled','ro');  %, 'MarkerEdgeColor', 'w'))
axis(ax,[120 180 -29 -2]);
% ax.YLabel.String = 'Along down-dip (km)';
% ax.YLabel.String = 'Along strike (km)';
hold(ax,'off');

ax = f6.ax(10);
hold(ax,'on');
scatter(ax,hfreg2(ih205,24),hfreg2(ih205,2), msizehf, 'filled','o');  %, 'MarkerEdgeColor', 'w')
scatter(ax,lfreg2(il205,24),lfreg2(il205,2), msizehf, 'filled','ro');  %, 'MarkerEdgeColor', 'w'))
axis(ax,[180 240 -29 -2]);
% ax.YLabel.String = 'Along down-dip (km)';
% ax.YLabel.String = 'Along strike (km)';
ax.XLabel.String = 'Time (hr)';
hold(ax,'off');

print(f6.fig,'-dpdf',strcat(rstpath,'/reg2strike.pdf'));


%% try the entire catalog
hfall = hfrelatime;
lfall = lfrelatime;
% now the 1,2 col is changing from E(043) and N(043), to down-dip and strike
[hfall(:,1),hfall(:,2)] = coordinate_rot(hfall(:,1),hfall(:,2),rotang,0,0);
[lfall(:,1),lfall(:,2)] = coordinate_rot(lfall(:,1),lfall(:,2),rotang,0,0);

% sort them according to their occurring time
hfall = sortrows(hfall,[13,15]);
lfall = sortrows(lfall,[13,15]);

% obtain the relative time in hr
ih03 = find(hfall(:,13) < 2004*1000);
hfall(ih03,24) = (hfall(ih03,13)-2003060).*24+hfall(ih03,15)./3600;
ih04 = find(hfall(:,13) < 2005*1000 & hfall(:,13) > 2004*1000);
hfall(ih04,24) = (hfall(ih04,13)-2004194).*24+hfall(ih04,15)./3600;
ih05 = find(hfall(:,13) > 2005*1000);
hfall(ih05,24) = (hfall(ih05,13)-2005254).*24+hfall(ih05,15)./3600;

il03 = find(lfall(:,13) < 2004*1000);
lfall(il03,24) = (lfall(il03,13)-2003060).*24+lfall(il03,15)./3600;
il04 = find(lfall(:,13) < 2005*1000 & lfall(:,13) > 2004*1000);
lfall(il04,24) = (lfall(il04,13)-2004194).*24+lfall(il04,15)./3600;
il05 = find(lfall(:,13) > 2005*1000);
lfall(il05,24) = (lfall(il05,13)-2005254).*24+lfall(il05,15)./3600;


%% all catalog, along down-dip
f7.fig = figure;
f7.fig.Renderer = 'painters';
widin = 8;  % maximum width allowed is 8.5 inches
htin = 16;   % maximum height allowed is 11 inches
set(f7.fig,'Position',[1*scrsz(3)/20 scrsz(4)/10 widin*res htin*res]);
nrow = 10;
ncol = 1;
for isub = 1:nrow*ncol
    f7.ax(isub) = subplot(nrow,ncol,isub);
end

% reposition
for isub = 1:nrow*ncol
    %     if mod(isub,ncol)==0
    %         icol = ncol;
    %         irow = floor(isub/ncol);
    %     else
    %         icol = mod(isub,ncol);
    %         irow = floor(isub/ncol)+1;
    %     end
    %     set(f7.ax(isub),'Position', [(icol-0.85)/ncol, 1-irow/nrow*0.9, 1/ncol*0.8 1/nrow*0.8]);
    f7.ax(isub).Box = 'on';
    grid(f7.ax(isub), 'on');
    f7.ax(isub).GridLineStyle = '--';
end

% marker size
msizehf = 2;

ax = f7.ax(1);
hold(ax,'on');
scatter(ax,hfall(ih03,24),hfall(ih03,1), msizehf, 'filled','o');  %, 'MarkerEdgeColor', 'w')
scatter(ax,lfall(il03,24),lfall(il03,1), msizehf, 'filled','ro');  %, 'MarkerEdgeColor', 'w'))
text(ax,0.05,0.93,'2003','FontSize',13,'unit','normalized','horizontalalignment','center',...
    'EdgeColor','k','Margin',2);
axis(ax,[0 60 -18 24]);
% axis(ax,[0 60 5 25]);
ax.YLabel.String = 'Along down-dip (km)';
% ax.YLabel.String = 'Along strike (km)';
hold(ax,'off');

ax = f7.ax(2);
hold(ax,'on');
scatter(ax,hfall(ih03,24),hfall(ih03,1), msizehf, 'filled','o');  %, 'MarkerEdgeColor', 'w')
scatter(ax,lfall(il03,24),lfall(il03,1), msizehf, 'filled','ro');  %, 'MarkerEdgeColor', 'w'))
axis(ax,[60 120 -18 24]);
% axis(ax,[60 120 5 25]);
% ax.YLabel.String = 'Along down-dip (km)';
% ax.YLabel.String = 'Along strike (km)';
hold(ax,'off');

ax = f7.ax(3);
hold(ax,'on');
scatter(ax,hfall(ih04,24),hfall(ih04,1), msizehf, 'filled','o');  %, 'MarkerEdgeColor', 'w')
scatter(ax,lfall(il04,24),lfall(il04,1), msizehf, 'filled','ro');  %, 'MarkerEdgeColor', 'w'))
text(ax,0.05,0.93,'2004','FontSize',13,'unit','normalized','horizontalalignment','center',...
    'EdgeColor','k','Margin',2);
axis(ax,[0 60 -18 24]);
% axis(ax,[0 60 5 25]);
% ax.YLabel.String = 'Along down-dip (km)';
% ax.YLabel.String = 'Along strike (km)';
hold(ax,'off');

ax = f7.ax(4);
hold(ax,'on');
scatter(ax,hfall(ih04,24),hfall(ih04,1), msizehf, 'filled','o');  %, 'MarkerEdgeColor', 'w')
scatter(ax,lfall(il04,24),lfall(il04,1), msizehf, 'filled','ro');  %, 'MarkerEdgeColor', 'w'))
axis(ax,[60 120 -18 24]);
% axis(ax,[60 120 5 25]);
% ax.YLabel.String = 'Along down-dip (km)';
% ax.YLabel.String = 'Along strike (km)';
hold(ax,'off');

ax = f7.ax(5);
hold(ax,'on');
scatter(ax,hfall(ih04,24),hfall(ih04,1), msizehf, 'filled','o');  %, 'MarkerEdgeColor', 'w')
scatter(ax,lfall(il04,24),lfall(il04,1), msizehf, 'filled','ro');  %, 'MarkerEdgeColor', 'w'))
axis(ax,[120 180 -18 24]);
% axis(ax,[120 180 5 25]);
% ax.YLabel.String = 'Along down-dip (km)';
% ax.YLabel.String = 'Along strike (km)';
hold(ax,'off');

ax = f7.ax(6);
hold(ax,'on');
scatter(ax,hfall(ih04,24),hfall(ih04,1), msizehf, 'filled','o');  %, 'MarkerEdgeColor', 'w')
scatter(ax,lfall(il04,24),lfall(il04,1), msizehf, 'filled','ro');  %, 'MarkerEdgeColor', 'w'))
axis(ax,[180 240 -18 24]);
% axis(ax,[180 240 5 25]);
% ax.YLabel.String = 'Along down-dip (km)';
% ax.YLabel.String = 'Along strike (km)';
hold(ax,'off');

ax = f7.ax(7);
hold(ax,'on');
scatter(ax,hfall(ih05,24),hfall(ih05,1), msizehf, 'filled','o');  %, 'MarkerEdgeColor', 'w')
scatter(ax,lfall(il05,24),lfall(il05,1), msizehf, 'filled','ro');  %, 'MarkerEdgeColor', 'w'))
text(ax,0.05,0.93,'2005','FontSize',13,'unit','normalized','horizontalalignment','center',...
    'EdgeColor','k','Margin',2);
axis(ax,[0 60 -18 24]);
% axis(ax,[0 60 5 25]);
% ax.YLabel.String = 'Along down-dip (km)';
% ax.YLabel.String = 'Along strike (km)';
hold(ax,'off');

ax = f7.ax(8);
hold(ax,'on');
scatter(ax,hfall(ih05,24),hfall(ih05,1), msizehf, 'filled','o');  %, 'MarkerEdgeColor', 'w')
scatter(ax,lfall(il05,24),lfall(il05,1), msizehf, 'filled','ro');  %, 'MarkerEdgeColor', 'w'))
axis(ax,[60 120 -18 24]);
% axis(ax,[60 120 5 25]);
% ax.YLabel.String = 'Along down-dip (km)';
% ax.YLabel.String = 'Along strike (km)';
hold(ax,'off');

ax = f7.ax(9);
hold(ax,'on');
scatter(ax,hfall(ih05,24),hfall(ih05,1), msizehf, 'filled','o');  %, 'MarkerEdgeColor', 'w')
scatter(ax,lfall(il05,24),lfall(il05,1), msizehf, 'filled','ro');  %, 'MarkerEdgeColor', 'w'))
axis(ax,[120 180 -18 24]);
% axis(ax,[120 180 5 25]);
% ax.YLabel.String = 'Along down-dip (km)';
% ax.YLabel.String = 'Along strike (km)';
hold(ax,'off');

ax = f7.ax(10);
hold(ax,'on');
scatter(ax,hfall(ih05,24),hfall(ih05,1), msizehf, 'filled','o');  %, 'MarkerEdgeColor', 'w')
scatter(ax,lfall(il05,24),lfall(il05,1), msizehf, 'filled','ro');  %, 'MarkerEdgeColor', 'w'))
axis(ax,[180 240 -18 24]);
% axis(ax,[180 240 5 25]);
% ax.YLabel.String = 'Along down-dip (km)';
% ax.YLabel.String = 'Along strike (km)';
ax.XLabel.String = 'Time (hr)';
hold(ax,'off');

print(f7.fig,'-depsc2',strcat(rstpath,'/alldowndip.eps'));


%% all catalog, along strike
f8.fig = figure;
f8.fig.Renderer = 'painters';
widin = 8;  % maximum width allowed is 8.5 inches
htin = 16;   % maximum height allowed is 11 inches
set(f8.fig,'Position',[1*scrsz(3)/20 scrsz(4)/10 widin*res htin*res]);
nrow = 10;
ncol = 1;
for isub = 1:nrow*ncol
    f8.ax(isub) = subplot(nrow,ncol,isub);
end

% reposition
for isub = 1:nrow*ncol
    %     if mod(isub,ncol)==0
    %         icol = ncol;
    %         irow = floor(isub/ncol);
    %     else
    %         icol = mod(isub,ncol);
    %         irow = floor(isub/ncol)+1;
    %     end
    %     set(f8.ax(isub),'Position', [(icol-0.85)/ncol, 1-irow/nrow*0.9, 1/ncol*0.8 1/nrow*0.8]);
    f8.ax(isub).Box = 'on';
    grid(f8.ax(isub), 'on');
    f8.ax(isub).GridLineStyle = '--';
end

% marker size
msizehf = 2;

ax = f8.ax(1);
hold(ax,'on');
scatter(ax,hfall(ih03,24),hfall(ih03,2), msizehf, 'filled','o');  %, 'MarkerEdgeColor', 'w')
scatter(ax,lfall(il03,24),lfall(il03,2), msizehf, 'filled','ro');  %, 'MarkerEdgeColor', 'w'))
text(ax,0.05,0.93,'2003','FontSize',13,'unit','normalized','horizontalalignment','center',...
    'EdgeColor','k','Margin',2);
axis(ax,[0 60 -36 35]);
% ax.YLabel.String = 'Along down-dip (km)';
ax.YLabel.String = 'Along strike (km)';
hold(ax,'off');

ax = f8.ax(2);
hold(ax,'on');
scatter(ax,hfall(ih03,24),hfall(ih03,2), msizehf, 'filled','o');  %, 'MarkerEdgeColor', 'w')
scatter(ax,lfall(il03,24),lfall(il03,2), msizehf, 'filled','ro');  %, 'MarkerEdgeColor', 'w'))
axis(ax,[60 120 -36 35]);
% ax.YLabel.String = 'Along down-dip (km)';
% ax.YLabel.String = 'Along strike (km)';
hold(ax,'off');

ax = f8.ax(3);
hold(ax,'on');
scatter(ax,hfall(ih04,24),hfall(ih04,2), msizehf, 'filled','o');  %, 'MarkerEdgeColor', 'w')
scatter(ax,lfall(il04,24),lfall(il04,2), msizehf, 'filled','ro');  %, 'MarkerEdgeColor', 'w'))
text(ax,0.05,0.93,'2004','FontSize',13,'unit','normalized','horizontalalignment','center',...
    'EdgeColor','k','Margin',2);
axis(ax,[0 60 -36 35]);
% ax.YLabel.String = 'Along down-dip (km)';
% ax.YLabel.String = 'Along strike (km)';
hold(ax,'off');

ax = f8.ax(4);
hold(ax,'on');
scatter(ax,hfall(ih04,24),hfall(ih04,2), msizehf, 'filled','o');  %, 'MarkerEdgeColor', 'w')
scatter(ax,lfall(il04,24),lfall(il04,2), msizehf, 'filled','ro');  %, 'MarkerEdgeColor', 'w'))
axis(ax,[60 120 -36 35]);
% ax.YLabel.String = 'Along down-dip (km)';
% ax.YLabel.String = 'Along strike (km)';
hold(ax,'off');

ax = f8.ax(5);
hold(ax,'on');
scatter(ax,hfall(ih04,24),hfall(ih04,2), msizehf, 'filled','o');  %, 'MarkerEdgeColor', 'w')
scatter(ax,lfall(il04,24),lfall(il04,2), msizehf, 'filled','ro');  %, 'MarkerEdgeColor', 'w'))
axis(ax,[120 180 -36 35]);
% ax.YLabel.String = 'Along down-dip (km)';
% ax.YLabel.String = 'Along strike (km)';
hold(ax,'off');

ax = f8.ax(6);
hold(ax,'on');
scatter(ax,hfall(ih04,24),hfall(ih04,2), msizehf, 'filled','o');  %, 'MarkerEdgeColor', 'w')
scatter(ax,lfall(il04,24),lfall(il04,2), msizehf, 'filled','ro');  %, 'MarkerEdgeColor', 'w'))
axis(ax,[180 240 -36 35]);
% ax.YLabel.String = 'Along down-dip (km)';
% ax.YLabel.String = 'Along strike (km)';
hold(ax,'off');

ax = f8.ax(7);
hold(ax,'on');
scatter(ax,hfall(ih05,24),hfall(ih05,2), msizehf, 'filled','o');  %, 'MarkerEdgeColor', 'w')
scatter(ax,lfall(il05,24),lfall(il05,2), msizehf, 'filled','ro');  %, 'MarkerEdgeColor', 'w'))
text(ax,0.05,0.93,'2005','FontSize',13,'unit','normalized','horizontalalignment','center',...
    'EdgeColor','k','Margin',2);
axis(ax,[0 60 -36 35]);
% ax.YLabel.String = 'Along down-dip (km)';
% ax.YLabel.String = 'Along strike (km)';
hold(ax,'off');

ax = f8.ax(8);
hold(ax,'on');
scatter(ax,hfall(ih05,24),hfall(ih05,2), msizehf, 'filled','o');  %, 'MarkerEdgeColor', 'w')
scatter(ax,lfall(il05,24),lfall(il05,2), msizehf, 'filled','ro');  %, 'MarkerEdgeColor', 'w'))
axis(ax,[60 120 -36 35]);
% ax.YLabel.String = 'Along down-dip (km)';
% ax.YLabel.String = 'Along strike (km)';
hold(ax,'off');

ax = f8.ax(9);
hold(ax,'on');
scatter(ax,hfall(ih05,24),hfall(ih05,2), msizehf, 'filled','o');  %, 'MarkerEdgeColor', 'w')
scatter(ax,lfall(il05,24),lfall(il05,2), msizehf, 'filled','ro');  %, 'MarkerEdgeColor', 'w'))
axis(ax,[120 180 -36 35]);
% ax.YLabel.String = 'Along down-dip (km)';
% ax.YLabel.String = 'Along strike (km)';
hold(ax,'off');

ax = f8.ax(10);
hold(ax,'on');
scatter(ax,hfall(ih05,24),hfall(ih05,2), msizehf, 'filled','o');  %, 'MarkerEdgeColor', 'w')
scatter(ax,lfall(il05,24),lfall(il05,2), msizehf, 'filled','ro');  %, 'MarkerEdgeColor', 'w'))
axis(ax,[180 240 -36 35]);
% ax.YLabel.String = 'Along down-dip (km)';
% ax.YLabel.String = 'Along strike (km)';
ax.XLabel.String = 'Time (hr)';
hold(ax,'off');

print(f8.fig,'-depsc2',strcat(rstpath,'/allstrike.eps'));



%% plot the time evolution

% SUFFIXhf = strcat('hf.time.',num2str(winlenhf),'_',num2str(winlenlf),'.',num2str(lofflf),...
%                   '.',num2str(ccminlf));
% fname = strcat(rstpath, '/evtloc.allfam.',SUFFIXhf);
% hfmaptime = load(fname);
% % 13 cols, format is:
% %   lon lat dep off12 off13 off12sec off13sec date strongest_arr_time center_of_win cc fam


% select a date to plot the general migration pattern
dateall = unique(hfrelatime(:,13));
refyr = 2005;
refjday = 254;
tmp = jul2dat(refyr,refjday);   % return month, day, year
refday = tmp(2);
daymax = 15;
selhf = hfrelatime(hfrelatime(:,13)>=refyr*1000, :);
mindate = dateall(find(dateall>=refyr*1000, 1,'first'));
timeint = selhf(:,13)-(refyr*1000+refjday)+refday;
timedec = selhf(:,13)/3600/24;
time = timeint+timedec;
selhf(:,24) = time;

sellf = lfrelatime(lfrelatime(:,13)>=refyr*1000, :);
timeint = sellf(:,13)-(refyr*1000+refjday)+refday;
timedec = sellf(:,13)/3600/24;
time = timeint+timedec;
sellf(:,24) = time;

selhfpeng = selhf(selhf(:,24)<=daymax,:);
sellfpeng = sellf(sellf(:,24)<=daymax,:);

%%% define and position the figure frame and axes of each plot
f2.fig = figure;
widin = 8;  % maximum width allowed is 8.5 inches
htin = 9;   % maximum height allowed is 11 inches
set(f2.fig,'Position',[1*scrsz(3)/20 scrsz(4)/10 widin*res htin*res]);
nrow = 2;
ncol = 2;
for isub = 1:nrow*ncol
    f2.ax(isub) = subplot(nrow,ncol,isub);
end

%%% reposition
for isub = 1:nrow*ncol
    if mod(isub,ncol)==0
        icol = ncol;
        irow = floor(isub/ncol);
    else
        icol = mod(isub,ncol);
        irow = floor(isub/ncol)+1;
    end
    set(f2.ax(isub),'Position', [(icol-0.9)/ncol, 1-irow/nrow*0.95, 1/ncol*0.8 1/nrow*0.8]);
end

% marker size
msizehf = 4;
msizelf = 8;

%%% deal with each subplot
% subplot 1 of figure 1
hold(f2.ax(1),'on');
dumhf = sortrows(selhfpeng,-13);
scatter(f2.ax(1),dumhf(:,1),dumhf(:,2), msizehf, dumhf(:,end), 'filled','o');
scatter(f2.ax(1),relacont(:,1),relacont(:,2),20,'k','filled','o');
colormap(f2.ax(1),'jet');
c=colorbar(f2.ax(1),'NorthOutside');
c.Label.String = 'Day in Sept. 2005';
c.Label.FontSize = 13;
caxis(f2.ax(1),[refday daymax]);
plot(f2.ax(1),[-100 100],[0 0],'k--');
plot(f2.ax(1),[0 0],[-100 100],'k--');
text(f2.ax(1),0.05,0.93,'a','FontSize',13,'unit','normalized','horizontalalignment','center',...
    'EdgeColor','k','Margin',2);
text(f2.ax(1),0.93,0.93,'LZB','FontSize',13,'unit','normalized','horizontalalignment','center',...
    'EdgeColor','k','Margin',2);
text(f2.ax(1),0.93,0.1,'HF','FontSize',15,'unit','normalized','horizontalalignment','center');
f2.ax(1).Box = 'on';
grid(f2.ax(1), 'on');
f2.ax(1).GridLineStyle = '--';
axis(f2.ax(1), 'equal');
axis(f2.ax(1),[-30 50 -30 40]);
f2.ax(1).XLabel.String = 'E (km)';
f2.ax(1).YLabel.String = 'N (km)';
hold(f2.ax(1),'off');

% subplot 2 of figure 1
hold(f2.ax(2),'on');
scatter(f2.ax(2),selhfpeng(:,1),selhfpeng(:,2), msizehf, selhfpeng(:,end), 'filled','o');  %, 'MarkerEdgeColor', 'w')
scatter(f2.ax(2),relacont(:,1),relacont(:,2),20,'k','filled','o');
colormap(f2.ax(2),'jet');
c=colorbar(f2.ax(2),'NorthOutside');
c.Label.String = 'Day in Sept. 2005';
c.Label.FontSize = 13;
caxis(f2.ax(2),[refday daymax]);
plot(f2.ax(2),[-100 100],[0 0],'k--');
plot(f2.ax(2),[0 0],[-100 100],'k--');
text(f2.ax(2),0.05,0.93,'b','FontSize',13,'unit','normalized','horizontalalignment','center',...
    'EdgeColor','k','Margin',2);
text(f2.ax(2),0.93,0.93,'LZB','FontSize',13,'unit','normalized','horizontalalignment','center',...
    'EdgeColor','k','Margin',2);
text(f2.ax(2),0.93,0.1,'HF','FontSize',15,'unit','normalized','horizontalalignment','center');
f2.ax(2).Box = 'on';
grid(f2.ax(2), 'on');
f2.ax(2).GridLineStyle = '--';
axis(f2.ax(2), 'equal');
axis(f2.ax(2),[-30 50 -30 40]);
f2.ax(2).XLabel.String = 'E (km)';
f2.ax(2).YLabel.String = 'N (km)';
hold(f2.ax(2),'off');

% subplot 3 of figure 1
hold(f2.ax(3),'on');
dumlf = sortrows(sellfpeng,-13);
scatter(f2.ax(3),dumlf(:,1),dumlf(:,2), msizelf, dumlf(:,end), 'filled','o');  %, 'MarkerEdgeColor', 'w')
scatter(f2.ax(3),relacont(:,1),relacont(:,2),20,'k','filled','o');
colormap(f2.ax(3),'jet');
colorbar(f2.ax(3),'NorthOutside');
caxis(f2.ax(3),[refday daymax]);
plot(f2.ax(3),[-100 100],[0 0],'k--');
plot(f2.ax(3),[0 0],[-100 100],'k--');
text(f2.ax(3),0.05,0.93,'c','FontSize',13,'unit','normalized','horizontalalignment','center',...
    'EdgeColor','k','Margin',2);
text(f2.ax(3),0.93,0.93,'LZB','FontSize',13,'unit','normalized','horizontalalignment','center',...
    'EdgeColor','k','Margin',2);
text(f2.ax(3),0.93,0.1,'LF','FontSize',15,'unit','normalized','horizontalalignment','center');
f2.ax(3).Box = 'on';
grid(f2.ax(3), 'on');
f2.ax(3).GridLineStyle = '--';
axis(f2.ax(3), 'equal');
axis(f2.ax(3),[-30 50 -30 40]);
f2.ax(3).XLabel.String = 'E (km)';
f2.ax(3).YLabel.String = 'N (km)';
hold(f2.ax(3),'off');

% subplot 4 of figure 1
hold(f2.ax(4),'on');
scatter(f2.ax(4),sellfpeng(:,1),sellfpeng(:,2), msizelf, sellfpeng(:,end), 'filled','o');  %, 'MarkerEdgeColor', 'w')
scatter(f2.ax(4),relacont(:,1),relacont(:,2),20,'k','filled','o');
colormap(f2.ax(4),'jet');
colorbar(f2.ax(4),'NorthOutside');
caxis(f2.ax(4),[refday daymax]);
plot(f2.ax(4),[-100 100],[0 0],'k--');
plot(f2.ax(4),[0 0],[-100 100],'k--');
text(f2.ax(4),0.05,0.93,'d','FontSize',13,'unit','normalized','horizontalalignment','center',...
    'EdgeColor','k','Margin',2);
text(f2.ax(4),0.93,0.93,'LZB','FontSize',13,'unit','normalized','horizontalalignment','center',...
    'EdgeColor','k','Margin',2);
text(f2.ax(4),0.93,0.1,'LF','FontSize',15,'unit','normalized','horizontalalignment','center');
f2.ax(4).Box = 'on';
grid(f2.ax(4), 'on');
f2.ax(4).GridLineStyle = '--';
axis(f2.ax(4), 'equal');
axis(f2.ax(4),[-30 50 -30 40]);
f2.ax(4).XLabel.String = 'E (km)';
f2.ax(4).YLabel.String = 'N (km)';
hold(f2.ax(4),'off');

%%% save figure
% shrink(f2.fig,1.25,1.25);
print(f2.fig,'-dpdf',strcat(rstpath,'/allfam.hflf.nodcut.tevo',num2str(refyr),'.',num2str(winlenhf),'_',...
    num2str(winlenlf),'.',num2str(lofflf),'.',num2str(ccminlf),'.rela.pdf'));
% close(f2.fig)














