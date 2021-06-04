% function plot_dcut_rela_LZB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function is used to plot the locations of detections that pass the
% distance cutoff & the duplicates removal from different fams  
%
%   REMEMBER to check the input format !!! 
% 
% Looks like figure 4 in Peng et al. 2015
%
% Modified by Chao Song, chaosong@princeton.edu
% First created date:   2019/11/26
% Last modified date:   2020/06/04
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Initialization
format short e   % Set the format to 5-digit floating point
clear
close all

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

relacont = loccont;
[dx, dy] = absloc2relaloc(loccont(:,1),loccont(:,2),loccont(2,1),loccont(2,2));
relacont(:,1) = dx;
relacont(:,2) = dy;

% load the time results after dist cutoff and duplicates removal
SUFFIXhf = strcat('up.hf.time.',num2str(winlenhf),'_',num2str(winlenlf),'.',num2str(lofflf),...
                  '.',num2str(ccminlf));
fname = strcat(rstpath, '/evtloc.allfam.dcutnodou.',SUFFIXhf);
hfrelatime = load(fname);
% 25 cols, format is:
%   E(043) N(043) E(own) N(own) dist(own) lon lat dep off12 off13 off12sec off13sec date
%   main_arrival_time  cent_of_win  avecc ccadd1 ccadd2 ccadd3 ccadd4 ampmax ampmin eng normeng famnum

SUFFIXlf = strcat('lf.time.',num2str(winlenhf),'_',num2str(winlenlf),'.',num2str(lofflf),...
                  '.',num2str(ccminlf));
fname = strcat(rstpath, '/evtloc.allfam.dcutnodou.',SUFFIXlf);
lfrelatime = load(fname);

% sort the results according to their absolute locations
hftemp = sortrows(hfrelatime, [5, 6]);
lftemp = sortrows(lfrelatime, [5, 6]);


% hf, obtain the counts at same position and combine them from diff fams
[hftempuniq,iuniq,~] = unique(hftemp(:,5:6),'rows','stable');
for i = 1:size(hftempuniq,1)
    [idup,~] = find(hftemp(:,5)==hftempuniq(i,1) & hftemp(:,6)==hftempuniq(i,2));
    ctothf(i,1) = length(idup);   % total count of hf
end
if sum(ctothf) == size(hftemp,1)
   hfrelacount = [hftemp(iuniq,1:2) hftemp(iuniq,5:7) ctothf];
else
    disp('Inconsistency of hf total counts number')
end

% lf, obtain the counts at same position and combine them from diff fams
[lftempuniq,iuniq,~] = unique(lftemp(:,5:6),'rows','stable');
for i = 1:size(lftempuniq,1)
    [idup,~] = find(lftemp(:,5)==lftempuniq(i,1) & lftemp(:,6)==lftempuniq(i,2));
    ctotlf(i,1) = length(idup);   % total count of lf
end
if sum(ctotlf) == size(lftemp,1)
   lfrelacount = [lftemp(iuniq,1:2) lftemp(iuniq,5:7) ctotlf];
else
    disp('Inconsistency of lf total counts number')
end


%% read the detections checked by additional fam, choose KLNB only
SUFFIXhf = strcat('up.hf.time.',num2str(winlenhf),'_',num2str(winlenlf),'.',num2str(lofflf),...
    '.',num2str(ccminlf));
fname = strcat(rstpath, '/evtloc.KLNBcheck.allfam.dcutnodou.',SUFFIXhf);
hfaddcheck = load(fname);
% 25 cols, format is:
%   E(043) N(043) E(own) N(own) dist(own) lon lat dep off12 off13 off12sec off13sec date
%   main_arrival_time  cent_of_win  avecc ccadd1 ccadd2 ccadd3 ccadd4 ampmax ampmin eng normeng famnum

SUFFIXlf = strcat('lf.time.',num2str(winlenhf),'_',num2str(winlenlf),'.',num2str(lofflf),...
    '.',num2str(ccminlf));
fname = strcat(rstpath, '/evtloc.KLNBcheck.allfam.dcutnodou.',SUFFIXlf);
lfaddcheck = load(fname);



%% plot the acculative hit count
%%% figure 1, plot the overall hf and lf detections
% define and position the figure frame and axes of each plot 
f.fig = figure;
widin = 8;  % maximum width allowed is 8.5 inches
htin = 4;   % maximum height allowed is 11 inches
set(f.fig,'Position',[1*scrsz(3)/20 scrsz(4)/10 widin*res htin*res]);
nrow = 1;
ncol = 2;
for isub = 1:nrow*ncol
    f.ax(isub) = subplot(nrow,ncol,isub);
end 

% reposition
for isub = 1:nrow*ncol
    if mod(isub,ncol)==0
        icol = ncol;
        irow = floor(isub/ncol);
    else
        icol = mod(isub,ncol);
        irow = floor(isub/ncol)+1;
    end
    set(f.ax(isub),'Position', [(icol-0.85)/ncol, 1-irow/nrow*0.9, 1/ncol*0.8 1/nrow*0.8]);
    f.ax(isub).Box = 'on';
    grid(f.ax(isub), 'on');
    f.ax(isub).GridLineStyle = '--';
end

% marker size
% msizehf = 4;
% msizelf = 8;
msizehf = 2;
msizelf = 8;

% subplot 1 of figure 1
hold(f.ax(1),'on');
dumhf = hfrelacount;
dumhf(dumhf(:,6)>1, :) = [];
scatter(f.ax(1),dumhf(:,1),dumhf(:,2), msizehf, log10(dumhf(:,6)),'o','linew',0.2);  %, 'MarkerEdgeColor', 'w')
dumhf = sortrows(hfrelacount,6);
dumhf(dumhf(:,6)==1, :) = [];
scatter(f.ax(1),dumhf(:,1),dumhf(:,2), msizehf, log10(dumhf(:,6)), 'filled','o');  %, 'MarkerEdgeColor', 'w')
scatter(f.ax(1),relacont(:,1),relacont(:,2),20,'k','o','LineWidth',1);
plot(f.ax(1),[-100 100],[0 0],'k--');
plot(f.ax(1),[0 0],[-100 100],'k--');
text(f.ax(1),0.05,0.93,'a','FontSize',12,'unit','normalized','horizontalalignment','center',...
     'EdgeColor','k','Margin',2);
text(f.ax(1),0.93,0.93,'LZB','FontSize',12,'unit','normalized','horizontalalignment','center',...
     'EdgeColor','k','Margin',2);
text(f.ax(1),0.93,0.1,'HF','FontSize',15,'unit','normalized','horizontalalignment','center');
% colormap(f.ax(1),flipud(hot));
colormap(f.ax(1),'jet');
c=colorbar(f.ax(1),'NorthOutside');
c.Label.String = 'log_{10}(hits)';
c.Label.FontSize = 12;
% caxis(f.ax(1),[0 1.5]);
axis(f.ax(1), 'equal');
axis(f.ax(1),[-20 35 -25 20]);
f.ax(1).XLabel.String = 'E (km)';
f.ax(1).YLabel.String = 'N (km)';
hold(f.ax(1),'off');

% subplot 2 of figure 1
hold(f.ax(2),'on');
dumlf = lfrelacount;
dumlf(dumlf(:,6)>1, :) = [];
scatter(f.ax(2),dumlf(:,1),dumlf(:,2), msizelf, log10(dumlf(:,6)),'o','linew',0.2);   %, 'MarkerEdgeColor', 'w')
dumlf = sortrows(lfrelacount,6);
dumlf(dumlf(:,6)==1, :) = [];
scatter(f.ax(2),dumlf(:,1),dumlf(:,2), msizelf, log10(dumlf(:,6)), 'filled','o');  %, 'MarkerEdgeColor', 'w')
scatter(f.ax(2),relacont(:,1),relacont(:,2),20,'k','o','LineWidth',1);
plot(f.ax(2),[-100 100],[0 0],'k--');
plot(f.ax(2),[0 0],[-100 100],'k--');
text(f.ax(2),0.05,0.93,'b','FontSize',12,'unit','normalized','horizontalalignment','center',...
     'EdgeColor','k','Margin',2);
text(f.ax(2),0.93,0.93,'LZB','FontSize',12,'unit','normalized','horizontalalignment','center',...
     'EdgeColor','k','Margin',2);
text(f.ax(2),0.93,0.1,'LF','FontSize',15,'unit','normalized','horizontalalignment','center');
% colormap(f.ax(2),flipud(hot));
colormap(f.ax(2),'jet');
c=colorbar(f.ax(2),'NorthOutside');
c.Label.String = 'log_{10}(hits)';
c.Label.FontSize = 12;
% caxis(f.ax(2),[0 1.7]);
axis(f.ax(2), 'equal');
axis(f.ax(2),[-20 35 -25 20]);
f.ax(2).XLabel.String = 'E (km)';
f.ax(2).YLabel.String = 'N (km)';
hold(f.ax(2),'off');

print(f.fig,'-dpdf',strcat(rstpath,'/asd.rela.pdf'));

%%% save figure to file
% print(f.fig,'-dpdf',strcat(rstpath,'/allfam.hflf.alld.accuhit.',num2str(winlenhf),'_',num2str(winlenlf),...
%       '.',num2str(lofflf),'.',num2str(ccminlf),'.rela.pdf'));
% keyboard

%% try the entire catalog
hfall = hfrelatime;
lfall = lfrelatime;
rotang = 30;
% now the 1,2 col is changing from E(043) and N(043), to down-dip and strike
[hfall(:,1),hfall(:,2)] = coordinate_rot(hfall(:,1),hfall(:,2),rotang,0,0);
[lfall(:,1),lfall(:,2)] = coordinate_rot(lfall(:,1),lfall(:,2),rotang,0,0);

% sort them according to their occurring time
hfall = sortrows(hfall,[12,14]);
lfall = sortrows(lfall,[12,14]);

% obtain the relative time in hr
ih03 = find(hfall(:,12) < 2004*1000);
hfall(ih03,17) = (hfall(ih03,12)-2003060).*24+hfall(ih03,14)./3600;
ih04 = find(hfall(:,12) < 2005*1000 & hfall(:,12) > 2004*1000);
hfall(ih04,17) = (hfall(ih04,12)-2004194).*24+hfall(ih04,14)./3600;
ih05 = find(hfall(:,12) > 2005*1000);
hfall(ih05,17) = (hfall(ih05,12)-2005254).*24+hfall(ih05,14)./3600;

il03 = find(lfall(:,12) < 2004*1000);
lfall(il03,17) = (lfall(il03,12)-2003060).*24+lfall(il03,14)./3600;
il04 = find(lfall(:,12) < 2005*1000 & lfall(:,12) > 2004*1000);
lfall(il04,17) = (lfall(il04,12)-2004194).*24+lfall(il04,14)./3600;
il05 = find(lfall(:,12) > 2005*1000);
lfall(il05,17) = (lfall(il05,12)-2005254).*24+lfall(il05,14)./3600;


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
scatter(ax,hfall(ih03,17),hfall(ih03,1), msizehf, 'filled','o');  %, 'MarkerEdgeColor', 'w')
scatter(ax,lfall(il03,17),lfall(il03,1), msizehf, 'filled','ro');  %, 'MarkerEdgeColor', 'w'))
text(ax,0.05,0.93,'2003','FontSize',12,'unit','normalized','horizontalalignment','center',...
     'EdgeColor','k','Margin',2);
axis(ax,[0 60 -18 24]);
% axis(ax,[0 60 5 25]);
ax.YLabel.String = 'Along down-dip (km)';
% ax.YLabel.String = 'Along strike (km)';
hold(ax,'off');

ax = f7.ax(2);
hold(ax,'on');
scatter(ax,hfall(ih03,17),hfall(ih03,1), msizehf, 'filled','o');  %, 'MarkerEdgeColor', 'w')
scatter(ax,lfall(il03,17),lfall(il03,1), msizehf, 'filled','ro');  %, 'MarkerEdgeColor', 'w'))
axis(ax,[60 120 -18 24]);
% axis(ax,[60 120 5 25]);
% ax.YLabel.String = 'Along down-dip (km)';
% ax.YLabel.String = 'Along strike (km)';
hold(ax,'off');

ax = f7.ax(3);
hold(ax,'on');
scatter(ax,hfall(ih04,17),hfall(ih04,1), msizehf, 'filled','o');  %, 'MarkerEdgeColor', 'w')
scatter(ax,lfall(il04,17),lfall(il04,1), msizehf, 'filled','ro');  %, 'MarkerEdgeColor', 'w'))
text(ax,0.05,0.93,'2004','FontSize',12,'unit','normalized','horizontalalignment','center',...
     'EdgeColor','k','Margin',2);
axis(ax,[0 60 -18 24]);
% axis(ax,[0 60 5 25]);
% ax.YLabel.String = 'Along down-dip (km)';
% ax.YLabel.String = 'Along strike (km)';
hold(ax,'off');

ax = f7.ax(4);
hold(ax,'on');
scatter(ax,hfall(ih04,17),hfall(ih04,1), msizehf, 'filled','o');  %, 'MarkerEdgeColor', 'w')
scatter(ax,lfall(il04,17),lfall(il04,1), msizehf, 'filled','ro');  %, 'MarkerEdgeColor', 'w'))
axis(ax,[60 120 -18 24]);
% axis(ax,[60 120 5 25]);
% ax.YLabel.String = 'Along down-dip (km)';
% ax.YLabel.String = 'Along strike (km)';
hold(ax,'off');

ax = f7.ax(5);
hold(ax,'on');
scatter(ax,hfall(ih04,17),hfall(ih04,1), msizehf, 'filled','o');  %, 'MarkerEdgeColor', 'w')
scatter(ax,lfall(il04,17),lfall(il04,1), msizehf, 'filled','ro');  %, 'MarkerEdgeColor', 'w'))
axis(ax,[120 180 -18 24]);
% axis(ax,[120 180 5 25]);
% ax.YLabel.String = 'Along down-dip (km)';
% ax.YLabel.String = 'Along strike (km)';
hold(ax,'off');

ax = f7.ax(6);
hold(ax,'on');
scatter(ax,hfall(ih04,17),hfall(ih04,1), msizehf, 'filled','o');  %, 'MarkerEdgeColor', 'w')
scatter(ax,lfall(il04,17),lfall(il04,1), msizehf, 'filled','ro');  %, 'MarkerEdgeColor', 'w'))
axis(ax,[180 240 -18 24]);
% axis(ax,[180 240 5 25]);
% ax.YLabel.String = 'Along down-dip (km)';
% ax.YLabel.String = 'Along strike (km)';
hold(ax,'off');

ax = f7.ax(7);
hold(ax,'on');
scatter(ax,hfall(ih05,17),hfall(ih05,1), msizehf, 'filled','o');  %, 'MarkerEdgeColor', 'w')
scatter(ax,lfall(il05,17),lfall(il05,1), msizehf, 'filled','ro');  %, 'MarkerEdgeColor', 'w'))
text(ax,0.05,0.93,'2005','FontSize',12,'unit','normalized','horizontalalignment','center',...
     'EdgeColor','k','Margin',2);
axis(ax,[0 60 -18 24]);
% axis(ax,[0 60 5 25]);
% ax.YLabel.String = 'Along down-dip (km)';
% ax.YLabel.String = 'Along strike (km)';
hold(ax,'off');

ax = f7.ax(8);
hold(ax,'on');
scatter(ax,hfall(ih05,17),hfall(ih05,1), msizehf, 'filled','o');  %, 'MarkerEdgeColor', 'w')
scatter(ax,lfall(il05,17),lfall(il05,1), msizehf, 'filled','ro');  %, 'MarkerEdgeColor', 'w'))
axis(ax,[60 120 -18 24]);
% axis(ax,[60 120 5 25]);
% ax.YLabel.String = 'Along down-dip (km)';
% ax.YLabel.String = 'Along strike (km)';
hold(ax,'off');

ax = f7.ax(9);
hold(ax,'on');
scatter(ax,hfall(ih05,17),hfall(ih05,1), msizehf, 'filled','o');  %, 'MarkerEdgeColor', 'w')
scatter(ax,lfall(il05,17),lfall(il05,1), msizehf, 'filled','ro');  %, 'MarkerEdgeColor', 'w'))
axis(ax,[120 180 -18 24]);
% axis(ax,[120 180 5 25]);
% ax.YLabel.String = 'Along down-dip (km)';
% ax.YLabel.String = 'Along strike (km)';
hold(ax,'off');

ax = f7.ax(10);
hold(ax,'on');
scatter(ax,hfall(ih05,17),hfall(ih05,1), msizehf, 'filled','o');  %, 'MarkerEdgeColor', 'w')
scatter(ax,lfall(il05,17),lfall(il05,1), msizehf, 'filled','ro');  %, 'MarkerEdgeColor', 'w'))
axis(ax,[180 240 -18 24]);
% axis(ax,[180 240 5 25]);
% ax.YLabel.String = 'Along down-dip (km)';
% ax.YLabel.String = 'Along strike (km)';
ax.XLabel.String = 'Time (hr)';
hold(ax,'off');

print(f7.fig,'-depsc2',strcat(rstpath,'/alldcutdowndip.eps'));


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
scatter(ax,hfall(ih03,17),hfall(ih03,2), msizehf, 'filled','o');  %, 'MarkerEdgeColor', 'w')
scatter(ax,lfall(il03,17),lfall(il03,2), msizehf, 'filled','ro');  %, 'MarkerEdgeColor', 'w'))
text(ax,0.05,0.93,'2003','FontSize',12,'unit','normalized','horizontalalignment','center',...
     'EdgeColor','k','Margin',2);
axis(ax,[0 60 -36 35]);
% ax.YLabel.String = 'Along down-dip (km)';
ax.YLabel.String = 'Along strike (km)';
hold(ax,'off');

ax = f8.ax(2);
hold(ax,'on');
scatter(ax,hfall(ih03,17),hfall(ih03,2), msizehf, 'filled','o');  %, 'MarkerEdgeColor', 'w')
scatter(ax,lfall(il03,17),lfall(il03,2), msizehf, 'filled','ro');  %, 'MarkerEdgeColor', 'w'))
axis(ax,[60 120 -36 35]);
% ax.YLabel.String = 'Along down-dip (km)';
% ax.YLabel.String = 'Along strike (km)';
hold(ax,'off');

ax = f8.ax(3);
hold(ax,'on');
scatter(ax,hfall(ih04,17),hfall(ih04,2), msizehf, 'filled','o');  %, 'MarkerEdgeColor', 'w')
scatter(ax,lfall(il04,17),lfall(il04,2), msizehf, 'filled','ro');  %, 'MarkerEdgeColor', 'w'))
text(ax,0.05,0.93,'2004','FontSize',12,'unit','normalized','horizontalalignment','center',...
     'EdgeColor','k','Margin',2);
axis(ax,[0 60 -36 35]);
% ax.YLabel.String = 'Along down-dip (km)';
% ax.YLabel.String = 'Along strike (km)';
hold(ax,'off');

ax = f8.ax(4);
hold(ax,'on');
scatter(ax,hfall(ih04,17),hfall(ih04,2), msizehf, 'filled','o');  %, 'MarkerEdgeColor', 'w')
scatter(ax,lfall(il04,17),lfall(il04,2), msizehf, 'filled','ro');  %, 'MarkerEdgeColor', 'w'))
axis(ax,[60 120 -36 35]);
% ax.YLabel.String = 'Along down-dip (km)';
% ax.YLabel.String = 'Along strike (km)';
hold(ax,'off');

ax = f8.ax(5);
hold(ax,'on');
scatter(ax,hfall(ih04,17),hfall(ih04,2), msizehf, 'filled','o');  %, 'MarkerEdgeColor', 'w')
scatter(ax,lfall(il04,17),lfall(il04,2), msizehf, 'filled','ro');  %, 'MarkerEdgeColor', 'w'))
axis(ax,[120 180 -36 35]);
% ax.YLabel.String = 'Along down-dip (km)';
% ax.YLabel.String = 'Along strike (km)';
hold(ax,'off');

ax = f8.ax(6);
hold(ax,'on');
scatter(ax,hfall(ih04,17),hfall(ih04,2), msizehf, 'filled','o');  %, 'MarkerEdgeColor', 'w')
scatter(ax,lfall(il04,17),lfall(il04,2), msizehf, 'filled','ro');  %, 'MarkerEdgeColor', 'w'))
axis(ax,[180 240 -36 35]);
% ax.YLabel.String = 'Along down-dip (km)';
% ax.YLabel.String = 'Along strike (km)';
hold(ax,'off');

ax = f8.ax(7);
hold(ax,'on');
scatter(ax,hfall(ih05,17),hfall(ih05,2), msizehf, 'filled','o');  %, 'MarkerEdgeColor', 'w')
scatter(ax,lfall(il05,17),lfall(il05,2), msizehf, 'filled','ro');  %, 'MarkerEdgeColor', 'w'))
text(ax,0.05,0.93,'2005','FontSize',12,'unit','normalized','horizontalalignment','center',...
     'EdgeColor','k','Margin',2);
axis(ax,[0 60 -36 35]);
% ax.YLabel.String = 'Along down-dip (km)';
% ax.YLabel.String = 'Along strike (km)';
hold(ax,'off');

ax = f8.ax(8);
hold(ax,'on');
scatter(ax,hfall(ih05,17),hfall(ih05,2), msizehf, 'filled','o');  %, 'MarkerEdgeColor', 'w')
scatter(ax,lfall(il05,17),lfall(il05,2), msizehf, 'filled','ro');  %, 'MarkerEdgeColor', 'w'))
axis(ax,[60 120 -36 35]);
% ax.YLabel.String = 'Along down-dip (km)';
% ax.YLabel.String = 'Along strike (km)';
hold(ax,'off');

ax = f8.ax(9);
hold(ax,'on');
scatter(ax,hfall(ih05,17),hfall(ih05,2), msizehf, 'filled','o');  %, 'MarkerEdgeColor', 'w')
scatter(ax,lfall(il05,17),lfall(il05,2), msizehf, 'filled','ro');  %, 'MarkerEdgeColor', 'w'))
axis(ax,[120 180 -36 35]);
% ax.YLabel.String = 'Along down-dip (km)';
% ax.YLabel.String = 'Along strike (km)';
hold(ax,'off');

ax = f8.ax(10);
hold(ax,'on');
scatter(ax,hfall(ih05,17),hfall(ih05,2), msizehf, 'filled','o');  %, 'MarkerEdgeColor', 'w')
scatter(ax,lfall(il05,17),lfall(il05,2), msizehf, 'filled','ro');  %, 'MarkerEdgeColor', 'w'))
axis(ax,[180 240 -36 35]);
% ax.YLabel.String = 'Along down-dip (km)';
% ax.YLabel.String = 'Along strike (km)';
ax.XLabel.String = 'Time (hr)';
hold(ax,'off');

print(f8.fig,'-depsc2',strcat(rstpath,'/alldcutstrike.eps'));


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


% events & tremor inside bound 1
[is,ion] = inpolygon(hfrelatime(:,1),hfrelatime(:,2),reg2(:,1),reg2(:,2));
isinreg2 = is | ion;
hfreg2 = hfrelatime(isinreg2 == 1, :);

[is,ion] = inpolygon(lfrelatime(:,1),lfrelatime(:,2),reg2(:,1),reg2(:,2));
isinreg2 = is | ion;
lfreg2 = lfrelatime(isinreg2 == 1, :);

% sort them according to their occurring time
hfreg1 = sortrows(hfreg1,[12,14]);
hfreg2 = sortrows(hfreg2,[12,14]);
lfreg1 = sortrows(lfreg1,[12,14]);
lfreg2 = sortrows(lfreg2,[12,14]);

dateall = unique(hfrelatime(:,12));

% obtain the relative time in hr
ih103 = find(hfreg1(:,12) < 2004*1000);
hfreg1(ih103,17) = (hfreg1(ih103,12)-2003060).*24+hfreg1(ih103,14)./3600;
ih104 = find(hfreg1(:,12) < 2005*1000 & hfreg1(:,12) > 2004*1000);
hfreg1(ih104,17) = (hfreg1(ih104,12)-2004194).*24+hfreg1(ih104,14)./3600;
ih105 = find(hfreg1(:,12) > 2005*1000);
hfreg1(ih105,17) = (hfreg1(ih105,12)-2005254).*24+hfreg1(ih105,14)./3600;

il103 = find(lfreg1(:,12) < 2004*1000);
lfreg1(il103,17) = (lfreg1(il103,12)-2003060).*24+lfreg1(il103,14)./3600;
il104 = find(lfreg1(:,12) < 2005*1000 & lfreg1(:,12) > 2004*1000);
lfreg1(il104,17) = (lfreg1(il104,12)-2004194).*24+lfreg1(il104,14)./3600;
il105 = find(lfreg1(:,12) > 2005*1000);
lfreg1(il105,17) = (lfreg1(il105,12)-2005254).*24+lfreg1(il105,14)./3600;

ih203 = find(hfreg2(:,12) < 2004*1000);
hfreg2(ih203,17) = (hfreg2(ih203,12)-2003060).*24+hfreg2(ih203,14)./3600;
ih204 = find(hfreg2(:,12) < 2005*1000 & hfreg2(:,12) > 2004*1000);
hfreg2(ih204,17) = (hfreg2(ih204,12)-2004194).*24+hfreg2(ih204,14)./3600;
ih205 = find(hfreg2(:,12) > 2005*1000);
hfreg2(ih205,17) = (hfreg2(ih205,12)-2005254).*24+hfreg2(ih205,14)./3600;

il203 = find(lfreg2(:,12) < 2004*1000);
lfreg2(il203,17) = (lfreg2(il203,12)-2003060).*24+lfreg2(il203,14)./3600;
il204 = find(lfreg2(:,12) < 2005*1000 & lfreg2(:,12) > 2004*1000);
lfreg2(il204,17) = (lfreg2(il204,12)-2004194).*24+lfreg2(il204,14)./3600;
il205 = find(lfreg2(:,12) > 2005*1000);
lfreg2(il205,17) = (lfreg2(il205,12)-2005254).*24+lfreg2(il205,14)./3600;

f3.fig = figure;
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

tmr = lfreg2;
i03 = il203;
i04 = il204;
i05 = il205;

ax = f3.ax(1);
hold(ax,'on');
scatter(ax,tmr(i03,17),tmr(i03,1), msizehf, 'filled','o');  %, 'MarkerEdgeColor', 'w')
text(ax,0.05,0.93,'2003','FontSize',12,'unit','normalized','horizontalalignment','center',...
     'EdgeColor','k','Margin',2);
axis(ax,[0 60 -20 0]);
axis(ax,[0 60 5 25]);
ax.YLabel.String = 'E (km)';
hold(ax,'off');

ax = f3.ax(2);
hold(ax,'on');
scatter(ax,tmr(i03,17),tmr(i03,1), msizehf, 'filled','o');  %, 'MarkerEdgeColor', 'w')
axis(ax,[60 120 -20 0]);
axis(ax,[60 120 5 25]);
ax.YLabel.String = 'E (km)';
hold(ax,'off');

ax = f3.ax(3);
hold(ax,'on');
scatter(ax,tmr(i04,17),tmr(i04,1), msizehf, 'filled','o');  %, 'MarkerEdgeColor', 'w')
text(ax,0.05,0.93,'2004','FontSize',12,'unit','normalized','horizontalalignment','center',...
     'EdgeColor','k','Margin',2);
axis(ax,[0 60 -20 0]);
axis(ax,[0 60 5 25]);
ax.YLabel.String = 'E (km)';
hold(ax,'off');

ax = f3.ax(4);
hold(ax,'on');
scatter(ax,tmr(i04,17),tmr(i04,1), msizehf, 'filled','o');  %, 'MarkerEdgeColor', 'w')
axis(ax,[60 120 -20 0]);
axis(ax,[60 120 5 25]);
ax.YLabel.String = 'E (km)';
hold(ax,'off');

ax = f3.ax(5);
hold(ax,'on');
scatter(ax,tmr(i04,17),tmr(i04,1), msizehf, 'filled','o');  %, 'MarkerEdgeColor', 'w')
axis(ax,[120 180 -20 0]);
axis(ax,[120 180 5 25]);
ax.YLabel.String = 'E (km)';
hold(ax,'off');

ax = f3.ax(6);
hold(ax,'on');
scatter(ax,tmr(i04,17),tmr(i04,1), msizehf, 'filled','o');  %, 'MarkerEdgeColor', 'w')
axis(ax,[180 240 -20 0]);
axis(ax,[180 240 5 25]);
ax.YLabel.String = 'E (km)';
hold(ax,'off');

ax = f3.ax(7);
hold(ax,'on');
scatter(ax,tmr(i05,17),tmr(i05,1), msizehf, 'filled','o');  %, 'MarkerEdgeColor', 'w')
text(ax,0.05,0.93,'2005','FontSize',12,'unit','normalized','horizontalalignment','center',...
     'EdgeColor','k','Margin',2);
axis(ax,[0 60 -20 0]);
axis(ax,[0 60 5 25]);
ax.YLabel.String = 'E (km)';
hold(ax,'off');

ax = f3.ax(8);
hold(ax,'on');
scatter(ax,tmr(i05,17),tmr(i05,1), msizehf, 'filled','o');  %, 'MarkerEdgeColor', 'w')
axis(ax,[60 120 -20 0]);
axis(ax,[60 120 5 25]);
ax.YLabel.String = 'E (km)';
hold(ax,'off');

ax = f3.ax(9);
hold(ax,'on');
scatter(ax,tmr(i05,17),tmr(i05,1), msizehf, 'filled','o');  %, 'MarkerEdgeColor', 'w')
axis(ax,[120 180 -20 0]);
axis(ax,[120 180 5 25]);
ax.YLabel.String = 'E (km)';
hold(ax,'off');

ax = f3.ax(10);
hold(ax,'on');
scatter(ax,tmr(i05,17),tmr(i05,1), msizehf, 'filled','o');  %, 'MarkerEdgeColor', 'w')
% axis(ax,[180 240 -20 0]);
axis(ax,[180 240 5 25]);
ax.YLabel.String = 'E (km)';
ax.XLabel.String = 'Time (hr)';
hold(ax,'off');

print(f3.fig,'-dpdf',strcat(rstpath,'/lf2E.rela.pdf'));





%% select a date to plot the time evolution 

% SUFFIXhf = strcat('hf.time.',num2str(winlenhf),'_',num2str(winlenlf),'.',num2str(lofflf),...
%                   '.',num2str(ccminlf));
% fname = strcat(rstpath, '/evtloc.allfam.',SUFFIXhf);
% hfmaptime = load(fname);
% % 12 cols, format is:
% %   lon lat dep off12 off13 off12sec off13sec date strongest_arr_time center_of_win cc fam


% select a date to plot the general migration pattern
dateall = unique(hfrelatime(:,12));
refyr = 2005;
refjday = 254;
tmp = jul2dat(refyr,refjday);   % return month, day, year
refday = tmp(2);    
daymax = 15;
selhf = hfrelatime(hfrelatime(:,12)>=refyr*1000, :);
mindate = dateall(find(dateall>=refyr*1000, 1,'first'));
timeint = selhf(:,12)-(refyr*1000+refjday)+refday;
timedec = selhf(:,13)/3600/24;
time = timeint+timedec;
selhf(:,17) = time;

sellf = lfrelatime(lfrelatime(:,12)>=refyr*1000, :);
timeint = sellf(:,12)-(refyr*1000+refjday)+refday;
timedec = sellf(:,13)/3600/24;
time = timeint+timedec;
sellf(:,17) = time;

selhfpeng = selhf(selhf(:,17)<=daymax,:);
sellfpeng = sellf(sellf(:,17)<=daymax,:);

%%% define and position the figure frame and axes of each plot 
f2.fig = figure;
widin = 8;  % maximum width allowed is 8.5 inches
htin = 8;   % maximum height allowed is 11 inches
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
dumhf = sortrows(selhfpeng,-12);
scatter(f2.ax(1),dumhf(:,1),dumhf(:,2), msizehf, dumhf(:,end), 'filled','o');
scatter(f2.ax(1),relacont(:,1),relacont(:,2),20,'k','filled','o');
colormap(f2.ax(1),'jet');
c=colorbar(f2.ax(1),'NorthOutside');
c.Label.String = 'Day in Sept. 2005';
c.Label.FontSize = 12;
caxis(f2.ax(1),[refday daymax]);
plot(f2.ax(1),[-100 100],[0 0],'k--');
plot(f2.ax(1),[0 0],[-100 100],'k--');
text(f2.ax(1),0.05,0.93,'a','FontSize',12,'unit','normalized','horizontalalignment','center',...
     'EdgeColor','k','Margin',2);
text(f2.ax(1),0.93,0.93,'LZB','FontSize',12,'unit','normalized','horizontalalignment','center',...
     'EdgeColor','k','Margin',2);
text(f2.ax(1),0.93,0.1,'HF','FontSize',15,'unit','normalized','horizontalalignment','center');
f2.ax(1).Box = 'on';
grid(f2.ax(1), 'on');
f2.ax(1).GridLineStyle = '--';
axis(f2.ax(1), 'equal');
axis(f2.ax(1),[-20 35 -25 20]);
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
c.Label.FontSize = 12;
caxis(f2.ax(2),[refday daymax]);
plot(f2.ax(2),[-100 100],[0 0],'k--');
plot(f2.ax(2),[0 0],[-100 100],'k--');
text(f2.ax(2),0.05,0.93,'b','FontSize',12,'unit','normalized','horizontalalignment','center',...
     'EdgeColor','k','Margin',2);
text(f2.ax(2),0.93,0.93,'LZB','FontSize',12,'unit','normalized','horizontalalignment','center',...
     'EdgeColor','k','Margin',2);
text(f2.ax(2),0.93,0.1,'HF','FontSize',15,'unit','normalized','horizontalalignment','center');
f2.ax(2).Box = 'on';
grid(f2.ax(2), 'on');
f2.ax(2).GridLineStyle = '--';
axis(f2.ax(2), 'equal');
axis(f2.ax(2),[-20 35 -25 20]);
f2.ax(2).XLabel.String = 'E (km)';
f2.ax(2).YLabel.String = 'N (km)';
hold(f2.ax(2),'off');

% subplot 3 of figure 1
hold(f2.ax(3),'on');
dumlf = sortrows(sellfpeng,-12);
scatter(f2.ax(3),dumlf(:,1),dumlf(:,2), msizelf, dumlf(:,end), 'filled','o');  %, 'MarkerEdgeColor', 'w')
scatter(f2.ax(3),relacont(:,1),relacont(:,2),20,'k','filled','o');
colormap(f2.ax(3),'jet');
colorbar(f2.ax(3),'NorthOutside');
caxis(f2.ax(3),[refday daymax]);
plot(f2.ax(3),[-100 100],[0 0],'k--');
plot(f2.ax(3),[0 0],[-100 100],'k--');
text(f2.ax(3),0.05,0.93,'c','FontSize',12,'unit','normalized','horizontalalignment','center',...
     'EdgeColor','k','Margin',2);
text(f2.ax(3),0.93,0.93,'LZB','FontSize',12,'unit','normalized','horizontalalignment','center',...
     'EdgeColor','k','Margin',2);
text(f2.ax(3),0.93,0.1,'LF','FontSize',15,'unit','normalized','horizontalalignment','center');
f2.ax(3).Box = 'on';
grid(f2.ax(3), 'on');
f2.ax(3).GridLineStyle = '--';
axis(f2.ax(3), 'equal');
axis(f2.ax(3),[-20 35 -25 20]);
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
text(f2.ax(4),0.05,0.93,'d','FontSize',12,'unit','normalized','horizontalalignment','center',...
     'EdgeColor','k','Margin',2);
text(f2.ax(4),0.93,0.93,'LZB','FontSize',12,'unit','normalized','horizontalalignment','center',...
     'EdgeColor','k','Margin',2);
text(f2.ax(4),0.93,0.1,'LF','FontSize',15,'unit','normalized','horizontalalignment','center');
f2.ax(4).Box = 'on';
grid(f2.ax(4), 'on');
f2.ax(4).GridLineStyle = '--';
axis(f2.ax(4), 'equal');
axis(f2.ax(4),[-20 35 -25 20]);
f2.ax(4).XLabel.String = 'E (km)';
f2.ax(4).YLabel.String = 'N (km)';
hold(f2.ax(4),'off');
    
%%% save figure
% shrink(f2.fig,1.25,1.25);
print(f2.fig,'-dpdf',strcat(rstpath,'/allfam.hflf.tevo',num2str(refyr),'.',num2str(winlenhf),'_',...
      num2str(winlenlf),'.',num2str(lofflf),'.',num2str(ccminlf),'.rela.pdf'));
% close(f2.fig)    



















