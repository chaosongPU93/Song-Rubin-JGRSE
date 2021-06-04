% function plt_fam_detection
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% THIS is the script to plot the detections from each fam before merging them,
% or applying any distance cutoff. The distribution of CC can imply the overall
% reliability of the correction, rotation parameters.
%
% this is to get a sense how far can detections from a fam be trusted. For example, the same spot
% is detected by and closer to fam 1, but fam 2 generally has a higher CC at this location,
% though it is farther away. In other words, the correction, rotations based on fam 2 is more
% reliable here. Therefore, if a time window (one detection) at a spot is detected by fam 1 and
% closer to 1, but not detected by other families (no duplicates from other fams) we preserve it
% only when fam 1 generally has a highest CC here, otherwise we throw it away.
%
%
% Modified by Chao Song, chaosong@princeton.edu
% First created date:   2019/11/25
% Last modified date:   2019/11/25
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

% convert absolute loc to relative loc to its own lfe fam
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

%% load detections from all fams before merging and without distance cutoff
%%% ie. this might contain duplicates from other fams, and no distance cutoff has been applied
SUFFIXhf = strcat('up.hf.time.',num2str(winlenhf),'_',num2str(winlenlf),'.',num2str(lofflf),...
    '.',num2str(ccminlf));
fname = strcat(rstpath, '/evtloc.allfam.nodcut.',SUFFIXhf);
% fname = strcat('/home/data2/chaosong/Seisbasics/hypoinverse/testslab/evtloc.oslab.allfam.',SUFFIXhf);
hfnocut = load(fname);
% 23 cols, format is:
%   E(043) N(043) E(own) N(own) dist(own) lon lat dep off12 off13 off12sec off13sec date
%   main_arrival_time  cent_of_win  avecc ccadd1 ccadd2 ccadd3 ccadd4 ampmax ampmin famnum


SUFFIXlf = strcat('lf.time.',num2str(winlenhf),'_',num2str(winlenlf),'.',num2str(lofflf),...
    '.',num2str(ccminlf));
% fname = strcat(rstpath, '/evtloc.allfam.',SUFFIXlf);
fname = strcat(rstpath, '/evtloc.allfam.nodcut.',SUFFIXlf);
% fname = strcat('/home/data2/chaosong/Seisbasics/hypoinverse/testslab/evtloc.oslab.allfam.',SUFFIXlf);
lfnocut = load(fname);

%% plot the original accumulative detections of each fam, no dist cutoff, color-coded by the median CC
%%% this is to get a sense how far can detections from a fam be trusted. For example, the same spot
%%% is detected by and closer to fam 1, but fam 2 generally has a higher CC at this location,
%%% though it is farther away. In other words, the correction, rotations based on fam 2 is more
%%% reliable here. Therefore, if a time window (one detection) at a spot is detected by fam 1 and
%%% closer to 1, but not detected by other families (no duplicates from other fams) we preserve it
%%% only when fam 1 generally has a highest CC here, otherwise we throw it away.
for ifam = 1: size(nfampool,1)
    fam = nfampool(ifam,:);
    hffam = hfnocut(hfnocut(:,end)==str2double(fam),:);
    lffam = lfnocut(lfnocut(:,end)==str2double(fam),:);
    
   
    %%% plot the accumulative hit count in absolute locations inverted by hypoinverse
    %%% this would contain the effect from velocity model and slab model
    % sort the results according to their absolute locations
    hftemp = sortrows(hffam, [6, 7]);
    lftemp = sortrows(lffam, [6, 7]);
    
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
        hfcount = [hftemp(iuniq,1:8) ctothf medcchf];     % this would give 8+1+5=14 cols
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
        lfcount = [lftemp(iuniq,1:8) ctotlf medcclf];   % this would give 8+1+5=14 cols
    else
        disp('Inconsistency of lf total counts number')
    end
    
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
    msizelf = 8;
    
    % subplot 1 of figure 1
    hold(f.ax(1),'on');
    dumhf = hfcount;
    dumhf(dumhf(:,9)>1, :) = [];
    scatter(f.ax(1),dumhf(:,1),dumhf(:,2), msizehf, log10(dumhf(:,9)),'o','linew',0.2);  %, 'MarkerEdgeColor', 'w')
    dumhf = sortrows(hfcount,9);
    dumhf(dumhf(:,9)==1, :) = [];
    scatter(f.ax(1),dumhf(:,1),dumhf(:,2), msizehf, log10(dumhf(:,9)), 'filled','o');  %, 'MarkerEdgeColor', 'w')
    scatter(f.ax(1),relacont(:,1),relacont(:,2),20,[0.6 0.6 0.6],'o','LineWidth',1);
    scatter(f.ax(1),relacont(ifam,1),relacont(ifam,2),20,'k','o','LineWidth',1);
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
    dumlf = lfcount;
    dumlf(dumlf(:,9)>1, :) = [];
    scatter(f.ax(2),dumlf(:,1),dumlf(:,2), msizelf, log10(dumlf(:,9)),'o','linew',0.2);   %, 'MarkerEdgeColor', 'w')
    dumlf = sortrows(lfcount,9);
    dumlf(dumlf(:,9)==1, :) = [];
    scatter(f.ax(2),dumlf(:,1),dumlf(:,2), msizelf, log10(dumlf(:,9)), 'filled','o');  %, 'MarkerEdgeColor', 'w')
    scatter(f.ax(2),relacont(:,1),relacont(:,2),20,[0.6 0.6 0.6],'o','LineWidth',1);
    scatter(f.ax(2),relacont(ifam,1),relacont(ifam,2),20,'k','o','LineWidth',1);
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
    dumhf = sortrows(hfcount,10);
    scatter(f.ax(3),dumhf(:,1),dumhf(:,2), msizehf, dumhf(:,10), 'filled','o');  %, 'MarkerEdgeColor', 'w')
    scatter(f.ax(3),relacont(:,1),relacont(:,2),20,[0.6 0.6 0.6],'o','LineWidth',1);
    scatter(f.ax(3),relacont(ifam,1),relacont(ifam,2),20,'k','o','LineWidth',1);
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
    dumlf = sortrows(lfcount,10);
    scatter(f.ax(4),dumlf(:,1),dumlf(:,2), msizelf, dumlf(:,10), 'filled','o');  %, 'MarkerEdgeColor', 'w')
    scatter(f.ax(4),relacont(:,1),relacont(:,2),20,[0.6 0.6 0.6],'o','LineWidth',1);
    scatter(f.ax(4),relacont(ifam,1),relacont(ifam,2),20,'k','o','LineWidth',1);
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

    supertit(f.ax,strcat(fam,'--absloc--no dist cutoff'));
    
    print(f.fig,'-dpdf',strcat(rstpath,'/',fam,'.absloc.ownaccu.nodcut.pdf'));
    
    
    
    %%% plot the accumulative hit count in time offset, should be strictly uniform grid
    % sort the results according to their absolute locations
    hftemp2 = sortrows(hffam, [11, 12]);
    lftemp2 = sortrows(lffam, [11, 12]);
    
    % hf, obtain the counts at same position and combine them from diff fams
    [hftempuniq,iuniq,~] = unique(hftemp2(:,11:12),'rows','stable');
    ctothf2 = zeros(size(hftempuniq,1),1);
    medcchf2 = zeros(size(hftempuniq,1),5);
    for i = 1:size(hftempuniq,1)
        [idup,~] = find(hftemp2(:,11)==hftempuniq(i,1) & hftemp2(:,12)==hftempuniq(i,2));
        ctothf2(i,1) = length(idup);   % total count of hf
        medcchf2(i,1:5) = median(hftemp2(idup,16:20));   % this is median CC from trio and additional sta
    end
    if sum(ctothf2) == size(hftemp2,1)
        hfcount2 = [hftemp2(iuniq,11:12) ctothf2 medcchf2];     % this would give 8+1+5=14 cols
    else
        disp('Inconsistency of hf total counts number')
        disp(fam)
        break
    end
    
    % lf, obtain the counts at same position and combine them from diff fams
    [lftempuniq,iuniq,~] = unique(lftemp2(:,11:12),'rows','stable');
    ctotlf2 = zeros(size(lftempuniq,1),1);
    medcclf2 = zeros(size(lftempuniq,1),5);
    for i = 1:size(lftempuniq,1)
        [idup,~] = find(lftemp2(:,11)==lftempuniq(i,1) & lftemp2(:,12)==lftempuniq(i,2));
        ctotlf2(i,1) = length(idup);   % total count of hf
        medcclf2(i,1:5) = median(lftemp2(idup,16:20));   % this is median CC from trio and additional sta
    end
    if sum(ctotlf2) == size(lftemp2,1)
        lfcount2 = [lftemp2(iuniq,11:12) ctotlf2 medcclf2];   % this would give 8+1+5=14 cols
    else
        disp('Inconsistency of lf total counts number')
        disp(fam)
        break
    end
    
    % define and position the figure frame and axes of each plot
    f2.fig = figure;
    widin = 8;  % maximum width allowed is 8.5 inches
    htin = 10;   % maximum height allowed is 11 inches
    set(f2.fig,'Position',[1*scrsz(3)/20 scrsz(4)/10 widin*res htin*res]);
    nrow = 2;
    ncol = 2;
    for isub = 1:nrow*ncol
        f2.ax(isub) = subplot(nrow,ncol,isub);
        f2.ax(isub).Box = 'on';
        grid(f2.ax(isub), 'on');
        f2.ax(isub).GridLineStyle = '--';
        axis(f2.ax(isub), 'equal');
    end

    % reposition
    set(f2.ax(3),'Position', [0.1, 0.1, 0.4, 0.4]);
    set(f2.ax(4),'Position', [0.55, 0.1, 0.4, 0.4]);
    set(f2.ax(1),'Position', [0.1, 0.55, 0.4, 0.4]);
    set(f2.ax(2),'Position', [0.55, 0.55, 0.4, 0.4]);

    % marker size
    msizehf = 1;
    msizelf = 8;
    
    % subplot 1 of figure 1
    hold(f2.ax(1),'on');
    dumhf = hfcount2;
    dumhf(dumhf(:,3)>1, :) = [];
    scatter(f2.ax(1),dumhf(:,1),dumhf(:,2), msizehf, log10(dumhf(:,3)),'o','linew',0.2);  %, 'MarkerEdgeColor', 'w')
    dumhf = sortrows(hfcount2,3);
    dumhf(dumhf(:,3)==1, :) = [];
    scatter(f2.ax(1),dumhf(:,1),dumhf(:,2), msizehf, log10(dumhf(:,3)), 'filled','o');  %, 'MarkerEdgeColor', 'w')
    scatter(f2.ax(1),0,0,20,'k','o','LineWidth',1);
    plot(f2.ax(1),[-100 100],[0 0],'k--');
    plot(f2.ax(1),[0 0],[-100 100],'k--');
    text(f2.ax(1),0.05,0.93,'a','FontSize',12,'unit','normalized','horizontalalignment','center',...
        'EdgeColor','k','Margin',2);
    text(f2.ax(1),0.93,0.93,'LZB','FontSize',12,'unit','normalized','horizontalalignment','center',...
        'EdgeColor','k','Margin',2);
    text(f2.ax(1),0.93,0.1,'HF','FontSize',15,'unit','normalized','horizontalalignment','center');
    colormap(f2.ax(1),'jet');
    c=colorbar(f2.ax(1),'SouthOutside');
    c.Label.String = 'log_{10}(hits)';
    c.Label.FontSize = 12;
    % caxis(f2.ax(1),[0 1.5]);
    f2.ax(1).YLabel.String = 'STA13 offset (s)';
    axis(f2.ax(1),[-0.5 0.5 -0.5 0.5]);
    hold(f2.ax(1),'off');
    
    % subplot 2 of figure 1
    hold(f2.ax(2),'on');
    dumlf = lfcount2;
    dumlf(dumlf(:,3)>1, :) = [];
    scatter(f2.ax(2),dumlf(:,1),dumlf(:,2), msizelf, log10(dumlf(:,3)),'o','linew',0.2);   %, 'MarkerEdgeColor', 'w')
    dumlf = sortrows(lfcount2,3);
    dumlf(dumlf(:,3)==1, :) = [];
    scatter(f2.ax(2),dumlf(:,1),dumlf(:,2), msizelf, log10(dumlf(:,3)), 'filled','o');  %, 'MarkerEdgeColor', 'w')
    scatter(f2.ax(2),0,0,20,'k','o','LineWidth',1);
    plot(f2.ax(2),[-100 100],[0 0],'k--');
    plot(f2.ax(2),[0 0],[-100 100],'k--');
    plot(f2.ax(2),[-.35 0 .35 .35 0 -.35 -.35],[-.35 -.35 0 .35 .35 0 -.35],'--','color',...
         [.5 .5 .5],'linew',1);
    text(f2.ax(2),0.05,0.93,'b','FontSize',12,'unit','normalized','horizontalalignment','center',...
        'EdgeColor','k','Margin',2);
    text(f2.ax(2),0.93,0.93,'LZB','FontSize',12,'unit','normalized','horizontalalignment','center',...
        'EdgeColor','k','Margin',2);
    text(f2.ax(2),0.93,0.1,'LF','FontSize',15,'unit','normalized','horizontalalignment','center');
    colormap(f2.ax(2),'jet');
    c=colorbar(f2.ax(2),'SouthOutside');
    c.Label.String = 'log_{10}(hits)';
    c.Label.FontSize = 12;
%     caxis(f2.ax(2),[0 1.7]);
    axis(f2.ax(2),[-1 1 -1 1]);
    hold(f2.ax(2),'off');
    
    % subplot 3 of figure 1
    hold(f2.ax(3),'on');
    dumhf = sortrows(hfcount2,4);
    scatter(f2.ax(3),dumhf(:,1),dumhf(:,2), msizehf, dumhf(:,4), 'filled','o');  %, 'MarkerEdgeColor', 'w')
    scatter(f2.ax(3),0,0,20,'k','o','LineWidth',1);
    plot(f2.ax(3),[-100 100],[0 0],'k--');
    plot(f2.ax(3),[0 0],[-100 100],'k--');
    text(f2.ax(3),0.05,0.93,'c','FontSize',12,'unit','normalized','horizontalalignment','center',...
        'EdgeColor','k','Margin',2);
    text(f2.ax(3),0.93,0.93,'LZB','FontSize',12,'unit','normalized','horizontalalignment','center',...
        'EdgeColor','k','Margin',2);
    text(f2.ax(3),0.93,0.1,'HF','FontSize',15,'unit','normalized','horizontalalignment','center');
    colormap(f2.ax(3),'jet');
    c=colorbar(f2.ax(3),'NorthOutside');
    c.Label.String = 'Median CC';
    c.Label.FontSize = 12;
    caxis(f2.ax(3),[0.35 0.65]);
    f2.ax(3).XLabel.String = 'STA12 offset (s)';
    f2.ax(3).YLabel.String = 'STA13 offset (s)';
    axis(f2.ax(3),[-0.5 0.5 -0.5 0.5]);
    hold(f2.ax(3),'off');
    
    % subplot 4 of figure 1
    hold(f2.ax(4),'on');
    dumlf = sortrows(lfcount2,4);
    scatter(f2.ax(4),dumlf(:,1),dumlf(:,2), msizelf, dumlf(:,4), 'filled','o');  %, 'MarkerEdgeColor', 'w')
    scatter(f2.ax(4),0,0,20,'k','o','LineWidth',1);
    plot(f2.ax(4),[-100 100],[0 0],'k--');
    plot(f2.ax(4),[0 0],[-100 100],'k--');
    plot(f2.ax(4),[-.35 0 .35 .35 0 -.35 -.35],[-.35 -.35 0 .35 .35 0 -.35],'--','color',...
         [.5 .5 .5],'linew',1);
    text(f2.ax(4),0.05,0.93,'d','FontSize',12,'unit','normalized','horizontalalignment','center',...
        'EdgeColor','k','Margin',2);
    text(f2.ax(4),0.93,0.93,'LZB','FontSize',12,'unit','normalized','horizontalalignment','center',...
        'EdgeColor','k','Margin',2);
    text(f2.ax(4),0.93,0.1,'LF','FontSize',15,'unit','normalized','horizontalalignment','center');
    colormap(f2.ax(4),'jet');
    c=colorbar(f2.ax(4),'NorthOutside');
    c.Label.String = 'Median CC';
    c.Label.FontSize = 12;
    caxis(f2.ax(4),[0.35 0.65]);
    f2.ax(4).XLabel.String = 'STA12 offset (s)';
    axis(f2.ax(4),[-1 1 -1 1]);
    hold(f2.ax(4),'off');

    supertit(f2.ax,strcat(fam,'--offtime--no dist cutoff'));
    
    print(f2.fig,'-dpdf',strcat(rstpath,'/',fam,'.offtime.ownaccu.nodcut.pdf'));
    
    
end



















