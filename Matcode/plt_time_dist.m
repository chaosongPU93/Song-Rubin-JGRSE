function [f] = plt_time_dist(hftime,hfdist,lftime,lfdist,rtmt,nrtmt,msizehf,msizelf,dt,mint,maxt,yran)
% 
% This function aims to ease the plot of time-distance of tremor
% detections in LF and HF in one ETS episode, in order to find 
% spatiotemporally clustered tremor bursts that would potentially
% be RTMs
%       
%
%
% By Chao Song, chaosong@princeton.edu
% First created date:   2020/12/23
% Last modified date:   2020/12/23

[scrsz, res] = pixelperinch(1);

f.fig = figure;
f.fig.Renderer = 'painters';
widin = 11;  % maximum width allowed is 8.5 inches
htin = 10;   % maximum height allowed is 11 inches
set(f.fig,'Position',[1*scrsz(3)/20 scrsz(4)/10 widin*res htin*res]);
% maxt = ceil(max(hftime));
% nrow = ceil((maxt-mint)/dt)*2;
nrow = 6;
ncol = 1;
for isub = 1:nrow*ncol
    f.ax(isub) = subplot(nrow,ncol,isub);
end

% reposition
for isub = 1:nrow*ncol
    f.ax(isub).Box = 'on';
    grid(f.ax(isub), 'on');
    f.ax(isub).GridLineStyle = '--';
end

% marker size
% msizehf = 2;
% msizelf = 4;

    
for i = 1: round(nrow/2)
    ax = f.ax(2*i-1);
    hold(ax,'on');
    for j = 1: size(rtmt,1)
        patarea = [rtmt(j,1) -100;
                   rtmt(j,2) -100;
                   rtmt(j,2) 100;
                   rtmt(j,1) 100;
                   rtmt(j,1) -100];
        patch(ax,patarea(:,1),patarea(:,2),'k','Facealpha',0.3,'edgecolor','none');
    end
    for j = 1: size(nrtmt,1)
        patarea = [nrtmt(j,1) -100;
                   nrtmt(j,2) -100;
                   nrtmt(j,2) 100;
                   nrtmt(j,1) 100;
                   nrtmt(j,1) -100];
        patch(ax,patarea(:,1),patarea(:,2),'b','Facealpha',0.3,'edgecolor','none');
    end
    scatter(ax, hftime, hfdist, msizehf, 'filled','ro');  %, 'MarkerEdgeColor', 'w')
    text(ax,0.05,0.9,'HF','FontSize',10,'unit','normalized','horizontalalignment','left',...
        'EdgeColor','k','Margin',2);
    xlim(ax,[mint+(i-1)*dt mint+i*dt]);
    ylim(ax, yran);
%     ax.YLabel.String = 'Along strike (km)';
    hold(ax,'off');

end

for i = 1: round(nrow/2)
    ax = f.ax(2*i);
    hold(ax,'on');
    for j = 1: size(rtmt,1)
        patarea = [rtmt(j,1) -100;
                   rtmt(j,2) -100;
                   rtmt(j,2) 100;
                   rtmt(j,1) 100;
                   rtmt(j,1) -100];
        patch(ax,patarea(:,1),patarea(:,2),'k','Facealpha',0.3,'edgecolor','none');
    end
    for j = 1: size(nrtmt,1)
        patarea = [nrtmt(j,1) -100;
                   nrtmt(j,2) -100;
                   nrtmt(j,2) 100;
                   nrtmt(j,1) 100;
                   nrtmt(j,1) -100];
        patch(ax,patarea(:,1),patarea(:,2),'b','Facealpha',0.3,'edgecolor','none');
    end
    scatter(ax, lftime, lfdist, msizelf, 'filled','ro');  %, 'MarkerEdgeColor', 'w'))
    text(ax,0.05,0.9,'LF','FontSize',10,'unit','normalized','horizontalalignment','left',...
        'EdgeColor','k','Margin',2);
    xlim(ax,[mint+(i-1)*dt mint+i*dt]);
    ylim(ax, yran);
%     ax.YLabel.String = 'Along strike (km)';
    hold(ax,'off');

end

% keyboard