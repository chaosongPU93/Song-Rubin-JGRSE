function [f,barelse,pdfxlocelse,pdfvalelse] = ...
    plt_hist_combined_RTMs_otherrtm(indelse,ranvechf98,medallhf,angbest,xran,yran,vdistelse,wtelse,...
                           binw,conf,edgetype)
% this function is to simplify the plotting of RTM coverage, LF-HF offset distribution
% of combined RTMs and its corresponding confidence interval and its zoom-in, for two
% opposite propagation groups at one region
%
% 
%
% By Chao Song, chaosong@princeton.edu
% First created date:   2020/11/22
% Last modified date:   2020/11/22 

[scrsz, res] = pixelperinch(1);

f.fig=figure;
widin = 8;  % maximum width allowed is 8.5 inches
htin = 5;   % maximum height allowed is 11 inches
set(f.fig,'Position',[scrsz(1)+1*scrsz(3)/10 scrsz(2)+scrsz(4)/20 widin*res htin*res]);
nrow = 2;
ncol = 2;
for isub = 1:nrow*ncol-1
    f.ax(isub) = subplot(nrow,ncol,isub);
end

%%% reposition
set(f.ax(1), 'position', [ 0.1, 0.1, 0.4375, 0.85]);    % 0.35

% subplot 1
hold(f.ax(1),'on');
ax = f.ax(1);
plot(ax,[-100 100],[0 0],'k--');
plot(ax,[0 0],[-100 100],'k--');
ax.FontSize = 9;
axis(ax, 'equal');
for i = 1: length(indelse)
    [rotx1, roty1] = complex_rot(0, ranvechf98(indelse(i),2)-medallhf(indelse(i),3), ...
                                 -angbest(indelse(i)));
    [rotx2, roty2] = complex_rot(0, medallhf(indelse(i),3)-ranvechf98(indelse(i),1), ...
                                 -angbest(indelse(i)));
    xvect = [medallhf(indelse(i),1)-rotx2 medallhf(indelse(i),1)+rotx1];
    yvect = [medallhf(indelse(i),2)-roty2 medallhf(indelse(i),2)+roty1];
    drawArrow(ax,xvect,yvect,xran,yran,'linewidth',1,'linestyle','-','color',[0.7 0.7 0.7]);
end
for i = 1: length(indelse)
    scatter(ax,medallhf(indelse(i),1),medallhf(indelse(i),2), 30, indelse(i), 'filled','o');  %, 'MarkerEdgeColor', 'w')
end
colormap(ax,'jet');
c=colorbar(ax,'SouthOutside');
caxis(ax,[1 size(ranvechf98,1)]);
c.Label.String = strcat('RTM number in chronological order');
c.Label.FontSize = 11;
% text(ax,0.85,0.1,'HF','FontSize',12,'unit','normalized');
% text(ax,0.04,0.95,'a','FontSize',11,'unit','normalized','EdgeColor','k','Margin',2);
% text(ax,0.95,0.95,'LZB','FontSize',11,'unit','normalized','EdgeColor','k','Margin',2,...
%     'horizontalalignment','right');
ax.Box = 'on';
grid(ax, 'on');
ax.GridLineStyle = '--';
% ax.XAxisLocation = 'top';
xticks(ax,xran(1):5:xran(2));
yticks(ax,yran(1):5:yran(2));
axis(ax,[xran yran]);
xlabel(ax,'E (km)','fontsize',11);
ylabel(ax,'N (km)','fontsize',11);
text(ax,0.5,0.06,'All detections','FontSize',12,'unit','normalized',...
     'horizontalalignment','center');
text(ax,0.5,0.12,'Other RTMs,','FontSize',12,'unit','normalized',...
     'horizontalalignment','center'); 
text(ax,0.04,0.94,'a','FontSize',11,'unit','normalized','EdgeColor','k','Margin',2);
text(ax,0.95,0.94,'TWKB','FontSize',11,'unit','normalized','EdgeColor','k','Margin',2,...
    'horizontalalignment','right');
hold(ax,'off');
ax.Position(2) = ax.Position(2)-0.05;
% c.Position(2) = c.Position(2)-0.1;

% subplot 2 
% binw = 1;

ywid = f.ax(1).Position(4)+f.ax(1).Position(2)-(c.Position(2)) - 0.068;
% ywid = f.ax(1).Position(4)-f.ax(1).Position(2)+(c.Position(2)) - 0.06;
set(f.ax(2), 'position', [ 0.62, c.Position(2), 0.35, ywid*2/3 ]);    % 0.55


meanelse = wt_mean(vdistelse,wtelse);
sigmaelse = sqrt(wt_var(vdistelse,wtelse,2));
neff = sum(wtelse);
% conf = 99;
CIelse = confidence_interval(meanelse,sigmaelse,neff,conf);

hold(f.ax(2),'on');
ax = f.ax(2);
patarea = [0 0;
           -100 0;
           -100 1;
            0 1;
            0 0];
patch(ax,patarea(:,1),patarea(:,2),'k','Facealpha',0.05,'edgecolor','none');
[ax, barelse, pdfxlocelse, ~, ~, pdfvalelse] = plt_weighted_dist(ax, vdistelse, wtelse, binw,edgetype);
barelse(1).FaceColor = [0.7 0.7 0.7];
plot(ax,[meanelse meanelse],[-100 100],'k--','linew',1.5);    % median of wt dist
errorbar(ax,meanelse,0.17,CIelse(1)-meanelse,CIelse(2)-meanelse,'horizontal','o','markersize',3,'color',...
         'k','linew',1,'MarkerEdgeColor','k','MarkerFaceColor','k');
ylim(ax,[0 0.2]);
xlim(ax, [-12.5 12.5]); 
xvect = [-7 -12.5+0.5];
yvect = [0.015 0.015];
% drawArrow(ax,xvect,yvect,ax.XLim,ax.YLim,'linewidth',1.5,'linestyle','-',...
%           'color',[0.5 0.5 0.5]);
% arrow1 = annotation('arrow','linewidth',1.5,'linestyle','-','color',[0.5 0.5 0.5]);
% arrow1.Parent = f.ax(2);
% arrow1.Position = [-7, 0.015, -5, 0] ;
xvect = [7 12.5-0.5];
yvect = [0.015 0.015];
% drawArrow(ax,xvect,yvect,ax.XLim,ax.YLim,'linewidth',1.5,'linestyle','-',...
%           'color',[0.5 0.5 0.5]);
text(f.ax(2),0.05,0.15,strcat({'LF lags'}),'fontsize',11,'unit','normalized');
text(f.ax(2),0.96,0.15,strcat({'LF leads'}),'fontsize',11,'unit','normalized',...
     'horizontalalignment','right');

text(ax,0.85,0.62,num2str(length(vdistelse)),'fontsize',12,'unit','normalized',...
     'horizontalalignment','center','color','k');
text(ax,0.85,0.55,strcat({'detections'}),'fontsize',12,'unit','normalized',...
     'horizontalalignment','center','color','k');
% text(f.ax(2),0.04,0.93,'b','FontSize',11,'unit','normalized','EdgeColor','k','Margin',2);
% lgd = legend(ax,[barne,barsw],{'ENE propagation','WSW propagation'},'fontsize',7,...
%        'numcolumns',2,...
%        'Position',[0.62+0.35*0.048  c.Position(2)+ywid*2/3*0.93  0.35*0.9  ywid*2/3*0.05]);
hold(ax,'off');


set(f.ax(3), 'position', [ 0.62, c.Position(2)+ywid*2/3+0.05, 0.35, ywid*1/3 ]);    % 0.55
ax = f.ax(3);
hold(ax,'on');
patarea = [0 0;
           -100 0;
           -100 1;
            0 1;
            0 0];
patch(ax,patarea(:,1),patarea(:,2),'k','Facealpha',0.05,'edgecolor','none');
plot(ax,[0 0],[-100 100],'--','color',[0.6 0.6 0.6],'linew',1);
plot(ax,[meanelse meanelse],[-100 100],'k--','linew',2);    % median of wt dist
errorbar(ax,meanelse,0.17,CIelse(1)-meanelse,CIelse(2)-meanelse,'horizontal','o','markersize',6,'color',...
         'k','linewidth',1.5,'MarkerEdgeColor','k','MarkerFaceColor','k','CapSize',10);
text(ax,0.95,0.85,strcat({'10X'}),'fontsize',10,'unit','normalized','EdgeColor','k',...
     'Margin',2,'horizontalalignment','right');
text(ax,0.05,0.15,strcat({'99% CI'}),'fontsize',12,'unit','normalized','horizontalalignment','left'); 
text(ax,0.04,0.85,'b','FontSize',11,'unit','normalized','EdgeColor','k','Margin',2);
ylim(ax,[0.1 0.2]);
xticks(ax,-3: 0.5: 3);
yticks(ax,[0.1 0.15 0.2]);
% xlabel(ax,'LF-HF (km)','fontsize',11);
ylabel(ax,'PDF estimate','fontsize',11);
ax.Box = 'on';
grid(ax, 'on');
ax.GridLineStyle = '--';
hold(ax,'off');

xlimit = f.ax(2).XLim/10;
xlim(f.ax(3),xlimit); 
% convert data units to global figure units
[xs,ys] = ds2nfu(f.ax(3),xlimit(1),0.1);
[xe,ye] = ds2nfu(f.ax(2),xlimit(1),0.2);
% plot two dashed lines denoting the zoom-in effect
annotation('line',[xs,xe],[ys,ye],'color','k','linestyle','--');
% convert data units to global figure units
[xs,ys] = ds2nfu(f.ax(3),xlimit(2),0.1);
[xe,ye] = ds2nfu(f.ax(2),xlimit(2),0.2);
% plot two dashed lines denoting the zoom-in effect
annotation('line',[xs,xe],[ys,ye],'color','k','linestyle','--');









