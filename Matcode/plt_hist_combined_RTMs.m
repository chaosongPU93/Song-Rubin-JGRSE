function [f,barsw,pdfxlocsw,pdfvalsw,barne,pdfxlocne,pdfvalne,lgd] = ...
    plt_hist_combined_RTMs(indne,indsw,ranvechf98,medallhf,angbest,xran,yran,vdistsw,wtsw,...
                           vdistne,wtne,binw,conf,edgetype)
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
htin = 6;   % maximum height allowed is 11 inches
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
% xran = [-15 5];
% yran = [-15 15];
% indplt = union(indsw,indne);    % index to plot
% indplt = indne;    % index to plot
[ax,c] = plt_mig_prop_onmap(ax, indne,indsw, ranvechf98, medallhf, angbest, xran,yran);
% text(f.ax(1),0.5,0.06,'All detections','FontSize',12,'unit','normalized',...
%      'horizontalalignment','center');
% text(f.ax(1),0.5,0.12,'Main LZB region,','FontSize',12,'unit','normalized',...
%      'horizontalalignment','center'); 
text(ax,0.04,0.95,'a','FontSize',11,'unit','normalized','EdgeColor','k','Margin',2);
text(ax,0.95,0.95,'TWKB','FontSize',11,'unit','normalized','EdgeColor','k','Margin',2,...
    'horizontalalignment','right');
hold(ax,'off');
ax.Position(2) = ax.Position(2)-0.05;
% c.Position(2) = c.Position(2)-0.1;

% subplot 2 
% binw = 1;

ywid = f.ax(1).Position(4)+f.ax(1).Position(2)-(c.Position(2)) - 0.05;
% ywid = f.ax(1).Position(4)-f.ax(1).Position(2)+(c.Position(2)) - 0.06;
set(f.ax(2), 'position', [ 0.62, c.Position(2), 0.35, ywid*2/3 ]);    % 0.55

meansw = wt_mean(vdistsw,wtsw);
sigmasw = sqrt(wt_var(vdistsw,wtsw,2));
neff = sum(wtsw);
% conf = 99;
CIsw = confidence_interval(meansw,sigmasw,neff,conf);
meanne = wt_mean(vdistne,wtne);
sigmane = sqrt(wt_var(vdistne,wtne,2));
neff = sum(wtne);
% conf = 99;
CIne = confidence_interval(meanne,sigmane,neff,conf);

hold(f.ax(2),'on');
ax = f.ax(2);
[ax, barsw, pdfxlocsw, ~, ~, pdfvalsw] = plt_weighted_dist(ax, vdistsw, wtsw, binw,edgetype);
barsw(1).FaceAlpha = 1;
[ax, barne, pdfxlocne, ~, ~, pdfvalne] = plt_weighted_dist(ax, vdistne, wtne, binw,edgetype);
barne(1).FaceColor = [0 1 1];
plot(ax,[meansw meansw],[-100 100],'r--','linew',1.5);    % median of wt dist
plot(ax,[meanne meanne],[-100 100],'b--','linew',1.5);    % median of wt dist
errorbar(ax,meansw,0.26,CIsw(1)-meansw,CIsw(2)-meansw,'horizontal','o','markersize',3,'color',...
         'r','linew',1,'MarkerEdgeColor','r','MarkerFaceColor','r');
errorbar(ax,meanne,0.26,CIne(1)-meanne,CIne(2)-meanne,'horizontal','o','markersize',3,'color',...
         'b','linew',1,'MarkerEdgeColor','b','MarkerFaceColor','b');
% xran1 = max(abs(pdfxlocsw(1)-binw/2-1), abs(pdfxlocsw(end)+binw/2)+1);
% xran2 = max(abs(pdfxlocne(1)-binw/2-1), abs(pdfxlocne(end)+binw/2)+1);
% xran = max(xran1, xran2);
% xlim(ax, [-xran, xran]);
% xvect = [-xran+3.5 -xran+0.5];
% yvect = [0.015 0.015];
% drawArrow(ax,xvect,yvect,ax.XLim,ax.YLim,'linewidth',1.5,'linestyle','-',...
%           'color',[0.5 0.5 0.5]);
% arrow1 = annotation('arrow','linewidth',1.5,'linestyle','-','color',[0.5 0.5 0.5]);
% arrow1.Parent = f.ax(2);
% arrow1.Position = [-7, 0.015, -5, 0] ;
% xvect = [xran-3.5 xran-0.5];
% yvect = [0.015 0.015];
% drawArrow(ax,xvect,yvect,ax.XLim,ax.YLim,'linewidth',1.5,'linestyle','-',...
%           'color',[0.5 0.5 0.5]);
% arrow2 = annotation('arrow','linewidth',1.5,'linestyle','-','color',[0.5 0.5 0.5]);
% arrow2.Parent = f.ax(2);
% arrow2.Position = [7, 0.015, 5, 0] ;
% text(f.ax(2),0.05,0.1,strcat({'WSW'}),'fontsize',11,'unit','normalized');
% text(f.ax(2),0.95,0.1,strcat({'ENE'}),'fontsize',11,'unit','normalized',...
%      'horizontalalignment','right');
text(ax,0.2,0.6,num2str(length(vdistne)),'fontsize',12,'unit','normalized',...
     'horizontalalignment','center','color','b');
text(ax,0.2,0.54,strcat({'detections'}),'fontsize',12,'unit','normalized',...
     'horizontalalignment','center','color','k');
text(ax,0.85,0.6,num2str(length(vdistsw)),'fontsize',12,'unit','normalized',...
     'horizontalalignment','center','color','r');
text(ax,0.85,0.54,strcat({'detections'}),'fontsize',12,'unit','normalized',...
     'horizontalalignment','center','color','k');
% text(f.ax(2),0.04,0.93,'b','FontSize',11,'unit','normalized','EdgeColor','k','Margin',2);
lgd = legend(ax,[barne,barsw],{'ENE propagation','WSW propagation'},'fontsize',7,...
       'numcolumns',2,...
       'Position',[0.62+0.35*0.048  c.Position(2)+ywid*2/3*0.93  0.35*0.9  ywid*2/3*0.05]);
hold(ax,'off');

set(f.ax(3), 'position', [ 0.62, c.Position(2)+ywid*2/3+0.05, 0.35, ywid*1/3 ]);    % 0.55
ax = f.ax(3);
hold(ax,'on');
patarea = [0 0;
           -100 0;
           -100 1;
            0 1;
            0 0];
% patch(ax,patarea(:,1),patarea(:,2),'k','Facealpha',0.1,'edgecolor','none');
plot(ax,[0 0],[-100 100],'--','color',[0.6 0.6 0.6],'linew',1);
plot(ax,[meansw meansw],[-100 100],'r--','linew',2);    % median of wt dist
plot(ax,[meanne meanne],[-100 100],'b--','linew',2);    % median of wt dist
errorbar(ax,meansw,0.26,CIsw(1)-meansw,CIsw(2)-meansw,'horizontal','o','markersize',6,'color',...
         'r','linewidth',1.5,'MarkerEdgeColor','r','MarkerFaceColor','r','CapSize',10);
errorbar(ax,meanne,0.26,CIne(1)-meanne,CIne(2)-meanne,'horizontal','o','markersize',6,'color',...
         'b','linewidth',1.5,'MarkerEdgeColor','b','MarkerFaceColor','b','CapSize',10);
% text(ax,0.95,0.86,strcat({'10X'}),'fontsize',10,'unit','normalized','EdgeColor','k',...
%      'Margin',2,'horizontalalignment','right');
% text(ax,0.05,0.15,strcat({'99% CI'}),'fontsize',12,'unit','normalized','horizontalalignment','left'); 
text(ax,0.04,0.86,'b','FontSize',11,'unit','normalized','EdgeColor','k','Margin',2);
ylim(ax,[0.2 0.3]);
xticks(ax,-3: 0.5: 3);
yticks(ax,[0.2 0.25 0.3]);
% xlabel(ax,'LF-HF (km)','fontsize',11);
ylabel(ax,'PDF estimate','fontsize',11);
ax.Box = 'on';
grid(ax, 'on');
ax.GridLineStyle = '--';
hold(ax,'off');

% xlimit = f.ax(2).XLim/10;
% xlim(ax,xlimit); 
% % convert data units to global figure units
% [xs,ys] = ds2nfu(f.ax(3),xlimit(1),0.2);
% [xe,ye] = ds2nfu(f.ax(2),xlimit(1),0.3);
% % plot two dashed lines denoting the zoom-in effect
% annotation('line',[xs,xe],[ys,ye],'color','k','linestyle','--');
% % convert data units to global figure units
% [xs,ys] = ds2nfu(f.ax(3),xlimit(2),0.2);
% [xe,ye] = ds2nfu(f.ax(2),xlimit(2),0.3);
% % plot two dashed lines denoting the zoom-in effect
% annotation('line',[xs,xe],[ys,ye],'color','k','linestyle','--');

%%% save figure
% print(f.fig,'-dpdf',strcat(rstpath,'/vdist_unshifted_mainSW+NE_all',num2str(nfam),'.pdf'));