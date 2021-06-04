function [ax, barHdl, pdfxloc, ctval, probval, pdfval] = plt_weighted_dist(ax, vdist, wt, binw, edgetype)
% This function is to plot the distribution normalized as PDF of distances
% and the corresponding weights. 
%   Served more as a simplification in the main code  
%
% By Chao Song, chaosong@princeton.edu
% First created date:   2020/09/30
% Last modified date:   2020/09/30

% compute the x loc and value for the bar plot normalized as the PDF estimate for the weighted data
[pdfxloc, ctval, probval, pdfval] = weighted_bar_pdf(vdist, wt, binw, edgetype);
binw = pdfxloc(2)-pdfxloc(1);
% xran = max(abs(pdfxloc(1)-binw/2-1), abs(pdfxloc(end)+binw/2)+1);   % show range of bar plot 

hold(ax,'on');
ax.FontSize = 9;
patarea = [0 0;
           -100 0;
           -100 1;
            0 1;
            0 0];
% patch(ax,patarea(:,1),patarea(:,2),'k','Facealpha',0.05,'edgecolor','none');
% plot(ax,[0 0],[-100 100],'--','color',[0.7 0.7 0.7]);

% plot the histogram normalized as pdf
barHdl = bar(ax,pdfxloc, pdfval,1,'stacked','facea',0.6);
barHdl(1).FaceColor = [1 96/255 0];
% plot(ax,[median(vdistsw) median(vdistsw)],[-100 100],'b--','linew',1.5);    % median of raw dist
% plot(ax,[mean(vdistsw) mean(vdistsw)],[-100 100],'b--','linew',1.5);    % mean of raw dist
% plot(ax,[median(tmp) median(tmp)],[-100 100],'b--','linew',1.5);    % median of nonzero-wt dist
% plot(ax,[mean(tmp) mean(tmp)],[-100 100],'b--','linew',1.5);    % mean of nonzero-wt dist
% med = wt_median(vdist,wt);
% plot(ax,[med med],[-100 100],'b--','linew',1.5);    % median of wt dist
% plot(ax,-50:0.05:50,pdffitsw,'r-','linew',1);     % pdf guassian fit curve
% plot(ax,[musw musw],[-100 100],'b-','linew',1.5);     % pdf guassian fit mean
% plot(ax,pdfxloc,pdfval,'b-','linew',1);     % pdf of data
ax.FontSize = 9;
% text(ax,0.04,0.93,'b','FontSize',11,'unit','normalized','EdgeColor','k','Margin',2);
ax.Box = 'on';
grid(ax, 'on');
ax.GridLineStyle = '--';
ylim(ax,[0 0.3]);
% xlim(ax,[-6 6]);
% xlim(ax, [-xran, xran]);
% xticks(ax,-4:1:4);
% yticks(ax,[0.05 0.1 0.15]);
xlabel(ax,'LF-HF (km)','fontsize',11);
ylabel(ax,'PDF estimate','fontsize',11);

% keyboard
