function [ax,c] = plt_mig_prop_onmap(ax, indne, indsw, ranvechf98, medallhf, angbest, xran,yran)
% This function to plot the migrations represented by its mdedian and propagation arrow
% bounded by its range on the map
%   Served more as a simplification in the main code  
%
% By Chao Song, chaosong@princeton.edu
% First created date:   2020/09/30
% Last modified date:   2020/09/30

hold(ax,'on');
plot(ax,[-100 100],[0 0],'k--');
plot(ax,[0 0],[-100 100],'k--');
ax.FontSize = 9;
axis(ax, 'equal');
for i = 1: length(indsw)
    [rotx1, roty1] = complex_rot(0, ranvechf98(indsw(i),2)-medallhf(indsw(i),3), ...
                                 -angbest(indsw(i)));
    [rotx2, roty2] = complex_rot(0, medallhf(indsw(i),3)-ranvechf98(indsw(i),1), ...
                                 -angbest(indsw(i)));
    xvect = [medallhf(indsw(i),1)-rotx2 medallhf(indsw(i),1)+rotx1];
    yvect = [medallhf(indsw(i),2)-roty2 medallhf(indsw(i),2)+roty1];
    drawArrow(ax,xvect,yvect,xran,yran,'linewidth',1,'linestyle','-','color',[0.7 0.7 0.7]);
%     arrow = annotation('arrow','linewidth',1,'linestyle','-','color',[0.7 0.7 0.7],'headwidth',8,...
%                        'headlength',8);
%     arrow.Parent = ax;
%     arrow.Position = [medallhf(indsw(i),1)-rotx2, medallhf(indsw(i),2)-roty2, ...
%                       rotx1+rotx2, roty1+roty2] ;
%     text(ax,0.5,0.1,num2str(indsw(i)),'FontSize',12,'unit','normalized');
end
for i = 1: length(indsw)
    scatter(ax,medallhf(indsw(i),1),medallhf(indsw(i),2), 30, indsw(i), 'filled','o');  %, 'MarkerEdgeColor', 'w')
end

for i = 1: length(indne)
    [rotx1, roty1] = complex_rot(0, ranvechf98(indne(i),2)-medallhf(indne(i),3), ...
                                 -angbest(indne(i)));
    [rotx2, roty2] = complex_rot(0, medallhf(indne(i),3)-ranvechf98(indne(i),1), ...
                                 -angbest(indne(i)));
    xvect = [medallhf(indne(i),1)-rotx2 medallhf(indne(i),1)+rotx1];
    yvect = [medallhf(indne(i),2)-roty2 medallhf(indne(i),2)+roty1];
    drawArrow(ax,xvect,yvect,xran,yran,'linewidth',1,'linestyle','-','color','k');
end
for i = 1: length(indne)
    scatter(ax,medallhf(indne(i),1),medallhf(indne(i),2), 30, indne(i), 'filled','o');  %, 'MarkerEdgeColor', 'w')
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
% hold(ax,'off');