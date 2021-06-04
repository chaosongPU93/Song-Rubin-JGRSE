function [pdfmisfit1, pdfmisfit2, fighdl] = compare_distribution2(randnum,muHat,sigmaHat,...
                                                                  vdist,wt,pdfxloc,pdfval)
% 
% This function is to get the distribution and corresponding pdf of a set of a numbers
% that are generated under a gaussian distribution, want to compare:
%   A. if it complies with parent gaussian
%   B. if it complies with parent data distribution that the gaussian is fitted from
%
%
% By Chao Song, chaosong@princeton.edu
% First created date:   2020/09/30
% Last modified date:   2020/09/30

% plot
[scrsz, res] = pixelperinch(1);

fighdl.fig = figure;
widin = 8;  % maximum width allowed is 8.5 inches
htin = 5;   % maximum height allowed is 11 inches
set(fighdl.fig,'Position',[scrsz(1)+1*scrsz(3)/10 scrsz(2)+scrsz(4)/20 widin*res htin*res]);

fighdl.ax(1) = subplot(121);
ax = gca;
hold(ax,'on');

% [ax, ~, ~, ~, pdfxloc, pdfval] = plt_weighted_dist(ax, vdist, wt, binw, edgetype);

% % trying to use the auto bin method histogram to have a sense of the reasonable bin width
% tmp = vdist(wt~=0);
% h=histogram(ax,tmp,'binmethod','auto','facec',[1 96/255 0],'facea',0.8);
% binw = h.BinWidth;
% delete(h);
% edgest = floor(min(vdist));
% edgeed = ceil(max(vdist));
% xloc = edgest+binw/2: binw: edgeed+binw/2;
% xloc = xloc';
% h = zeros(length(xloc),1);
% h2 = zeros(length(xloc),1);
% for i = 1: length(xloc)
% %     ibar = find(vdist>xloc(i)-binw/2 & vdist<=xloc(i)+binw/2);
%     wtobj = wt(vdist>xloc(i)-binw/2 & vdist<=xloc(i)+binw/2);
%     h(i) = sum(wtobj)/sum(wt);    % probability
%     h2(i) = h(i)/binw;  % probobility density estimate, pdf
% end
% indst = find(h2>0, 1, 'first');     % non-zero start index
% inded = find(h2>0, 1, 'last');      % non-zero end index
% xran1 = max(abs(floor(xloc(indst))), abs(ceil(xloc(inded))));   % show range of bar plot 
% % get the non-zero pdf and its corresponding x location
% pdfxloc = xloc(indst: inded);
% pdfval = h2(indst: inded);
% plot the histogram normalized as pdf
barHdl = bar(ax,pdfxloc, pdfval,1,'stacked','facea',0.8);
barHdl(1).FaceColor = [1 96/255 0];
plot(ax,pdfxloc,pdfval,'b-','linew',1.5);      % pdf of data 
ax.Box = 'on';
grid(ax, 'on');
ax.GridLineStyle = '--';
binw = pdfxloc(2)-pdfxloc(1);
xran = max(abs(pdfxloc(1)-binw/2), abs(pdfxloc(end)+binw/2));   % show range of bar plot 
xlim(ax, [-xran, xran]);
ylim(ax,[0 ax.YLim(2)*1.2]);
xx = -50:0.05:50;
pdffit = normpdf(xx, muHat, sigmaHat);
plot(ax,xx,pdffit,'k-','linew',1.5);     % pdf guassian fit curve
hold(ax,'off');
legend('Hist of data','PDF of data','PDF of fitted Gaussian');

fighdl.ax(2) = subplot(122);
ax = gca;
hold(ax,'on');
edge = pdfxloc(1)-binw/2: binw: pdfxloc(end)+binw/2;
hist1 = histogram(ax,randnum,'binedges',edge,'normalization','pdf','facec',[0.6 0.6 0.6],'facea',0.8);
plot(ax,pdfxloc,pdfval,'b-','linew',1.5);     % pdf of data
plot(ax,xx,pdffit,'k-','linew',1.5);     % pdf guassian fit curve
ylim(ax,[0 ax.YLim(2)*1.5]);
xlim(ax,fighdl.ax(1).XLim);
ax.Box = 'on';
grid(ax, 'on');
ax.GridLineStyle = '--';
plot(ax,[wt_mean(vdist, wt) wt_mean(vdist, wt)],ax.YLim,'b--','linew',1);     % pdf guassian fit curve
plot(ax,[mean(randnum) mean(randnum)],ax.YLim,'k--','linew',1);     % pdf guassian fit curve
disp(wt_mean(vdist, wt));
disp(mean(randnum));
disp(wt_mean(vdist, wt) - mean(randnum));
legend('Hist of generated numbers','PDF of data','PDF of fitted Gaussian','Mean of weighted data ',...
       'Mean of generated numbers');
hold off

%%% compute the misfit between data PDF and fitted Gaussian PDF
pdffit = normpdf(pdfxloc, muHat, sigmaHat);
pdfmisfit1 = sum((pdffit-pdfval).^2)/ length(pdfxloc);

%%% compute the misfit between data PDF and random numbers PDF 
% get the pdf and xloc of the random numbers from histogram
xrand = (hist1.BinEdges(1:end-1)+ hist1.BinEdges(2:end)) /2;
pdfrand = hist1.Values;

% interpolation of the randum numbers' pdf to the same x location as the data pdf
pdfi=interp1(xrand,pdfrand,pdfxloc,'linear');

% get the summed squares of the 2 pdfs as the misfits
nonnan = sum(~isnan(pdfi));
pdfmisfit2 = sum((pdfi(~isnan(pdfi)) - pdfval(~isnan(pdfi))).^2)/ nonnan;









