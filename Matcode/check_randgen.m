function [pdfmisfit, fighdl] = check_randgen(randnum,vdist,wt,pdfxloc,pdfval,refval,binw)
% 
% This function is to get the distribution and corresponding pdf of a set of a numbers
% and compare it with the desired pdf
%
%   [pdfmisfit, fighdl] = compare_distribution1(randnum, pdfxloc, pdfval, binw)
%   pdfmisfit = compare_distribution1(randnum, pdfxloc, pdfval, binw)
%
% By Chao Song, chaosong@princeton.edu
% First created date:   2020/09/30
% Last modified date:   2020/09/30

% computing the cumulative distribution function for input pdf
cdfval = cumtrapz(pdfxloc,pdfval);

% and cut out the parts
ind=[true; not(diff(cdfval)==0)];
cdfval=cdfval(ind);
pdfxloc=pdfxloc(ind);
pdfval=pdfval(ind);

% plot    
[scrsz, res] = pixelperinch(1);

fighdl.fig = figure;
widin = 8;  % maximum width allowed is 8.5 inches
htin = 5;   % maximum height allowed is 11 inches
set(fighdl.fig,'Position',[scrsz(1)+1*scrsz(3)/10 scrsz(2)+scrsz(4)/20 widin*res htin*res]);

fighdl.ax(1) = subplot(3,4,[1 2 5 6]);
hold on
edge = pdfxloc(1)-binw/2: binw: pdfxloc(end)+binw/2;
hist1 = histogram(randnum,'binedges',edge,'normalization','pdf','facec',[0 1 1],'facea',0.8);
plot(pdfxloc, pdfval./trapz(pdfxloc,pdfval),'r', 'linew', 1.5);
ax = gca;
ylim([0 ax.YLim(2)*1.5]);
grid on
plot(ax,[refval refval],ax.YLim,'r--','linew',1);     % expectation of pdf
plot(ax,[mean(randnum) mean(randnum)],ax.YLim,'k--','linew',1);     % pdf guassian fit curve
% disp(wt_mean(vdist, wt));
% disp(mean(randnum));
disp(refval - mean(randnum));

legend('PDF of generated numbers','Input PDF of data','Mean of weighted data ',...
       'Mean of generated numbers');
hold off

fighdl.ax(2) = subplot(3,4,[3 4 7 8]);
plot(pdfxloc, cdfval,'k', 'linew', 2);
ylim([0 1])
grid on
legend('Input cdf of data');

fighdl.ax(3) = subplot(3,4,[9:12]);
plot(randnum)
grid on
legend('Generated numbers');
    
% get the pdf and xloc of the random numbers from histogram
xrand = (hist1.BinEdges(1:end-1)+ hist1.BinEdges(2:end)) /2;
pdfrand = hist1.Values;

% interpolation of the randum numbers' pdf to the same x location as the data pdf
pdfi=interp1(xrand,pdfrand,pdfxloc,'linear');

% get the summed squares of the 2 pdfs as the misfits
nonnan = sum(~isnan(pdfi));
pdfmisfit = sum((pdfi(~isnan(pdfi)) - pdfval(~isnan(pdfi))).^2)/ nonnan;

% output graphs if 2 outputs besides misfit needed 
if nargout==1       % output misfit only
    close(fighdl.fig)
end


% keyboard











