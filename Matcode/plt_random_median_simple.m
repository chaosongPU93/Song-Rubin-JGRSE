function [prob,f] = plt_random_median_simple(medrand,refval)
% This function is to plt the distribution of the medians of the randoms
% generated multiple times and compare to a desired value. Output the 
% percentage the abs(medrand) >= abs(refval) 
%
%   Served more as a simplification in the main code  
%
% By Chao Song, chaosong@princeton.edu
% First created date:   2020/10/01
% Last modified date:   2020/10/01

prob = zeros(1,size(medrand,2));

[scrsz, res] = pixelperinch(1);

f.fig=figure;
widin = 8;  % maximum width allowed is 8.5 inches
htin = 8;   % maximum height allowed is 11 inches
set(f.fig,'Position',[scrsz(1)+1*scrsz(3)/10 scrsz(2)+scrsz(4)/20 widin*res htin*res]);
nrow = size(medrand,2);
ncol = 2;
for isub = 1:nrow*ncol
    f.ax(isub) = subplot(nrow,ncol,isub);
end

titlestr = {'rand gen from custom pdf';
            'rand gen from ksdensity';
            'rand gen from normal fitting';
           };
       
% subplot 1,3,5
for i = 1: size(medrand,2)
    ax = f.ax(2*i-1);
    hold(ax,'on');
    ax.Box = 'on';
    grid(ax, 'on');
    N = length(medrand(:,i));
    scatter(ax,1:N, medrand(:,i),15,'k','x');
    plot(ax,[1 N],[refval refval],'b--','linew',1);
    plot(ax,[1 N],[-refval -refval],'b--','linew',1);
    prob(i) = sum(abs(medrand(:,i))>=abs(refval)) /N * 100;
    text(ax,0.1,0.9,strcat({'Prob(abs>=med) = '},sprintf('%.2f',prob(i)),{'%'}),'fontsize',12,...
         'unit','normalized');
    xlabel(ax,'Sampling index');
    ylabel(ax,'Mean of randoms');
    title(ax,titlestr(i, :), 'fontsize',12);
    hold(ax,'off');
end

% subplot 2,4,6
for i = 1: size(medrand,2)
    ax = f.ax(2*i);
    hold(ax,'on');
    ax.Box = 'on';
    grid(ax, 'on');
    histogram(ax,medrand(:,i),'facec',[0.6 0.6 0.6],'facea',0.8);
    plot(ax,[refval refval],ax.YLim,'b--','linew',1);
    plot(ax,[-refval -refval],ax.YLim,'b--','linew',1);
    xlabel(ax,'Mean of randoms');
    ylabel(ax,'Counts');
    hold(ax,'off');
end