function [ax1, ax2] = plt_random_mean_simple(ax1, ax2, meanrand,refval,binw)
% This function is to plt the distribution of the medians of the randoms
% generated multiple times and compare to a desired value. Output the 
% percentage the abs(medrand) >= abs(refval) 
%
%   Served more as a simplification in the main code  
%
% By Chao Song, chaosong@princeton.edu
% First created date:   2020/10/01
% Last modified date:   2020/10/01

prob = zeros(1,size(meanrand,2));

% [scrsz, res] = pixelperinch(1);
% 
% f.fig=figure;
% widin = 8;  % maximum width allowed is 8.5 inches
% htin = 8;   % maximum height allowed is 11 inches
% set(f.fig,'Position',[scrsz(1)+1*scrsz(3)/10 scrsz(2)+scrsz(4)/20 widin*res htin*res]);
% nrow = 1;
% ncol = 2;
% for isub = 1:nrow*ncol
%     f.ax(isub) = subplot(nrow,ncol,isub);
% end

% titlestr = {'rand gen from normal fitting';
%             'rand gen from ksdensity';
%             'rand gen from custom pdf';
%            };

i=1;

% subplot 1
ax = ax1;
hold(ax,'on');
ax.Box = 'on';
grid(ax, 'on');
N = length(meanrand(:,i));
p1 = scatter(ax,1, meanrand(1,i),15,'k','+');
scatter(ax,1:N, meanrand(:,i),15,'k','+');
p2=plot(ax,[1 N],[refval(i) refval(i)],'b--','linew',1);
plot(ax,[1 N],[-refval(i) -refval(i)],'b--','linew',1);
prob(i) = sum(abs(meanrand(:,i))>=abs(refval(i))) /N ;
text(ax,0.6,0.95,strcat({'Prob(abs>=meandiff) = '},sprintf('%4.3e',prob(i))),'fontsize',10,...
    'unit','normalized','horizontalalignment','center');
xlabel(ax,'N\_th Sampling');
ylabel(ax,'Diff of mean of randoms');
% title(ax,titlestr(i, :), 'fontsize',12);
ylim(ax, [-1.5*abs(refval(i)) 1.5*abs(refval(i))]);
%     ylim(ax, [0 1.5*abs(refval)]);
legend([p1,p2],{'diff of mean of randoms','diff of mean of data'},...
    'location','south');
hold(ax,'off');

% subplot 2
mmeanrand = mean(meanrand,1);   % mean of diff of mean of random numbers
% binw = 0.1;
% following 2 lines could make the refval falls at the edge of the bin
nbin = ceil(2*refval/binw);
binw = 2*refval/nbin;

ax = ax2;
hold(ax,'on');
ax.Box = 'on';
grid(ax, 'on');
h=histogram(ax,meanrand(:,i),'binw',binw,'facec',[0.6 0.6 0.6],'facea',0.8);  %,'normalization','pdf'
p2 = plot(ax,[refval(i) refval(i)],[0 500],'b--','linew',1);
p3 = plot(ax,[mmeanrand(i) mmeanrand(i)],[0 200],'k--','linew',1);
plot(ax,[-refval(i) -refval(i)],[0 500],'b--','linew',1);
xlabel(ax,'Diff of mean of randoms');
ylabel(ax,'Counts');
xlim(ax, [-1.5*abs(refval(i)) 1.5*abs(refval(i))]);
ylim(ax, [0 max(h.Values)*5/4]);
xx = ax.XLim(1): 0.001: ax.XLim(2);
[muhat, sigmahat] = normfit(meanrand(:,i));
synpdf = normpdf(xx,muhat,sigmahat);    % synthetic PDF
synct = synpdf*binw*N;      % synthetic counts
plot(ax,xx, synct, 'r-', 'linew', 1.5);
%     xlim(ax, [0 1.5*abs(refval)]);
legend([p2,p3],{'diff of mean of data',...
    'mean of diff of mean of randoms'},'location','best');
hold(ax,'off');

% keyboard