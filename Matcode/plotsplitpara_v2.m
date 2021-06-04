%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% THIS is the script to plot the shear wave splitting parameter calculated
% by Chao (anyone)
%
%   Version 2, different from plotsplitpara, this script uses the lfe locations
%   inverted from hypoinverse, instead of from Bostocks's lfe fam
%
%   
%
% Modified by Chao Song, chaosong@princeton.edu
% First created date:   2019/09/02
% Last modified date:   2019/09/02
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
clc
close all

set(0,'DefaultFigureVisible','on');
scrsz=get(0,'ScreenSize');

% set path
workpath = getenv('ALLAN');
temppath = strcat(workpath, '/templates/');
if ~exist(temppath, 'dir')
    mkdir(temppath);
end

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
            '017';
            '006';
            '001';
            ];   
nfam = length(nfampool);

datapath = strcat(workpath,'/data-no-resp');
% split = load(strcat(datapath, '/split_chao/lzbtrio/','rot_para_11fam0.5-6.5_25_35LZB_40sps7stas'));

split = [];
lo = 0.5;
hi = 6.5;
sps = 40;
bef = 25;
aft = 35;
nsta = 7;
famp1 = ['002';'043';'141';'047';'010';'099';'017';'001';'019';'021';'045';'076';'176';'015';...
             '158';'231';'234';'006'];
famp2 = ['144';'068';'125';'147'];
for i = 1: nfam
    fam = nfampool(i,:);
    if ismember(fam,famp1,'rows')
        PREFIX = strcat(fam,'rot_para',num2str(lo),'-',num2str(hi),'_',num2str(bef),'_',...
            num2str(aft),'LZB');
    elseif ismember(fam,famp2,'rows')
        PREFIX = strcat(fam,'rot_para',num2str(lo),'-',num2str(hi),'_',num2str(bef),'_',...
            num2str(aft),'TWKB');
    end
    ROTS = load(strcat(datapath, '/split_chao/lzbtrio/',PREFIX,'_',num2str(sps),'sps',num2str(nsta),'stas'));
    split = [split; ROTS];
end    


%%% this is inverted from (0,0) of all fams, same order, location of control points
%%% inverted from /home/data2/chaosong/Seisbasics/hypoinverse/lzbrst/chat00p.reffam
contlzb = [-123.492667 48.451500 38.1400;
          -123.772167 48.493000 35.5900;
          -123.863167 48.528167 35.2100;
          -123.603333 48.440167 36.7100;
          -123.800167 48.408833 34.5200;
          -123.893333 48.536500 35.0700;
          -123.864500 48.498667 34.8800;
          -123.753333 48.525667 36.2000;
          -123.703667 48.502667 36.4100;
          -123.814333 48.538667 35.7900;
          -123.838500 48.544833 35.6600;
          -123.908000 48.494167 34.5100;       % 006
          -123.879667 48.446167 34.2600;       % 001
          ];

lfeloc = zeros(size(contlzb,1),2);
for i = 1: size(contlzb,1)
    [dx, dy] = absloc2relaloc(contlzb(i,1),contlzb(i,2),contlzb(2,1),contlzb(2,2));
    lfeloc(i,:) = [dx dy];
end

% station locations
dtdir = ('/home/data2/chaosong/matlab/allan/2004/JUL/');
flist= ('datalist');
outf = ('staloc.txt');
stainfo = getstainfo(dtdir, flist, outf);
ind=[3,11,4];
stainfo = stainfo(ind,:);

stalat = str2double(stainfo(:,2));
stalon = str2double(stainfo(:,3));
[staloc(:,1),staloc(:,2)] = absloc2relaloc(stalon,stalat,contlzb(2,1),contlzb(2,2));

% plot
[scrsz, res] = pixelperinch(1);
f.fig=figure;
f.fig.Renderer='Painters';
widin = 8;  % maximum width allowed is 8.5 inches
htin = 2.5;   % maximum height allowed is 11 inches
set(f.fig,'Position',[scrsz(1)+1*scrsz(3)/10 scrsz(2)+scrsz(4)/20 widin*res htin*res]);

nrow = 3;
ncol = 1;
for isub = 1:nrow*ncol
    f.ax(isub) = subplot(nrow,ncol,isub);
end

% reposition
set(f.ax(1), 'position', [ 0.1, 0.1, 0.25, 0.9]);
set(f.ax(2), 'position', [ 0.4, 0.1, 0.25, 0.9]);
set(f.ax(3), 'position', [ 0.7, 0.1, 0.25, 0.9]);

staind = [4,5,6];
for ii=1:3
    hold(f.ax(ii),'on');
    %%% plot the location of new fam pool
    scatter(f.ax(ii),lfeloc(:,1), lfeloc(:,2), 15,'filled','MarkerEdgeColor',[0.2 0.2 0.2],...
            'MarkerFaceColor',[0.2 0.2 0.2]);
    % text(lfeloc(ind,3),  lfeloc(ind,2), num2str(lfeloc(ind,1)), 'fontsize',12,'color','m');
    refsamp = 4;   % 0.1s in 40 sps
    for i = 1: nfam
        irow = (i-1)*7+staind(ii);
        angle = split(irow,2);
        offset = split(irow,1);
        x0 = lfeloc(i,1);
        y0 = lfeloc(i,2);
        reflen = 2.5;
        
        a = 0+1i*reflen*offset/refsamp;  % complex on positive y axis, reference vector points to N
        b = a*exp(1i * deg2rad(angle));   % rotate counter-clockwise by angle to the fast direction (axis)
        dx = real(b);
        dy = imag(b);
        
        x = [x0-0.5*dx; x0+0.5*dx];
        y = [y0-0.5*dy; y0+0.5*dy];
        
        plot(f.ax(ii),x,y,'color',[0.6 0.6 0.6],'linewidth',1.5);
        text(f.ax(ii),x0-1.5,y0+1.5,nfampool(i,:),'fontsize',7);
    end
    plot(f.ax(ii),[-10 -10+reflen],[19 19],'color',[0.6 0.6 0.6],'linewidth',1.5);
    text(f.ax(ii),-11,17.5,'0.1 s','fontsize',8);
    for jj=1:3
        if jj == 1
            plot(f.ax(ii),staloc(jj,1),staloc(jj,2),'^',...
                'MarkerSize',5,'MarkerEdgeColor','k');
        else
            plot(f.ax(ii),staloc(jj,1),staloc(jj,2),'s',...
                'MarkerSize',5,'MarkerEdgeColor','k');
        end
    end
    text(f.ax(ii),staloc(1,1)-1.5,staloc(1,2)+2,stainfo(1,1),...
        'fontsize',8);
    text(f.ax(ii),staloc(2,1)-2.5,staloc(2,2)+2,stainfo(2,1),...
        'fontsize',8);
    text(f.ax(ii),staloc(3,1)+1,staloc(3,2)+1,stainfo(3,1),...
        'fontsize',8);
    if ii == 1
        plot(f.ax(ii),staloc(ii,1),staloc(ii,2),'^','MarkerFaceColor','k',...
            'MarkerSize',5,'MarkerEdgeColor', 'k');
    else
        plot(f.ax(ii),staloc(ii,1),staloc(ii,2),'s','MarkerFaceColor','k',...
            'MarkerSize',5,'MarkerEdgeColor', 'k');
    end
    f.ax(ii).Box='on';
    grid(f.ax(ii),'on');
    f.ax(ii).GridLineStyle = '--';
    axis(f.ax(ii),'equal');
    axis(f.ax(ii),[-13 23 -13 22]);
    hold(f.ax(ii),'off');
end

% f.ax(2).YAxis.Visible = 'off';
% f.ax(3).YAxis.Visible = 'off';
xlabel(f.ax(1),'E (km)');
ylabel(f.ax(1),'N (km)');
xlabel(f.ax(2),'E (km)');
xlabel(f.ax(3),'E (km)');

print(f.fig,'-depsc2',strcat('/home/data2/chaosong/CurrentResearch/Song_Rubin_2020/Figures/splitpara.eps'));


% %%% how to plot a circle
% th = 0:pi/50:2*pi;
% xunit = reflen/6 * cos(th) + x0;
% yunit = reflen/6 * sin(th) + y0;
% plot(xunit, yunit, 'k','linewidth',1);
% 
% xunit = reflen/3 * cos(th) + x0;
% yunit = reflen/3 * sin(th) + y0;
% plot(xunit, yunit, 'k','linewidth',1.5);