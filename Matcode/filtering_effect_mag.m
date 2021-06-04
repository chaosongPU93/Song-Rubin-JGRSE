function filtering_effect_mag
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script is used to get the magnitude of the filtering correction in km
% for both PGC and LZB trio
%
% Chao Song, chaosong@princeton.edu
% First created date:   2019/12/05
% Last modified date:   2019/12/05
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Initialization
format short e   % Set the format to 5-digit floating point
clc
clear
close all

set(0,'DefaultFigureVisible','on');
% set(0,'DefaultFigureVisible','off');   % switch to show the plots or not

% set path
workpath = getenv('MHOME');
lzbpath = strcat(workpath, '/Seisbasics/hypoinverse/lzbrst');


%% PGC trio
%%% fc is inverted with (lf-hf_12, lf-hf_13) relative to fam 002
fcpgc = [-123.592000 48.436500 36.7900];

%%% this is inverted from (0,0) relative to fam 002, location of control points
contpgc = [-123.585000 48.436667 36.8800];

%%% convert to km
[dx, dy] = absloc2relaloc(fcpgc(1),fcpgc(2),contpgc(1),contpgc(2));
fcvecpgc = [dx dy];
fcmagpgc = sqrt(dx.^2+dy.^2);


%% LZB trio

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
%             '158';
%             '234';
            '006';
            '001';
            ];
nfam = length(nfampool);

%%% fc is inverted with (lf-hf_12, lf-hf_13) relative to all fams
%%% the same content as /home/data2/chaosong/Seisbasics/hypoinverse/lzbrst/LZBallfamfc, which is
%%% inverted by hypoinverse from /home/data2/chaosong/Seisbasics/hypoinverse/lzbrst/chat00p.fc
fclzb = [-123.509500 48.457000 38.0200;
    -123.801000 48.561333 36.2500;
    -123.896333 48.594000 35.7800;
    -123.638833 48.474000 36.7200;
    -123.797167 48.440333 34.8100;
    -123.925000 48.599333 35.5600;
    -123.898667 48.555833 35.2500;
    -123.772167 48.575333 36.7300;
    -123.734833 48.562667 36.8900;
    -123.837500 48.587500 36.3200;
    -123.867833 48.590000 36.0100;
%     -123.984500 48.498500 34.0300;       % 158
%     -123.974167 48.474333 33.9100;       % 234
    -123.930667 48.545167 34.8600;       % 006
    -123.892500 48.467333 34.3600;       % 001
    ];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% VERSION 1 for the locations of lfe families, by inversion from hypoinverse
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
%     -123.967000 48.458667 33.7800;       % 158
%     -123.960833 48.470000 33.9300;       % 234
    -123.908000 48.494167 34.5100;       % 006
    -123.879667 48.446167 34.2600;       % 001
    ];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%% VERSION 2 for the locations of lfe families, directly from Bostock's catalog
% bosdir = ('/home/data2/chaosong/matlab/allan/BOSTOCK');
% lfefnm = ('newlfeloc');
% contlzb = load(fullfile(bosdir, lfefnm));
% %%% lfeloc is loctions of all fams, but not necessarily each fam has detections
% [~,~,ind] = intersect(str2num(nfampool(:,:)), contlzb(:, 1),'stable');
% contlzb = contlzb(ind,:);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%% convert each filterring correction to km relative to its own fam
fcveclzb = zeros(size(fclzb,1),2);
fcmaglzb = zeros(size(fclzb,1),1);
for i = 1: size(fclzb,1)
    [dx, dy] = absloc2relaloc(fclzb(i,1),fclzb(i,2),contlzb(i,1),contlzb(i,2));
    fcveclzb(i,:) = [dx dy];
    fcmaglzb(i) = sqrt(dx.^2+dy.^2);
end

%%% convert lfe location to km relative to 043
loclfe = zeros(size(contlzb,1),2);
for i = 1: size(contlzb,1)
    [dx, dy] = absloc2relaloc(contlzb(i,1),contlzb(i,2),contlzb(2,1),contlzb(2,2));
    loclfe(i,:) = [dx dy];
end

dtdir = ('/home/data2/chaosong/matlab/allan/2004/JUL/');
flist= ('datalist');
outf = ('staloc.txt');
stainfo = getstainfo(dtdir, flist, outf);
ind=[3,11,4,5,8,6];
stainfo = stainfo(ind,:);

stalat = str2double(stainfo(:,2));
stalon = str2double(stainfo(:,3));
[staloc(:,1),staloc(:,2)] = absloc2relaloc(stalon,stalat,contlzb(2,1),contlzb(2,2));

%%% load the location resolution points, +-2 samples
reslzb = load(fullfile(lzbpath,'evtloc.13fam_locres_hf2spl'));
ifam = 2;   % choose fam 043 as the representative
reslzb = reslzb((ifam-1)*4+1: ifam*4, :);

% vecreslzb = zeros(size(reslzb,1),2);
[dx, dy] = absloc2relaloc(reslzb(:,1),reslzb(:,2),contlzb(ifam,1),contlzb(ifam,2));
vecreslzb = [dx dy];


%%
[scrsz, res] = pixelperinch(1);
f.fig=figure;
f.fig.Renderer='Painters';
widin = 6;  % maximum width allowed is 8.5 inches
htin = 5;   % maximum height allowed is 11 inches
set(f.fig,'Position',[scrsz(1)+1*scrsz(3)/10 scrsz(2)+scrsz(4)/20 widin*res htin*res]);

nrow = 1;
ncol = 2;
for isub = 1:nrow*ncol
    f.ax(isub) = subplot(nrow,ncol,isub);
end

% reposition
set(f.ax(1), 'position', [ 0.1, 0.1, 0.8, 0.8]);
% set(f.ax(2), 'position', [ 0.55, 0.1, 0.4, 0.88]);


ax = f.ax(1);

hold(ax,'on');
plot(ax,[-100 100],[0 0],'k--');
plot(ax,[0 0],[-100 100],'k--');
    
%%% plot the location of new fam pool
scatter(ax,loclfe(:,1), loclfe(:,2), 20,'filled','MarkerEdgeColor',[0.2 0.2 0.2],...
    'MarkerFaceColor',[0.2 0.2 0.2]);
% text fam name
text(ax,loclfe(:,1)-1.8, loclfe(:,2)+1.5,nfampool(:,:),'fontsize',8);
for jj=1:6
    if jj == 1
        plot(ax,staloc(jj,1),staloc(jj,2),'^',...
            'MarkerSize',7,'MarkerEdgeColor','k','markerfacec','k');
    elseif jj == 2 || jj == 3
        plot(ax,staloc(jj,1),staloc(jj,2),'s',...
            'MarkerSize',8,'MarkerEdgeColor','k','markerfacec','k');
    elseif jj == 5 || jj == 6
        plot(ax,staloc(jj,1),staloc(jj,2),'s',...
            'MarkerSize',8,'MarkerEdgeColor','k','markerfacec','r');
    else
        plot(ax,staloc(jj,1),staloc(jj,2),'^',...
            'MarkerSize',7,'MarkerEdgeColor','k','markerfacec','r');
    end
end

text(ax,staloc(1,1)-2,staloc(1,2)+2,stainfo(1,1),'fontsize',10);
text(ax,staloc(2,1)-3,staloc(2,2)+2,stainfo(2,1),'fontsize',10);
text(ax,staloc(3,1)+1,staloc(3,2)+1,stainfo(3,1),'fontsize',10);
text(ax,staloc(4,1)-2,staloc(4,2)+2,stainfo(4,1),'fontsize',10);
text(ax,staloc(5,1)-2.5,staloc(5,2)-2,stainfo(5,1),'fontsize',10);
text(ax,staloc(6,1)-2.5,staloc(6,2)+2,stainfo(6,1),'fontsize',10);

% filtering effect for lzb, invert ( (lf-hf)_12, (lf-hf)_13 ) to locations 
for i = 1: nfam
    x = [loclfe(i,1); loclfe(i,1)+fcveclzb(i,1)];
    y = [loclfe(i,2); loclfe(i,2)+fcveclzb(i,2)];
    plot(ax,x,y,'color',[0.6 0.6 0.6],'linewidth',2);
end

% filtering effect for pgc
x = [loclfe(1,1); loclfe(1,1)+fcvecpgc(1)];
y = [loclfe(1,2)-1; loclfe(1,2)+fcvecpgc(2)-1];
plot(ax,x,y,'color','b','linewidth',2);

% text(ax,loclfe(1,1), loclfe(1,2)-2,'002','fontsize',9);

% loc resolution indicated by +-2 sample cross
x = [10+ vecreslzb(1,1), 10+ vecreslzb(2,1)];
y = [25+ vecreslzb(1,2), 25+ vecreslzb(2,2)];
plot(ax,x,y,'color','r','linewidth',1.5);

x = [10+ vecreslzb(3,1), 10+ vecreslzb(4,1)];
y = [25+ vecreslzb(3,2), 25+ vecreslzb(4,2)];
plot(ax,x,y,'color','r','linewidth',1.5);
    
ax.Box='on';
grid(ax,'on');
ax.GridLineStyle = '--';
axis(ax,'equal');
% axis(ax,[-15 40 -13 32]);
axis(ax,[-15 55 -15 35]);
xlabel(ax,'E (km)');
ylabel(ax,'N (km)');

fcazilzb = zeros(size(fclzb,1),1);
for i = 1: nfam
    tmp = rad2deg(atan2(fcveclzb(i,1),fcveclzb(i,2)));
    if tmp < 0
        fcazilzb(i) = 360+tmp;
    else
        fcazilzb(i) = tmp;
    end
end

tmp = rad2deg(atan2(fcvecpgc(1),fcvecpgc(2)));
if tmp < 0
    fcazipgc = 360+tmp;
else
    fcazipgc = tmp;
end


% % move to a common center, excluding some different ones, i.e. 010, 047, 002, 001
% mind = [2,3,6,7,8,9,10,11,12];
% for i = 1: length(mind)
%     x = [30; 30+fcveclzb(mind(i),1)];
%     y = [0; 0+fcveclzb(mind(i),2)];
%     plot(ax,x,y,'color',[0.6 0.6 0.6],'linewidth',1);
% end
% maxazidiff = max(fcazilzb(mind)) - min(fcazilzb(mind));
% maxmagdiff = max(fcmaglzb(mind)) - min(fcmaglzb(mind));
% text(ax,25,-3,strcat({'mag diff: '}, sprintf('%.2f',maxmagdiff)));
% text(ax,25,-6,strcat({'azi diff: '}, sprintf('%.2f',maxazidiff)));

patarea = [22 -14;
           54 -14;
           54 10;
           22 10;
           22 -14;];
patch(ax,patarea(:,1),patarea(:,2),'w','edgecolor',[0.3 0.3 0.3]);
hold(ax,'off');


% a small inset
set(f.ax(2), 'position', [ 0.59, 0.23, 0.29, 0.25]);
ax = f.ax(2);
hold(ax,'on');
ax.Box='on';
% grid(ax,'on');
% ax.GridLineStyle = '--';
ax.FontSize = 7;

fcmaglzbprop = zeros(nfam,1);
fcmaglzbort = zeros(nfam,1);
loclfenew = zeros(nfam,2);
azi = 90-22.5;
for i = 1: nfam
    % decompose the absolute mag on its own azimuth to 2 components on ENE azi and its ortho
    fcmaglzbprop(i) = fcmaglzb(i)*cos(deg2rad(fcazilzb(i)-azi));
    fcmaglzbort(i) = fcmaglzb(i)*sin(deg2rad(fcazilzb(i)-azi));
    x0 = loclfe(i,1);
    y0 = loclfe(i,2);
    [loclfenew(i,1),loclfenew(i,2)] = coordinate_rot(x0,y0,-(azi-90),0,0);
end

famcheck = ['043';
            '068';
            '147';
            '141';
            '099';
            '006';
            '125';
            '017';
            '144';];
tmp1 = ['006';
        '144';
        '141';
        '147';
        '068';
        '125';];
       
for i = 1: size(famcheck,1)
    [~,ind] = ismember(famcheck(i,:),nfampool,'rows');
%     disp(fcazilzb(ind));    
    disp(fcmaglzbprop(ind));
    scatter(ax, loclfenew(ind,1), 0, 10,'filled','MarkerEdgeColor',[0.2 0.2 0.2],...
            'MarkerFaceColor',[0.2 0.2 0.2]);
    plot(ax,[loclfenew(ind,1) loclfenew(ind,1)], [0 fcmaglzbprop(ind)],'color',[0.6 0.6 0.6],...
         'linewidth',1.5);
    if ismember(nfampool(ind,:),tmp1,'rows')
        text(ax,loclfenew(ind,1), -0.08, nfampool(ind,:),'fontsize',6, 'horizontalalignment',...
             'center');
    else
        text(ax,loclfenew(ind,1), -0.16, nfampool(ind,:),'fontsize',6, 'horizontalalignment',...
             'center');
    end
end
axis(ax,[-10 6 -0.25 1]);
yticks(ax, -0.2: 0.2: 1);
xlabel(ax,'ENE (km)');
ylabel(ax,'Filtering effect along ENE (km)');


print(f.fig,'-depsc2',...
    strcat('/home/data2/chaosong/CurrentResearch/Song_Rubin_2020/Figure/filtercorrection.eps'));

print(f.fig,'-dpdf',...
    strcat('/home/data2/chaosong/CurrentResearch/Song_Rubin_2020/Figure/filtercorrection.pdf'));



keyboard











