% function identify_tremor_bursts_intertime
% This scripts is to generate time-sample plots and time-distance plots
% to visually identify the time periods of tremor bursts that are 
% potentially eligible RTMs which would then undergo regression analysis
% for further discards
%       
%
%
% By Chao Song, chaosong@princeton.edu
% First created date:   2020/12/23
% Last modified date:   2020/12/23

%% Initialization
format short e   % Set the format to 5-digit floating point
clear
close all
clc

% set(0,'DefaultFigureVisible','on');
set(0,'DefaultFigureVisible','off');   % switch to show the plots or not

% get the scrsz in pixels and number of pixels per inch of monitor 1
[scrsz, res] = pixelperinch(1);

workpath = getenv('MHOME');
rstpath = strcat(workpath, '/Seisbasics/hypoinverse/lzbrst');

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
%             '158';      % 158, 20200916,testing purpose
            ];     % 2020/09/07, i am adding a new family 006 to the pool, final version
        
nfam = size(nfampool,1);
disp(nfam); 

% load files
winlenhf = 4;

winlenlf = 16;
lofflf = 4;
ccminlf = 0.35;

SUFFIXhf = strcat('up.hf.time.',num2str(winlenhf),'_',num2str(winlenlf),'.',num2str(lofflf),...
                  '.',num2str(ccminlf));
fname = strcat(rstpath, '/evtloc.all',num2str(nfam),'fam.dcutnodou.',SUFFIXhf);
hftime = load(fname);
% 25 cols, format is:
%   E(043) N(043) E(own) N(own) dist(own) lon lat dep off12 off13 off12sec off13sec date
%   main_arrival_time  cent_of_win  avecc ccadd1 ccadd2 ccadd3 ccadd4 ampmax ampmin eng normeng famnum
%


SUFFIXlf = strcat('lf.time.',num2str(winlenhf),'_',num2str(winlenlf),'.',num2str(lofflf),...
                  '.',num2str(ccminlf));
fname = strcat(rstpath, '/evtloc.all',num2str(nfam),'fam.dcutnodou.',SUFFIXlf);
lftime = load(fname);

% this is inverted from (0,0) of all fams, same order, location of control points
loccont = [-123.492667 48.451500 38.1400; 
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
%            -123.967000 48.458667 33.7800;       % 158, 20200916,testing purpose
           ];
       

relacont = loccont;
[dx, dy] = absloc2relaloc(loccont(:,1),loccont(:,2),loccont(2,1),loccont(2,2));
relacont(:,1) = dx;
relacont(:,2) = dy;


% sort according to day, sec
hftime = sortrows(hftime, [13, 15]);
lftime = sortrows(lftime, [13, 15]);


% for 2005 only, rotate to the strike N30W and dip N60E direction, 2003 and 2004 have the main front
% migration direction to N
rotang = 30;  % counter-clockwise from y to strike, 2005

% obtain the relative time in days
ih03 = find(hftime(:,13) < 2004*1000);
hftime(ih03,26) = (hftime(ih03,13)-2003060)+hftime(ih03,15)./(3600.*24);
ih04 = find(hftime(:,13) < 2005*1000 & hftime(:,13) > 2004*1000);
hftime(ih04,26) = (hftime(ih04,13)-2004194)+hftime(ih04,15)./(3600.*24);
ih05 = find(hftime(:,13) > 2005*1000);
hftime(ih05,26) = (hftime(ih05,13)-2005254)+hftime(ih05,15)./(3600.*24);
% rotate to the strike N30W and dip N60E direction,now the 1,2 col is changing from E(043) and N(043),
% to down-dip and strike
[hftime(ih05,1),hftime(ih05,2)] = coordinate_rot(hftime(ih05,1),hftime(ih05,2),rotang,0,0);

il03 = find(lftime(:,13) < 2004*1000);
lftime(il03,26) = (lftime(il03,13)-2003060)+lftime(il03,15)./(3600.*24);
il04 = find(lftime(:,13) < 2005*1000 & lftime(:,13) > 2004*1000);
lftime(il04,26) = (lftime(il04,13)-2004194)+lftime(il04,15)./(3600.*24);
il05 = find(lftime(:,13) > 2005*1000);
lftime(il05,26) = (lftime(il05,13)-2005254)+lftime(il05,15)./(3600.*24);
% rotate to the strike N30W and dip N60E direction,now the 1,2 col is changing from E(043) and N(043),
% to down-dip and strike
[lftime(il05,1),lftime(il05,2)] = coordinate_rot(lftime(il05,1),lftime(il05,2),rotang,0,0);


%%
% obtain the separation in time between itself and its nearest
% detection in time for the catalog of each ETS.
hfnear03 = zeros(size(hftime(ih03,:),1),2);
% 1st col: occurence time
% 2nd col: inter-detection time to its nearst detection
hfnear03(:,1) = hftime(ih03,26);
for i = 1: size(hfnear03,1)
    if i == 1
        hfnear03(i,2) = hfnear03(2,1)-hfnear03(1,1);
    elseif i == size(hfnear03,1)
        hfnear03(i,2) = hfnear03(end,1)-hfnear03(end-1,1);
    else
        hfnear03(i,2) = min(hfnear03(i+1,1)-hfnear03(i,1), hfnear03(i,1)-hfnear03(i-1,1));
    end
end

hfnear04 = zeros(size(hftime(ih04,:),1),2);
hfnear04(:,1) = hftime(ih04,26);
for i = 1: size(hfnear04,1)
    if i == 1
        hfnear04(i,2) = hfnear04(2,1)-hfnear04(1,1);
    elseif i == size(hfnear04,1)
        hfnear04(i,2) = hfnear04(end,1)-hfnear04(end-1,1);
    else
        hfnear04(i,2) = min(hfnear04(i+1,1)-hfnear04(i,1), hfnear04(i,1)-hfnear04(i-1,1));
    end
end

hfnear05 = zeros(size(hftime(ih05,:),1),2);
hfnear05(:,1) = hftime(ih05,26);
for i = 1: size(hfnear05,1)
    if i == 1
        hfnear05(i,2) = hfnear05(2,1)-hfnear05(1,1);
    elseif i == size(hfnear05,1)
        hfnear05(i,2) = hfnear05(end,1)-hfnear05(end-1,1);
    else
        hfnear05(i,2) = min(hfnear05(i+1,1)-hfnear05(i,1), hfnear05(i,1)-hfnear05(i-1,1));
    end
end

% throw way 5% the catalog of each ETS, whose nearst time is highest. 
ind0395 = find(hfnear03(:,2)<=prctile(hfnear03(:,2),95));
hfnear0395 = hfnear03(ind0395, :);
ind0495 = find(hfnear04(:,2)<=prctile(hfnear04(:,2),95));
hfnear0495 = hfnear04(ind0495, :);
ind0595 = find(hfnear05(:,2)<=prctile(hfnear05(:,2),95));
hfnear0595 = hfnear05(ind0595, :);

% Or 5% of the entire catalog?
hfnear = [hfnear03; hfnear04; hfnear05];
ind0395a = find(hfnear03(:,2)<=prctile(hfnear(:,2),95));
hfnear0395a = hfnear03(ind0395a, :);
ind0495a = find(hfnear04(:,2)<=prctile(hfnear(:,2),95));
hfnear0495a = hfnear04(ind0495a, :);
ind0595a = find(hfnear05(:,2)<=prctile(hfnear(:,2),95));
hfnear0595a = hfnear05(ind0595a, :);

%% 
f.fig = figure;
f.fig.Renderer = 'painters';
widin = 12;  % maximum width allowed is 8.5 inches
htin = 10;   % maximum height allowed is 11 inches
set(f.fig,'Position',[1*scrsz(3)/20 scrsz(4)/10 widin*res htin*res]);

nrow = 3;
ncol = 1;

for isub = 1:nrow*ncol
    f.ax(isub) = subplot(nrow,ncol,isub);
    f.ax(isub).Box = 'on';
end

msizehf = 2;

ax = f.ax(1);
hold(ax,'on'); 
scatter(ax, hfnear03(:,1), hfnear03(:,2), msizehf, 'filled','ro');  %, 'MarkerEdgeColor', 'w')
plot(ax, ax.XLim, [prctile(hfnear03(:,2),95) prctile(hfnear03(:,2),95)], 'k--');
plot(ax, ax.XLim, [prctile(hfnear(:,2),95) prctile(hfnear(:,2),95)], 'b--');
plot(ax, ax.XLim, [1e-3 1e-3], 'g--');
text(ax,0.9,0.9,'2003','FontSize',10,'unit','normalized','horizontalalignment','left',...
    'EdgeColor','k','Margin',2);
ylabel(ax, 'Time (day) to neastest detection');
xlabel(ax, 'Time (day) since 2003060 Mar. 1, 2003');
ylim(ax,[0, 0.005]);
hold(ax,'off');

ax = f.ax(2);
hold(ax,'on'); 
scatter(ax, hfnear04(:,1), hfnear04(:,2), msizehf, 'filled','ro');  %, 'MarkerEdgeColor', 'w')
plot(ax, ax.XLim, [prctile(hfnear04(:,2),95) prctile(hfnear04(:,2),95)], 'k--');
plot(ax, ax.XLim, [prctile(hfnear(:,2),95) prctile(hfnear(:,2),95)], 'b--');
plot(ax, ax.XLim, [1e-3 1e-3], 'g--');
text(ax,0.9,0.9,'2004','FontSize',10,'unit','normalized','horizontalalignment','left',...
    'EdgeColor','k','Margin',2);
ylabel(ax, 'Time (day) to neastest detection');
xlabel(ax, 'Time (day) since 2004194 Jul. 12, 2004');
ylim(ax,[0, 0.005]);
hold(ax,'off');

ax = f.ax(3);
hold(ax,'on'); 
scatter(ax, hfnear05(:,1), hfnear05(:,2), msizehf, 'filled','ro');  %, 'MarkerEdgeColor', 'w')
plot(ax, ax.XLim, [prctile(hfnear05(:,2),95) prctile(hfnear05(:,2),95)], 'k--');
plot(ax, ax.XLim, [prctile(hfnear(:,2),95) prctile(hfnear(:,2),95)], 'b--');
plot(ax, ax.XLim, [1e-3 1e-3], 'g--');
text(ax,0.9,0.9,'2005','FontSize',10,'unit','normalized','horizontalalignment','left',...
    'EdgeColor','k','Margin',2);
ylabel(ax, 'Time (day) to neastest detection');
xlabel(ax, 'Time (day) since 2005254 Sep. 11, 2005');
ylim(ax,[0, 0.005]);
hold(ax,'off');

print(f.fig,'-depsc2',strcat(rstpath,'/nearest_time.eps'));  
    

%%
hfnear0395a = hfnear03;
hfnear0495a = hfnear04;
hfnear0595a = hfnear05;

% for the remaining detections, obtain the separation in time between itself and its preceding
% detection for the catalog of each ETS.
hfinter03 = zeros(size(hfnear0395a,1),2);
% 1st col: occurence time
% 2nd col: inter-detection time to its preceding detection
hfinter03(:,1) = hfnear0395a(:,1);
for i = 1: size(hfinter03,1)
    if i == 1
        hfinter03(i,2) = 0;
    elseif i == size(hfinter03,1)
        hfinter03(i,2) = hfinter03(end,1)-hfinter03(end-1,1);
    else
        hfinter03(i,2) = hfinter03(i,1)-hfinter03(i-1,1);
    end
end

hfinter04 = zeros(size(hfnear0495a,1),2);
hfinter04(:,1) = hfnear0495a(:,1);
for i = 1: size(hfinter04,1)
    if i == 1
        hfinter04(i,2) = 0;
    elseif i == size(hfinter04,1)
        hfinter04(i,2) = hfinter04(end,1)-hfinter04(end-1,1);
    else
        hfinter04(i,2) = hfinter04(i,1)-hfinter04(i-1,1);
    end
end

hfinter05 = zeros(size(hfnear0595a,1),2);
hfinter05(:,1) = hfnear0595a(:,1);
for i = 1: size(hfinter05,1)
    if i == 1
        hfinter05(i,2) = 0;
    elseif i == size(hfinter05,1)
        hfinter05(i,2) = hfinter05(end,1)-hfinter05(end-1,1);
    else
        hfinter05(i,2) = hfinter05(i,1)-hfinter05(i-1,1);
    end
end

hfinter = [hfinter03; hfinter04; hfinter05];

%% 
f.fig = figure;
f.fig.Renderer = 'painters';
widin = 8;  % maximum width allowed is 8.5 inches
htin = 8;   % maximum height allowed is 11 inches
set(f.fig,'Position',[1*scrsz(3)/20 scrsz(4)/10 widin*res htin*res]);

nrow = 3;
ncol = 1;

for isub = 1:nrow*ncol
    f.ax(isub) = subplot(nrow,ncol,isub);
    f.ax(isub).Box = 'on';
end

msizehf = 2;

% set a threshold of inter-detection time
ttol1 = 2e-3;
ttol2 = 1e-3;

ax = f.ax(1);
hold(ax,'on'); 
scatter(ax, hfinter03(:,1), hfinter03(:,2), msizehf, 'filled','ro');  %, 'MarkerEdgeColor', 'w')
% plot(ax, ax.XLim, [ttol1 ttol1], 'k--');
plot(ax, ax.XLim, [ttol2 ttol2], 'b--');
text(ax,0.92,0.9,'2003','FontSize',11,'unit','normalized','horizontalalignment','left',...
    'EdgeColor','k','Margin',2);
ylabel(ax, 'Time (day) to the preceding detection');
xlabel(ax, 'Time (day) since Mar. 1, 2003 (2003060)');
ylim(ax,[0, 0.002]);
hold(ax,'off');

ax = f.ax(2);
hold(ax,'on'); 
scatter(ax, hfinter04(:,1), hfinter04(:,2), msizehf, 'filled','ro');  %, 'MarkerEdgeColor', 'w')
% plot(ax, ax.XLim, [ttol1 ttol1], 'k--');
plot(ax, ax.XLim, [ttol2 ttol2], 'b--');
text(ax,0.92,0.9,'2004','FontSize',11,'unit','normalized','horizontalalignment','left',...
    'EdgeColor','k','Margin',2);
% ylabel(ax, 'Time (day) to the preceding detection');
xlabel(ax, 'Time (day) since Jul. 12, 2004 (2004194)');
ylim(ax,[0, 0.002]);
hold(ax,'off');

ax = f.ax(3);
hold(ax,'on'); 
scatter(ax, hfinter05(:,1), hfinter05(:,2), msizehf, 'filled','ro');  %, 'MarkerEdgeColor', 'w')
% plot(ax, ax.XLim, [ttol1 ttol1], 'k--');
plot(ax, ax.XLim, [ttol2 ttol2], 'b--');
text(ax,0.92,0.9,'2005','FontSize',11,'unit','normalized','horizontalalignment','left',...
    'EdgeColor','k','Margin',2);
% ylabel(ax, 'Time (day) to the preceding detection');
xlabel(ax, 'Time (day) since Sep. 11, 2005 (2005254)');
ylim(ax,[0, 0.002]);
hold(ax,'off');

print(f.fig,'-dpdf',strcat(rstpath,'/inter-detection_time.pdf'));  

%%
% set a threshold of inter-detection time
ttol1 = mean(hfinter(:,2));
ttol2 = 5e-4;
ttol3 = 1e-3;
ttol4 = 2e-3;
ttol5 = 3e-3;
ttol = ttol3;

% set a threshold of minimum number of detections in the burst period
% ntol = 15;
ntol = 20;


% group the indice of the tremor bursts accroding to the threshold on the inter-detection time and
% number of detections in the burst for both HF and LF
[burst03, nburstlf03] = group_tremor_burst(hfinter03,lftime(il03,26),ttol,ntol);
[burst04, nburstlf04] = group_tremor_burst(hfinter04,lftime(il04,26),ttol,ntol);
[burst05, nburstlf05] = group_tremor_burst(hfinter05,lftime(il05,26),ttol,ntol);

% recover the index of the detections to occurence times with the same format as trange
if ~isempty(burst03)
    newt03 = zeros(size(burst03,1), 2);
    ntran03 = zeros(size(burst03,1), 3);
    nhf03 = 0;
    for i = 1: size(burst03,1)
        nhf03 = nhf03+size(burst03{i},1);
        ind = burst03{i};
        newt03(i,1) = hfinter03(ind(1),1); %-0.5/60/24 
        newt03(i,2) = hfinter03(ind(end),1); %+0.5/60/24   
    end
    ntran03(:,1) = floor(newt03(:,1))+2003060;     
    ntran03(:,2:3) = (newt03(:,1:2)-floor(newt03(:,1))).*(3600.*24);
    perchf03 = nhf03/length(ih03)*100;
    nlf03 = sum(nburstlf03);
    perclf03 = nlf03/length(il03)*100;
end

if ~isempty(burst04)
    newt04 = zeros(size(burst04,1), 2);
    ntran04 = zeros(size(burst04,1), 3);
    nhf04 = 0;
    for i = 1: size(burst04,1)
        nhf04 = nhf04+size(burst04{i},1);
        ind = burst04{i};
        newt04(i,1) = hfinter04(ind(1),1);  %-0.5/60/24 add half minute as a taper
        newt04(i,2) = hfinter04(ind(end),1); % +0.5/60/24    
    end
    ntran04(:,1) = floor(newt04(:,1))+2004194;     
    ntran04(:,2:3) = (newt04(:,1:2)-floor(newt04(:,1))).*(3600.*24);
    perchf04 = nhf04/length(ih04)*100;
    nlf04 = sum(nburstlf04);
    perclf04 = nlf04/length(il04)*100;
end

if ~isempty(burst05)
    newt05 = zeros(size(burst05,1), 2);
    ntran05 = zeros(size(burst05,1), 3);
    nhf05 = 0;
    for i = 1: size(burst05,1)
        nhf05 = nhf05+size(burst05{i},1);
        ind = burst05{i};
        newt05(i,1) = hfinter05(ind(1),1);  %-0.5/60/24 add half minute as a taper
        newt05(i,2) = hfinter05(ind(end),1); %+0.5/60/24   
    end
    ntran05(:,1) = floor(newt05(:,1))+2005254;     
    ntran05(:,2:3) = (newt05(:,1:2)-floor(newt05(:,1))).*(3600.*24);
    perchf05 = nhf05/length(ih05)*100;
    nlf05 = sum(nburstlf05);
    perclf05 = nlf05/length(il05)*100;
end

perchf = (nhf03+nhf04+nhf05)/(length(ih03)+length(ih04)+length(ih05))*100;
perclf = (nlf03+nlf04+nlf05)/(length(il03)+length(il04)+length(il05))*100;


%%
f.fig = figure;
f.fig.Renderer = 'painters';
widin = 8;  % maximum width allowed is 8.5 inches
htin = 8;   % maximum height allowed is 11 inches
set(f.fig,'Position',[1*scrsz(3)/20 scrsz(4)/10 widin*res htin*res]);

nrow = 3;
ncol = 1;

for isub = 1:nrow*ncol
    f.ax(isub) = subplot(nrow,ncol,isub);
    f.ax(isub).Box = 'on';
end

msizehf = 2;

ax = f.ax(1);
hold(ax,'on');
if ~isempty(burst03)
    for j = 1: size(newt03,1)
        patarea = [newt03(j,1) -3e+4;
                   newt03(j,2) -3e+4;
                   newt03(j,2) 3e+4;
                   newt03(j,1) 3e+4;
                   newt03(j,1) -3e+4];
        patch(ax,patarea(:,1),patarea(:,2),'k','Facealpha',0.3,'edgecolor','none');
    end
end
stairs(ax, hfinter03(:,1), 1: length(hfinter03(:,1)),'r-');  %, 'MarkerEdgeColor', 'w')
stairs(ax, hftime(ih03,26), 1: length(hftime(ih03,26)),'b-');  %, 'MarkerEdgeColor', 'w')

text(ax,0.92,0.9,'2003','FontSize',11,'unit','normalized','horizontalalignment','left',...
    'EdgeColor','k','Margin',2);
text(ax,0.06,0.8,strcat(num2str(round(perchf03)),'%'),'FontSize',12,'unit','normalized','horizontalalignment',...
    'left');
ylabel(ax, 'Cumulative number of detections');
xlabel(ax, 'Time (day) since Mar. 1, 2003 (2003060)');
ylim(ax,[0, 6e+3]);
% ytickformat(ax,'%e');
hold(ax,'off');

ax = f.ax(2);
hold(ax,'on'); 
if ~isempty(burst04)
    for j = 1: size(newt04,1)
        patarea = [newt04(j,1) -3e+4;
                   newt04(j,2) -3e+4;
                   newt04(j,2) 3e+4;
                   newt04(j,1) 3e+4;
                   newt04(j,1) -3e+4];
        patch(ax,patarea(:,1),patarea(:,2),'k','Facealpha',0.3,'edgecolor','none');
    end
end
stairs(ax, hfinter04(:,1), (1: length(hfinter04(:,1)))','r-');  %, 'MarkerEdgeColor', 'w')
stairs(ax, hftime(ih04,26), 1: length(hftime(ih04,26)),'b-');  %, 'MarkerEdgeColor', 'w')
text(ax,0.92,0.9,'2004','FontSize',11,'unit','normalized','horizontalalignment','left',...
    'EdgeColor','k','Margin',2);
text(ax,0.06,0.8,strcat(num2str(round(perchf04)),'%'),'FontSize',12,'unit','normalized','horizontalalignment',...
    'left');
% ylabel(ax, 'Cumulative number of detections');
xlabel(ax, 'Time (day) since Jul. 12, 2004 (2004194)');
ylim(ax,[0, 2e+4]);
hold(ax,'off');

ax = f.ax(3);
hold(ax,'on'); 
if ~isempty(burst05)
    for j = 1: size(newt05,1)
        patarea = [newt05(j,1) -3e+4;
                   newt05(j,2) -3e+4;
                   newt05(j,2) 3e+4;
                   newt05(j,1) 3e+4;
                   newt05(j,1) -3e+4];
        patch(ax,patarea(:,1),patarea(:,2),'k','Facealpha',0.3,'edgecolor','none');
    end
end
stairs(ax, hfinter05(:,1), (1: length(hfinter05(:,1)))','r-');  %, 'MarkerEdgeColor', 'w')
stairs(ax, hftime(ih05,26), (1: length(hftime(ih05,26)))','b-');  %, 'MarkerEdgeColor', 'w')
text(ax,0.92,0.9,'2005','FontSize',11,'unit','normalized','horizontalalignment','left',...
    'EdgeColor','k','Margin',2);
text(ax,0.06,0.8,strcat(num2str(round(perchf05)),'%'),'FontSize',12,'unit','normalized','horizontalalignment',...
    'left');
% ylabel(ax, 'Cumulative number of detections');
xlabel(ax, 'Time (day) since Sep. 11, 2005 (2005254)');
ylim(ax,[0, 3e+4]);
hold(ax,'off');

print(f.fig,'-dpdf',strcat(rstpath,'/selected_bursts.pdf'));

% N = 10000;
% burst05 = cell(1,N);
% j = 1;
% % for j = 1: N
% tmp = [];
% for i = 1: length(hfinter04(:,2))
%     if hfinter04(i,2) <= ttol
%         tmp = [tmp; i];
%         
%     else
%         if length(tmp) >= ntol
%             burst05{j} = tmp;
%             j = j + 1;
%         end
%         tmp = [];
%         continue
%     end
%     if j > N
%         disp('array size is too small');
%         break
%     end
% end

%%
trange = ntran05;

angbestl1 = zeros(size(trange,1),1);
angbestl2 = zeros(size(trange,1),1);
angbestl3 = zeros(size(trange,1),1);
angbestl4 = zeros(size(trange,1),1);
angbestl5 = zeros(size(trange,1),1);

xran = [-15 25];
yran = [-20 20];

resprophf = nan(size(trange,1)+1,200);
resproplf = nan(size(trange,1)+1,50);

% create fit object with free constraints
fttpfree = fittype( @(a,b,x) a*x+b);

indcheck = 1: size(trange,1);
indplt = indcheck;

for i = 1: length(indplt)
% for i = 41: 63
    disp(trange(indplt(i),:));
    indhf = find(hftime(:,13)==trange(indplt(i),1) & hftime(:,15)>=trange(indplt(i),2) & ...
                 hftime(:,15)<=trange(indplt(i),3));
    mighf = hftime(indhf,:);

    indlf = find(lftime(:,13)==trange(indplt(i),1) & lftime(:,15)>=trange(indplt(i),2) & ...
                 lftime(:,15)<=trange(indplt(i),3));
    miglf = lftime(indlf,:);
    
    if size(mighf,1)< ntol || size(miglf,1) < ntol
        disp(strcat(num2str(indplt(i)),' not enough detections'));
    end
% end
    angle = 0:5:360;
    
%     l1normhf = zeros(length(angle),1);
%     l2normhf = zeros(length(angle),1);
    slopehf = zeros(length(angle),1);
%     ssehf = zeros(length(angle),1);
    rmsehf = zeros(length(angle),1);
    rsquarehf = zeros(length(angle),1);
    
%     l1normlf = zeros(length(angle),1);
%     l2normlf = zeros(length(angle),1);
    slopelf = zeros(length(angle),1);
%     sself = zeros(length(angle),1);
    rmself = zeros(length(angle),1);
    rsquarelf = zeros(length(angle),1);
    
    for iang = 1: length(angle)
        %%% propagation trial of hf
        mighfdum = mighf;
        for j = 1: size(mighf,1)
            x0 = mighfdum(j,1);
            y0 = mighfdum(j,2);
            [newx,newy] = coordinate_rot(x0,y0,-(angle(iang)-90),0,0);
            mighfdum(j,1) = newx;
            mighfdum(j,2) = newy;
        end        
        % linear robust least square
        [fitobj,gof,~] = fit(mighfdum(:,15)/3600, mighfdum(:,1),fttpfree,'Robust','Bisquare',...
                                  'StartPoint',[1 1]);
        % output fit parameters
        coef = coeffvalues(fitobj);
        slopehf(iang) = coef(1);
        fitprophf = feval(fitobj,mighfdum(:,15)/3600);
%         l1normhf(iang) = sum(abs(mighfdum(:,1)-fitprophf))/(length(mighfdum(:,1)));
%         l2normhf(iang) = sum((mighfdum(:,1)-fitprophf).^2)/(length(mighfdum(:,1)));
%         ssehf(iang) = gof.sse;
        rmsehf(iang) = gof.rmse;
        rsquarehf(iang) = gof.rsquare;    % r-square, defined as ratio of the sum of squares of the regression and the total sum of squares
    
        %%% propagation trial of lf
        miglfdum = miglf;
        for j = 1: size(miglf,1)
            x0 = miglfdum(j,1);
            y0 = miglfdum(j,2);
            [newx,newy] = coordinate_rot(x0,y0,-(angle(iang)-90),0,0);
            miglfdum(j,1) = newx;
            miglfdum(j,2) = newy;
        end
        % linear robust least square
        [fitobj,gof,~] = fit(miglfdum(:,15)/3600, miglfdum(:,1),fttpfree,'Robust','Bisquare',...
                                  'StartPoint',[0.5 5]);
        % output fit parameters
        coef = coeffvalues(fitobj);
        slopelf(iang) = coef(1);
        fitproplf = feval(fitobj,miglfdum(:,15)/3600);
%         l1normlf(iang) = sum(abs(miglfdum(:,1)-fitproplf))/(length(miglfdum(:,1)));
%         l2normlf(iang) = sum((miglfdum(:,1)-fitproplf).^2)/(length(miglfdum(:,1)));
%         sself(iang) = gof.sse;
        rmself(iang) = gof.rmse;
        rsquarelf(iang) = gof.rsquare;    % r-square, defined as ratio of the sum of squares of the regression and the total sum of squares
        
    end

    %%% best angle estimate from hf
    ind = find(slopehf>0);
%     ind1 = find(l1normhf(ind)==min(l1normhf(ind)));
%     angbestl1(indplt(i)) = angle(ind(ind1(1)));
%     ind2 = find(l2normhf(ind)==min(l2normhf(ind)));
%     angbestl2(indplt(i)) = angle(ind(ind2(1)));
    ind3 = find(rmsehf(ind)==min(rmsehf(ind)));     % rmse, Root Mean Squared Error, the smaller, the better
    if length(ind3) > 1
        disp(strcat({'multiple angles have same rmse: '},num2str(angle(ind3))));
    end
    angbestl3(indplt(i)) = angle(ind(ind3(1)));
%     ind4 = find(rsquarehf(ind)==max(rsquarehf(ind)));   % R-square, 0-1, the larger the better
%     angbestl4(indplt(i)) = angle(ind(ind4(1)));
%     ind5 = find(ssehf(ind)==min(ssehf(ind)));   % sum of square due to error, the smaller, the better
%     angbestl5(indplt(i)) = angle(ind(ind5(1)));
%     disp([angbestl3(indplt(i)),angbestl4(indplt(i)),angbestl5(indplt(i))])
%     angbesthf(indplt(i)) = median([angbestl3(indplt(i)),angbestl4(indplt(i)),angbestl5(indplt(i))]);

    ind6 = find(slopehf==max(slopehf)); % one with the largest slope, i.e., migrating speed
    if length(ind6) > 1
        disp(strcat({'multiple angles have same slope:'},num2str(angle(ind6))));
    end
    angbestl6(indplt(i)) = angle(ind6(1));
    
    angbesthf(indplt(i)) = angbestl3(indplt(i));
    
    
    %%% best angle estimate from hf
    ind = find(slopelf>0);
%     ind1 = find(l1normlf(ind)==min(l1normlf(ind)));
%     angbestl1(indplt(i)) = angle(ind(ind1(1)));
%     ind2 = find(l2normlf(ind)==min(l2normlf(ind)));
%     angbestl2(indplt(i)) = angle(ind(ind2(1)));
    ind3 = find(rmself(ind)==min(rmself(ind)));     % rmse, Root Mean Squared Error, the smaller, the better
    angbestl3(indplt(i)) = angle(ind(ind3(1)));
%     ind4 = find(rsquarelf(ind)==max(rsquarelf(ind)));   % R-square, 0-1, the larger the better
%     angbestl4(indplt(i)) = angle(ind(ind4(1)));
%     ind5 = find(sself(ind)==min(sself(ind)));   % sum of square due to error, the smaller, the better
%     angbestl5(indplt(i)) = angle(ind(ind5(1)));
%     disp([angbestl3(indplt(i)),angbestl4(indplt(i)),angbestl5(indplt(i))])
%     angbestlf(indplt(i)) = median([angbestl3(indplt(i)),angbestl4(indplt(i)),angbestl5(indplt(i))]);
    angbestlf(indplt(i)) = angbestl3(indplt(i));
    
    % determine if more inspection is needed
    if abs(angbesthf(indplt(i))-angbestlf(indplt(i))) > 20
        disp('Difference in propagation direction between hf & lf is too large!');
    end
    
end

angbesthf = angbesthf';
angbestlf = angbestlf';

% angbest = [angbesthf angbestlf];
angbest = angbesthf;

%
for i = 1: size(trange,1)
% for i = 1: 211
%     i=5;
    disp(trange(i,:));
    indhf = find(hftime(:,13)==trange(i,1) & hftime(:,15)>=trange(i,2) & ...
                 hftime(:,15)<=trange(i,3));
    mighf = hftime(indhf,:);
    
    mighfdum = mighf;
    for j = 1: size(mighf,1)
        x0 = mighfdum(j,1);
        y0 = mighfdum(j,2);
        [newx,newy] = coordinate_rot(x0,y0,-(angbest(i)-90),0,0);
        mighfdum(j,1) = newx;
        mighfdum(j,2) = newy;
    end
    
    indlf = find(lftime(:,13)==trange(i,1) & lftime(:,15)>=trange(i,2) & ...
                 lftime(:,15)<=trange(i,3));
    miglf = lftime(indlf,:); 
    
    miglfdum = miglf;
    for j = 1: size(miglf,1)
        x0 = miglfdum(j,1);
        y0 = miglfdum(j,2);
        [newx,newy] = coordinate_rot(x0,y0,-(angbest(i)-90),0,0);
        miglfdum(j,1) = newx;
        miglfdum(j,2) = newy;
    end
    
    %%% define and position the figure frame and axes of each plot
    f.fig=figure;
    widin = 8;  % maximum width allowed is 8.5 inches
    htin = 9;   % maximum height allowed is 11 inches
    set(f.fig,'Position',[scrsz(1)+1*scrsz(3)/10 scrsz(2)+scrsz(4)/20 widin*res htin*res]);
    nrow = 2;
    ncol = 2;    
    for isub = 1:nrow*ncol
        f.ax(isub) = subplot(nrow,ncol,isub);
    end

    %%% reposition
    set(f.ax(1), 'position', [ 0.1, 0.64, 0.32, 0.32]);
    set(f.ax(2), 'position', [ 0.5, 0.64, 0.32, 0.32]);
    set(f.ax(3), 'position', [ 0.1, 0.35, 0.32, 0.22]);
    set(f.ax(4), 'position', [ 0.5, 0.35, 0.32, 0.22]);
%     set(f.ax(5), 'position', [ 0.1, 0.078, 0.32, 0.22]);
%     set(f.ax(6), 'position', [ 0.5, 0.078, 0.32, 0.22]);
    
    % subplot 1 of figure i
    hold(f.ax(1),'on');
    plot(f.ax(1),[-100 100],[0 0],'k--');
    plot(f.ax(1),[0 0],[-100 100],'k--');
    f.ax(1).FontSize = 9;
    mighf = sortrows(mighf,-15);
    scatter(f.ax(1),mighf(:,1),mighf(:,2), 20, mighf(:,15)/3600, 'filled','o');  %, 'MarkerEdgeColor', 'w')
    colormap(f.ax(1),'jet');
    c=colorbar(f.ax(1),'SouthOutside');
    pos = f.ax(1).Position;
    c.Position = [pos(1), pos(2)-0.018, pos(3), 0.02];
%     c.TickLabels=[];
    juldate = num2str(trange(i,1));
    yr = str2double(juldate(1:4));
    date = str2double(juldate(5:end));
    a = jul2dat(yr,date);
    mo = a(1);
    if mo == 9 
        mo = {' Sep. '};
    elseif mo == 7
        mo = {' Jul. '};
    else
        mo = {' Mar. '};
    end
    day = num2str(a(2));
    yr = num2str(a(3));
    c.Label.String = strcat({'Time (hr) on '}, day, mo, yr);
%     c.Label.String = strcat(num2str(trange(i,1)),' of HF',' (hr)');
    c.Label.FontSize = 11;
    caxis(f.ax(1),[trange(i,2)/3600 trange(i,3)/3600])
    text(f.ax(1),0.85,0.1,'HF','FontSize',12,'unit','normalized');
    text(f.ax(1),0.04,0.93,'a','FontSize',11,'unit','normalized','EdgeColor','k','Margin',2);
    text(f.ax(1),0.04,0.1,'LZB','FontSize',11,'unit','normalized','EdgeColor','k','Margin',2);
    f.ax(1).Box = 'on';
    grid(f.ax(1), 'on');
    axis(f.ax(1), 'equal');
    f.ax(1).GridLineStyle = '--';
    f.ax(1).XAxisLocation = 'top';
    medxhf = median(mighf(:,1));
    medyhf = median(mighf(:,2));
    [rotx, roty] = complex_rot(0,5,-angbest(i));
    xvect = [medxhf-rotx medxhf+rotx];
    yvect = [medyhf-roty medyhf+roty];
    drawArrow(f.ax(1),xvect,yvect,xran,yran,'linewidth',1.5);
    [rotx, roty] = complex_rot(-5,0,-angbest(i));
    xvect = [medxhf-rotx medxhf+rotx];
    yvect = [medyhf-roty medyhf+roty];
    drawArrow(f.ax(1),xvect,yvect,xran,yran,'linewidth',1.5,'linestyle','--','color',[0.4 0.4 0.4]);
    xticks(f.ax(1),xran(1):5:xran(2));
    yticks(f.ax(1),yran(1):5:yran(2));    
    xlabel(f.ax(1),'E (km)','fontsize',11);
    ylabel(f.ax(1),'N (km)','fontsize',11);
    text(f.ax(1),0.5,0.1,num2str(i),'FontSize',12,'unit','normalized');
    text(f.ax(1),0.83,0.95,strcat(num2str(size(mighf,1)),{' detections'}),'FontSize',8,...
         'unit','normalized','horizontalalignment','center');
    text(f.ax(1),0.83,0.90,strcat({'in '},num2str(trange(i,3)-trange(i,2)),{' s'}),'FontSize',8,...
         'unit','normalized','horizontalalignment','center');
    rate = sprintf('%.3f',size(mighf,1)/(trange(i,3)-trange(i,2)));
    text(f.ax(1),0.83,0.85,strcat({'rate: '},rate),'FontSize',8,...
         'unit','normalized','horizontalalignment','center');
    text(f.ax(1),0.83,0.75,strcat(num2str(angbest(i)),{' ^o'}),'FontSize',10,...
         'unit','normalized','horizontalalignment','center'); 
    hold(f.ax(1),'off');

    % subplot 2 of figure i
    hold(f.ax(2),'on');
    plot(f.ax(2),[-100 100],[0 0],'k--');
    plot(f.ax(2),[0 0],[-100 100],'k--');
    f.ax(2).FontSize = 9;
    miglf = sortrows(miglf,-15);
    scatter(f.ax(2),miglf(:,1),miglf(:,2), 20, miglf(:,15)/3600, 'filled','o');  %, 'MarkerEdgeColor', 'w')
    colormap(f.ax(2),'jet');
    c=colorbar(f.ax(2),'SouthOutside');
    pos = f.ax(2).Position;
    c.Position = [pos(1), pos(2)-0.018, pos(3), 0.02];
%     c.TickLabels=[];
    c.Label.String = strcat({'Time (hr) on '}, day, mo, yr);
%     c.Label.String = strcat(num2str(trange(i,1)),' of LF',' (hr)');
    c.Label.FontSize = 11;
    caxis(f.ax(2),[trange(i,2)/3600 trange(i,3)/3600])
    text(f.ax(2),0.85,0.1,'LF','FontSize',12,'unit','normalized');
    text(f.ax(2),0.04,0.93,'b','FontSize',11,'unit','normalized','EdgeColor','k','Margin',2);
    text(f.ax(2),0.04,0.1,'LZB','FontSize',11,'unit','normalized','EdgeColor','k','Margin',2);
    f.ax(2).Box = 'on';
    grid(f.ax(2), 'on');
    axis(f.ax(2), 'equal');
    f.ax(2).GridLineStyle = '--';
    f.ax(2).XAxisLocation = 'top';
    medxlf = median(miglf(:,1));
    medylf = median(miglf(:,2));
    [rotx, roty] = complex_rot(0,5,-angbest(i));
    xvect = [medxlf-rotx medxlf+rotx];
    yvect = [medylf-roty medylf+roty];
    drawArrow(f.ax(2),xvect,yvect,xran,yran,'linewidth',1.5);
    [rotx, roty] = complex_rot(-5,0,-angbest(i));
    xvect = [medxlf-rotx medxlf+rotx];
    yvect = [medylf-roty medylf+roty];
    drawArrow(f.ax(2),xvect,yvect,xran,yran,'linewidth',1.5,'linestyle','--','color',[0.4 0.4 0.4]);
    xticks(f.ax(2),xran(1):5:xran(2));
    yticks(f.ax(2),yran(1):5:yran(2));    
    xlabel(f.ax(2),'E (km)','fontsize',11);
    ylabel(f.ax(2),'N (km)','fontsize',11);
    text(f.ax(2),0.5,0.1,num2str(i),'FontSize',12,'unit','normalized');
    text(f.ax(2),0.83,0.95,strcat(num2str(size(miglf,1)),{' detections'}),'FontSize',8,...
         'unit','normalized','horizontalalignment','center');
    text(f.ax(2),0.83,0.9,strcat({'in '},num2str(trange(i,3)-trange(i,2)),{' s'}),'FontSize',8,...
         'unit','normalized','horizontalalignment','center');
    rate = sprintf('%.3f',size(miglf,1)/(trange(i,3)-trange(i,2)));
    text(f.ax(2),0.83,0.85,strcat({'rate: '},rate),'FontSize',8,...
         'unit','normalized','horizontalalignment','center');
    text(f.ax(2),0.83,0.75,strcat(num2str(angbest(i)),{' ^o'}),'FontSize',10,...
         'unit','normalized','horizontalalignment','center');
    hold(f.ax(2),'off');
    
    % subplot 3 of figure i
    hold(f.ax(3),'on');
    f.ax(3).FontSize = 9;
%     scatter(f.ax(3),mighfdum(:,15)/3600,mighfdum(:,1),20,[1 0 0],'filled','o',...
%             'MarkerEdgeColor','k'); 
    scatter(f.ax(3),mighfdum(:,15)/3600,mighfdum(:,1),20,[0.6 0.6 0.6],'filled','o',...
            'MarkerEdgeColor','k');    
    scatter(f.ax(3),miglfdum(:,15)/3600,miglfdum(:,1),20,[0.6 1 1],'filled','o',...
            'MarkerEdgeColor','k');
    
%     indpool = [12;13;16;19;21;22;29;32;36;40];
%     indpool = [indpool; 3;5;7;17;20;23;25;26;30];
%     if ismember(i,indpool)
        mintime = max(min(mighfdum(:,15)), min(miglfdum(:,15)));
        maxtime = min(max(mighfdum(:,15)), max(miglfdum(:,15)));
        mighfdum2 = mighfdum(mighfdum(:,15)>=mintime & mighfdum(:,15)<=maxtime, :);
        miglfdum2 = miglfdum(miglfdum(:,15)>=mintime & miglfdum(:,15)<=maxtime, :);
%     else
%         mighfdum2 = mighfdum;
%         miglfdum2 = miglfdum;
%     end

    [fitobjhfprop,gofhfprop,outphfprop] = fit(mighfdum2(:,15)/3600, mighfdum2(:,1),fttpfree,'Robust',...
                                              'Bisquare','StartPoint',[1 1]);
    % output fit parameters
    coefprop = coeffvalues(fitobjhfprop);
    slopeprophf(i) = coefprop(1);
    intcptprophf(i) = coefprop(2);
    fitprophf = feval(fitobjhfprop,mighfdum2(:,15)/3600);
    plot(f.ax(3),mighfdum2(:,15)/3600,fitprophf,'-','linewidth',2,'color',[0.6 0.6 0.6]);
        
%     intcptproplf = ones(size(miglfdum(:,end)))\(miglfdum(:,1)-slopeprophf*miglfdum(:,end)/3600); % direct LS fit
%     fitproplf = slopeprophf*miglfdum(:,end)/3600+intcptproplf;
    a=slopeprophf(i);
    fttpfix = fittype( @(b,x) a*x+b);
    [fitobjlfprop,goflfprop,outplfprop] = fit(miglfdum2(:,15)/3600, miglfdum2(:,1),fttpfix,'Robust',...
                                              'Bisquare','StartPoint',1);
%     [fitobjlfprop,goflfprop,outplfprop] = fit(miglfdum(:,10)/3600, miglfdum(:,1),fttpfix,'Robust',...
%                                               'Bisquare','Weights',miglfdum(:,11));
    fitproplf = feval(fitobjlfprop,miglfdum2(:,15)/3600);
    intcptproplf(i) = coeffvalues(fitobjlfprop);
    offset(i) = intcptproplf(i)-intcptprophf(i);
    plot(f.ax(3),miglfdum2(:,15)/3600,fitproplf,'-','linewidth',2,'color',[0.6 1 1]);
    f.ax(3).Box = 'on';
    grid(f.ax(3), 'on');
    f.ax(3).GridLineStyle = '--';
    xran1 = [trange(i,2)/3600 trange(i,3)/3600];
    yran1 = [round(min([mighfdum(:,1);miglfdum(:,1)]))-1 ...
             round(max([mighfdum(:,1);miglfdum(:,1)]))+1];
    xlim(f.ax(3),xran1);
    ylim(f.ax(3),yran1);
    text(f.ax(3),0.04,0.91,'c','FontSize',11,'unit','normalized','EdgeColor','k','Margin',2);
    xlabel(f.ax(3),'Time (hr)','fontsize',11);
    ylabel(f.ax(3),'Dist. along prop. (km)','fontsize',11); 

    
    % compute the HF weights in robust linear regression, see NOTES above
    rhf = outphfprop.residuals;   % usual residuals
    x = [mighfdum2(:,15)/3600];
    hatmat = x*inv(x'*x)*x';
    h = zeros(size(hatmat,1),1);    % leverage of least square
    for jj = 1 : size(hatmat,1)
        h(jj) = hatmat(jj,jj);
    end    
    radj = rhf./sqrt(1-h);      % adjusted residuals
    K = 4.685;
    s = mad(rhf,1)/0.6745;
    u = radj/(K*s);
    wthf = zeros(length(u),1);    % rubust weight of next iteration
    for jj = 1 : length(u)
        if abs(u(jj)) < 1
            wthf(jj) = (1-(u(jj))^2)^2;
        else
            wthf(jj) = 0;
        end
    end
    
    % get the standard error of the estimated parameters, may indicate the compare the quality
    % of fitting cross different RTMs, NOTE that rmse, mse are definitely not a good measure
    % the obtained CI is comfirmed to be correct by comparing it with the 'fitobjhfprop'
    slopese = gofhfprop.rmse./sqrt(sum((x-mean(x)).^2));
    est = slopeprophf(i);
    slopeprophfCI(i,:) = confidence_interval_general(est,slopese,length(x)-2,95);
    interceptse = slopese.*sqrt(sum(x.^2)./length(x));
    est = intcptprophf(i);
    intcptprophfCI(i,:) = confidence_interval_general(est,interceptse,length(x)-2,95);
    sehf(i, :) = [slopese interceptse]; 
    
    x = mighfdum2(:,15)/3600;
    y = mighfdum2(:,1);
    x_bar = wt_mean(x,wthf);
    y_bar = wt_mean(y,wthf);
    x_var = sum(wthf.*(x-x_bar).^2) / sum(wthf);
    y_var = sum(wthf.*(y-y_bar).^2) / sum(wthf);
    xy_cov = sum(wthf.*(x-x_bar).*(y-y_bar)) / sum(wthf);
    pearwthf(i) = xy_cov / sqrt(x_var*y_var);
    
    % compute the LF weights in robust linear regression, see NOTES above
    rlf = outplfprop.residuals;   % usual residuals
    x = [miglfdum2(:,15)/3600];
    hatmat = x*inv(x'*x)*x';
    h = zeros(size(hatmat,1),1);    % leverage of least square
    for jj = 1 : size(hatmat,1)
        h(jj) = hatmat(jj,jj);
    end    
    radj = rlf./sqrt(1-h);      % adjusted residuals
    K = 4.685;
    s = mad(rlf,1)/0.6745;
    u = radj/(K*s);
    wtlf = zeros(length(u),1);    % rubust weight of next iteration
    for jj = 1 : length(u)
        if abs(u(jj)) < 1
            wtlf(jj) = (1-(u(jj))^2)^2;
        else
            wtlf(jj) = 0;
        end
    end
    
    text(f.ax(3),0.15,0.91,strcat({'Slope: '}, sprintf('%.1f',a)),'FontSize',8,...
         'unit','normalized','horizontalalignment','left');
    text(f.ax(3),0.9,0.91,strcat({'SE: '}, sprintf('%.1f',slopese)),'FontSize',8,...
         'unit','normalized','horizontalalignment','right'); 
    text(f.ax(3),0.15,0.85,strcat({'Pearson: '}, sprintf('%.2f',pearwthf(i))),'FontSize',8,...
         'unit','normalized','horizontalalignment','left');
    hold(f.ax(3),'off');
    
    
    % subplot 4 of figure i
    hold(f.ax(4),'on');
    f.ax(4).FontSize = 9;
    numhf(i) = length(mighfdum2(:,1));
    numlf(i) = length(miglfdum2(:,1));
    resprophf(i,1:numhf(i)) = outphfprop.residuals;
    resproplf(i,1:numlf(i)) = outplfprop.residuals+(intcptproplf(i)-intcptprophf(i));
    [xlochf, a, b, pdfvalhf] = weighted_bar_pdf(resprophf(i,1:numhf(i)), wthf, 0.5, 'int');
    hfHdl = bar(f.ax(4),xlochf, pdfvalhf,1,'stacked');
    hfHdl(1).FaceColor = [0.6 0.6 0.6];
    [xloclf, ~, ~, pdfvallf] = weighted_bar_pdf(resproplf(i,1:numlf(i)), wtlf, 0.5, 'int');
    lfHdl = bar(f.ax(4),xloclf, pdfvallf,1,'stacked','facea',0.6);
    lfHdl(1).FaceColor = [0.6 1 1];
    text(f.ax(4),0.04,0.91,'d','FontSize',11,'unit','normalized','EdgeColor','k','Margin',2);
    text(f.ax(4),0.4,0.92,sprintf('LF-HF = %.2f km',offset(i)),'FontSize',11,'unit','normalized');
    f.ax(4).Box = 'on';
    grid(f.ax(4), 'on');
    f.ax(4).GridLineStyle = '--';
    ymax = f.ax(4).YLim(2)+0.1;
    ylim(f.ax(4),[0 ymax]);
%     ylim(f.ax(4),[0 1]);
    xlabel(f.ax(4),'Residual in prop. (km)','fontsize',11);
    ylabel(f.ax(4),'PDF estimate','fontsize',11);
    hold(f.ax(4),'off');

%     if pearwthf(i) >= 0.5 && sehf(i,1)<3
%         print(f.fig,'-dpdf',strcat('/home/data2/chaosong/CurrentResearch/Song_Rubin_2020/Figures/nrtm03_',...
%             num2str(i),'.pdf'));
%     end

end

[tmp,tmpind] = sort(sehf(:,1),'descend');
sortsehf = [tmpind tmp];
[tmp,tmpind] = sort(pearwthf(:),'ascend');
sortpearwthf = [tmpind tmp];


tmp1 = sortpearwthf(sortpearwthf(:,2)>=0.5, 1);
tmp2 = sortsehf(sortsehf(:,2)<3, 1);

if issame(trange,ntran04)
    ind04 = intersect(tmp1,tmp2);
    acctran04 = ntran04(ind04,:);
elseif issame(trange,ntran05)
    ind05 = intersect(tmp1,tmp2);
    acctran05 = ntran05(ind05,:);
else
    ind03 = intersect(tmp1,tmp2);
    acctran03 = ntran03(ind03,:);
end




















































