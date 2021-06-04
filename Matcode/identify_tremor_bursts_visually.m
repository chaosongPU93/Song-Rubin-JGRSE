% function identify_tremor_bursts_visually
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

set(0,'DefaultFigureVisible','on');
% set(0,'DefaultFigureVisible','off');   % switch to show the plots or not

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

%% previously visually chosen RTM periods
trange = [
           2004195,4.15e+4,4.60e+4;
           2004196,7.55e+4,7.70e+4;
           2004197,0.03e+4,0.35e+4;
           2004197,0.50e+4,0.75e+4;
           2004197,4.60e+4,4.72e+4;
           2004197,8.388e+4,8.46e+4;
           2004197,8.488e+4,8.52e+4;
           2004197,8.55e+4,8.64e+4;
           2004198,0.58e+4,0.98e+4;
           2004198,1.90e+4,2.18e+4;
           2004198,5.40e+4,5.80e+4;
           2004198,7.35e+4,7.60e+4;
           2004198,8.35e+4,8.64e+4;
           2004199,0.50e+4,0.63e+4;
           2004199,4.61e+4,4.91e+4;
           2004199,4.98e+4,5.08e+4;
           2004199,8.10e+4,8.20e+4;
           2004200,1.38e+4,1.80e+4;
           2004200,1.85e+4,1.99e+4;
           2004203,1.60e+4,2.20e+4;
           2005255,3.42e+4,3.58e+4;
           2005255,6.70e+4,6.85e+4;
           2005255,7.50e+4,7.57e+4;
           2005255,8.42e+4,8.56e+4;
           2005256,0.35e+4,0.50e+4;
           2005256,0.80e+4,0.98e+4;
           2005256,2.15e+4,2.23e+4;
           2005256,2.82e+4,2.95e+4;
           2005256,34776,3.53e+4;
           2005256,5.20e+4,5.26e+4;
           2005256,7.275e+4,7.40e+4;
           2005256,7.60e+4,7.80e+4;
           2005257,0.60e+4,0.70e+4;
           2005257,2.10e+4,2.50e+4;
           2005257,3.90e+4,4.40e+4;
           2005257,6.10e+4,6.40e+4;
           2005257,7.36e+4,7.80e+4;
           2005258,3.52e+4,3.66e+4;
           2005258,3.80e+4,4.00e+4;
           2005259,2340,3312;
           2005259,0.50e+4,0.77e+4;
           2005259,3.70e+4,4.26e+4;
           2005259,7.20e+4,7.58e+4;
           2005260,0.36e+4,0.43e+4;
           2005260,0.43e+4,0.57e+4;
           2005260,5.63e+4,5.82e+4;
           2005261,0.82e+4,1.10e+4;];

it03 = find(trange(:,1) < 2004*1000);
rtmt(it03,1:2) = (trange(it03,1)-2003060)+trange(it03,2:3)./(3600.*24);
% rtmt(it03,2) = (trange(it03,1)-2003060)+trange(it03,3)./(3600.*24);
it04 = find(trange(:,1) < 2005*1000 & trange(:,1) > 2004*1000);
rtmt(it04,1:2) = (trange(it04,1)-2004194)+trange(it04,2:3)./(3600.*24);
% rtmt(it04,2) = (trange(it04,1)-2004194)+trange(it04,3)./(3600.*24);
it05 = find(trange(:,1) > 2005*1000);
rtmt(it05,1:2) = (trange(it05,1)-2005254)+trange(it05,2:3)./(3600.*24);


%% possbile new time periods for visual check
newt04 = [
          1.33 1.43;
          2.35 2.45;
          2.5 2.6;
          2.82 2.88;
          3.09 3.1;
          3.13 3.17;
          3.21 3.23;
          3.28 3.31;
          3.34 3.38;
          3.385 3.41;
          3.41 3.44;
          3.47 3.51;
          3.515 3.53;
          3.58 3.61;
          3.63 3.66;
          3.81 3.85;
          3.9 3.92;
          3.92 3.95;
          4.11 4.13;
          4.25 4.28;
          4.37 4.41;
          4.485 4.51;
          4.7 4.75;
          5.03 5.05;
          5.075 5.09;
          5.2 5.22;
          5.37 5.4;
          5.41 5.43;
          5.49 5.51;
          5.66 5.71;
          5.77 5.79;
          5.95 5.97;
          6.1 6.14;
          6.31 6.34;
          6.48 6.53;
          6.55 6.61;
          7 7.05;
          7.05 7.11;
          7.28 7.31;
          8.17 8.2;
          8.53 8.56;
          8.585 8.61;
          9.68 9.73;
          9.75 9.79;
         ];

ntran04(:,1) = floor(newt04(:,1))+2004194;     
ntran04(:,2:3) = (newt04(:,1:2)-floor(newt04(:,1))).*(3600.*24);


newt05 = [
          0.23 0.25;
          0.6 0.63;
          0.78 0.82;
          1.05 1.07;
          1.2 1.22;
          1.255 1.27;
          1.5 1.53;
          1.59 1.615;
          1.675 1.69;
          1.69 1.74;
          1.83 1.85;
          2.11 2.19;
          2.21 2.23;
          2.33 2.37;
          2.37 2.39;
          2.41 2.43;
          2.52 2.54;
          2.805 2.815;
          2.9 2.925;
          2.935 2.98;
          3 3.2;
          3.295 3.305;
          4.285 4.31;
          6 6.03;
          6.1 6.11;
          6.46 6.5;
          7.115 7.15;
         ];

ntran05(:,1) = floor(newt05(:,1))+2005254;     
ntran05(:,2:3) = (newt05(:,1:2)-floor(newt05(:,1))).*(3600.*24);

modt04 = [
          1.43 1.55;
         ];
mtran04(:,1) = floor(modt04(:,1))+2004194;     
mtran04(:,2:3) = (modt04(:,1:2)-floor(modt04(:,1))).*(3600.*24);

modt05 = [
          2.82 2.856;
          3.451 3.72;
          7.095 7.115;
         ];
mtran05(:,1) = floor(modt05(:,1))+2005254;     
mtran05(:,2:3) = (modt05(:,1:2)-floor(modt05(:,1))).*(3600.*24);


%% Along strike, 2003, hf & lf

msizehf = 2;
msizelf = 4;
dt = 2;
mint = 0;
maxt = ceil(max(hftime(ih03,26)));
yran = [-20 20];
[f] = plt_time_dist(hftime(ih03,26),hftime(ih03,2),lftime(il03,26),lftime(il03,2),rtmt(it03,:),...
                    [],msizehf,msizelf,dt,mint,maxt,yran);

text(f.ax(1),0.9,0.9,'2003','FontSize',10,'unit','normalized','horizontalalignment','left',...
        'EdgeColor','k','Margin',2);
ylabel(f.ax(end), 'Along strike, N (km)');
xlabel(f.ax(end), 'Time (day) since 2003060 Mar. 1, 2003');

% print(f.fig,'-depsc2',strcat(rstpath,'/dcutstrikehflf2003.eps'));
print(f.fig,'-depsc2',strcat(rstpath,'/dcutstrikehflf2003_new.eps'));


%% Along dip, 2003, hf & lf

msizehf = 2;
msizelf = 4;
dt = 2;
mint = 0;
maxt = ceil(max(hftime(ih03,26)));
yran = [-20 20];
[f] = plt_time_dist(hftime(ih03,26),hftime(ih03,1),lftime(il03,26),lftime(il03,1),rtmt(it03,:),...
                    [],msizehf,msizelf,dt,mint,maxt,yran);

text(f.ax(1),0.9,0.9,'2003','FontSize',10,'unit','normalized','horizontalalignment','left',...
        'EdgeColor','k','Margin',2);
ylabel(f.ax(end), 'Along dip, E (km)');
xlabel(f.ax(end), 'Time (day) since 2003060 Mar. 1, 2003');

% print(f.fig,'-depsc2',strcat(rstpath,'/dcutdiphflf2003.eps'));
print(f.fig,'-depsc2',strcat(rstpath,'/dcutdiphflf2003_new.eps'));


%% Along strike, 2004, hf & lf

msizehf = 2;
msizelf = 4;
dt = 2;
mint = 0;
maxt = 6;
yran = [-30 30];
[f] = plt_time_dist(hftime(ih04,26),hftime(ih04,2),lftime(il04,26),lftime(il04,2),rtmt(it04,:),...
                    newt04,msizehf,msizelf,dt,mint,maxt,yran);

text(f.ax(1),0.9,0.9,'2004','FontSize',10,'unit','normalized','horizontalalignment','left',...
        'EdgeColor','k','Margin',2);
ylabel(f.ax(end), 'Along strike, N (km)');
xlabel(f.ax(end), 'Time (day) since 2004194 Jul. 12, 2004');

% print(f.fig,'-depsc2',strcat(rstpath,'/dcutstrikehflf2004p1.eps'));
print(f.fig,'-depsc2',strcat(rstpath,'/dcutstrikehflf2004p1_new.eps'));

mint = 6;
maxt = ceil(max(hftime(ih04,26)));
[f] = plt_time_dist(hftime(ih04,26),hftime(ih04,2),lftime(il04,26),lftime(il04,2),rtmt(it04,:),...
                    newt04,msizehf,msizelf,dt,mint,maxt,yran);

text(f.ax(1),0.9,0.9,'2004','FontSize',10,'unit','normalized','horizontalalignment','left',...
        'EdgeColor','k','Margin',2);
ylabel(f.ax(end), 'Along strike, N (km)');
xlabel(f.ax(end), 'Time (day) since 2004194 Jul. 12, 2004');

% print(f.fig,'-depsc2',strcat(rstpath,'/dcutstrikehflf2004p2.eps'));
print(f.fig,'-depsc2',strcat(rstpath,'/dcutstrikehflf2004p2_new.eps'));


%% Along dip, 2004, hf & lf

msizehf = 2;
msizelf = 4;
dt = 2;
mint = 0;
maxt = 6;
yran = [-30 30];
[f] = plt_time_dist(hftime(ih04,26),hftime(ih04,1),lftime(il04,26),lftime(il04,1),rtmt(it04,:),...
                    newt04,msizehf,msizelf,dt,mint,maxt,yran);

text(f.ax(1),0.9,0.9,'2004','FontSize',10,'unit','normalized','horizontalalignment','left',...
        'EdgeColor','k','Margin',2);
ylabel(f.ax(end), 'Along dip, E (km)');
xlabel(f.ax(end), 'Time (day) since 2004194 Jul. 12, 2004');

% print(f.fig,'-depsc2',strcat(rstpath,'/dcutdiphflf2004p1.eps'));
print(f.fig,'-depsc2',strcat(rstpath,'/dcutdiphflf2004p1_new.eps'));

mint = 6;
maxt = ceil(max(hftime(ih04,26)));
[f] = plt_time_dist(hftime(ih04,26),hftime(ih04,1),lftime(il04,26),lftime(il04,1),rtmt(it04,:),...
                    newt04,msizehf,msizelf,dt,mint,maxt,yran);

text(f.ax(1),0.9,0.9,'2004','FontSize',10,'unit','normalized','horizontalalignment','left',...
        'EdgeColor','k','Margin',2);
ylabel(f.ax(end), 'Along dip, E (km)');
xlabel(f.ax(end), 'Time (day) since 2004194 Jul. 12, 2004');

% print(f.fig,'-depsc2',strcat(rstpath,'/dcutdiphflf2004p2.eps'));
print(f.fig,'-depsc2',strcat(rstpath,'/dcutdiphflf2004p2_new.eps'));


%% Along strike, 2005, hf & lf

msizehf = 2;
msizelf = 4;
dt = 2;
mint = 0;
maxt = 6;
yran = [-30 30];
[f] = plt_time_dist(hftime(ih05,26),hftime(ih05,2),lftime(il05,26),lftime(il05,2),rtmt(it05,:),...
                    newt05,msizehf,msizelf,dt,mint,maxt,yran);

text(f.ax(1),0.9,0.9,'2005','FontSize',10,'unit','normalized','horizontalalignment','left',...
        'EdgeColor','k','Margin',2);
ylabel(f.ax(end), 'Along strike, N30W (km)');
xlabel(f.ax(end), 'Time (day) since 2005254 Sep. 11, 2005');

% print(f.fig,'-depsc2',strcat(rstpath,'/dcutstrikehflf2005p1.eps'));
print(f.fig,'-depsc2',strcat(rstpath,'/dcutstrikehflf2005p1_new.eps'));

mint = 6;
maxt = ceil(max(hftime(ih05,26)));
[f] = plt_time_dist(hftime(ih05,26),hftime(ih05,2),lftime(il05,26),lftime(il05,2),rtmt(it05,:),...
                    newt05,msizehf,msizelf,dt,mint,maxt,yran);

text(f.ax(1),0.9,0.9,'2005','FontSize',10,'unit','normalized','horizontalalignment','left',...
        'EdgeColor','k','Margin',2);
ylabel(f.ax(end), 'Along strike, N30W (km)');
xlabel(f.ax(end), 'Time (day) since 2005254 Sep. 11, 2005');

% print(f.fig,'-depsc2',strcat(rstpath,'/dcutstrikehflf2005p2.eps'));
print(f.fig,'-depsc2',strcat(rstpath,'/dcutstrikehflf2005p2_new.eps'));


%% Along dip, 2005, hf & lf

msizehf = 2;
msizelf = 4;
dt = 2;
mint = 0;
maxt = 6;
yran = [-30 30];
[f] = plt_time_dist(hftime(ih05,26),hftime(ih05,1),lftime(il05,26),lftime(il05,1),rtmt(it05,:),...
                    newt05,msizehf,msizelf,dt,mint,maxt,yran);

text(f.ax(1),0.9,0.9,'2005','FontSize',10,'unit','normalized','horizontalalignment','left',...
        'EdgeColor','k','Margin',2);
ylabel(f.ax(end), 'Along dip, N60E (km)');
xlabel(f.ax(end), 'Time (day) since 2005254 Sep. 11, 2005');

% print(f.fig,'-depsc2',strcat(rstpath,'/dcutdiphflf2005p1.eps'));
print(f.fig,'-depsc2',strcat(rstpath,'/dcutdiphflf2005p1_new.eps'));

mint = 6;
maxt = ceil(max(hftime(ih05,26)));
[f] = plt_time_dist(hftime(ih05,26),hftime(ih05,1),lftime(il05,26),lftime(il05,1),rtmt(it05,:),...
                    newt05,msizehf,msizelf,dt,mint,maxt,yran);

text(f.ax(1),0.9,0.9,'2005','FontSize',10,'unit','normalized','horizontalalignment','left',...
        'EdgeColor','k','Margin',2);
ylabel(f.ax(end), 'Along dip, N60E (km)');
xlabel(f.ax(end), 'Time (day) since 2005254 Sep. 11, 2005');

% print(f.fig,'-depsc2',strcat(rstpath,'/dcutdiphflf2005p2.eps'));
print(f.fig,'-depsc2',strcat(rstpath,'/dcutdiphflf2005p2_new.eps'));












