% function mig_LZB_ortho
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script is used to use the poropagation direction estimate fit the hf
% the migration linearly to get the propagation slope and error distribution.
%
% Chao Song, chaosong@princeton.edu
% First created date:   2020/08/13
% Last modified date:   2020/08/13
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Initialization
format short e   % Set the format to 5-digit floating point
clear
close all
clc

set(0,'DefaultFigureVisible','on');
% set(0,'DefaultFigureVisible','off');   % switch to show the plots or not

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
            '017'];

%%% corresponds to PERMROTS
PERMSTA=['PGC'        % permanent station names
         'LZB'];
POLSTA=['SSIB '           % polaris station names
        'SILB '
        'KLNB '
        'MGCB '
        'TWKB '];
    
stas=['TWKB '
      'LZB  '
      'MGCB '];     % determine the trio and order         


nfam = length(nfampool);

if nfam == 11
    %%% this is inverted from (0,0) of all fams, same order, location of control points
    %%% inverted from /home/data2/chaosong/Seisbasics/hypoinverse/lzbrst/chat00p.reffam
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
        -123.838500 48.544833 35.6600];
    
elseif nfam == 12
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
        -123.879667 48.446167 34.2600];
end

%%% convert lfe location to km relative to 043
relacont = loccont;
[dx, dy] = absloc2relaloc(loccont(:,1),loccont(:,2),loccont(2,1),loccont(2,2));
relacont(:,1) = dx;
relacont(:,2) = dy;

% load files
winlenhf = 4;

winlenlf = 16;
lofflf = 4;
ccminlf = 0.35;

SUFFIXhf = strcat('up.hf.time.',num2str(winlenhf),'_',num2str(winlenlf),'.',num2str(lofflf),...
                  '.',num2str(ccminlf));
fname = strcat(rstpath, '/evtloc.allfam.dcutnodou.',SUFFIXhf);
hftime = load(fname);
% 25 cols, format is:
%   E(043) N(043) E(own) N(own) dist(own) lon lat dep off12 off13 off12sec off13sec date
%   main_arrival_time  cent_of_win  avecc ccadd1 ccadd2 ccadd3 ccadd4 ampmax ampmin eng normeng famnum
%


SUFFIXlf = strcat('lf.time.',num2str(winlenhf),'_',num2str(winlenlf),'.',num2str(lofflf),...
                  '.',num2str(ccminlf));
fname = strcat(rstpath, '/evtloc.allfam.dcutnodou.',SUFFIXlf);
lftime = load(fname);


% sort according to day, sec
hftime = sortrows(hftime, [13, 15]);
lftime = sortrows(lftime, [13, 15]);


trange = [
           2004195,4.15e+4,4.60e+4;
           2004196,7.55e+4,7.70e+4;
           2004197,0.03e+4,0.35e+4;
           2004197,0.50e+4,0.75e+4;
           2004197,4.60e+4,4.72e+4;
           2004197,8.35e+4,8.46e+4;
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
           2005256,3.45e+4,3.53e+4;
           2005256,5.20e+4,5.26e+4;
           2005256,7.24e+4,7.40e+4;
           2005256,7.60e+4,7.80e+4;
           2005257,0.60e+4,0.70e+4;
           2005257,2.10e+4,2.50e+4;
           2005257,3.90e+4,4.40e+4;
           2005257,6.10e+4,6.40e+4;
           2005257,7.36e+4,7.80e+4;
           2005258,3.52e+4,3.66e+4;
           2005258,3.80e+4,4.00e+4;
           2005259,0.20e+4,0.40e+4;
           2005259,0.50e+4,0.77e+4;
           2005259,3.70e+4,4.26e+4;
           2005259,7.20e+4,7.58e+4;
           2005260,0.36e+4,0.43e+4;
           2005260,0.43e+4,0.57e+4;
           2005260,5.63e+4,5.82e+4;
           2005261,0.82e+4,1.10e+4;];
        
xran = [-15 25];
yran = [-20 20];

angbest = [265;
215;
250;
255;
235;
220;
245;
250;
235;
230;
75;
235;
240;
250;
120;
195;
85;
115;
255;
115;
60;
35;  % weird, check timing, lf is also bad, use neither
70;
240;
215;
215;
115;
110;
250;
260;
270;
250;
195;
245;
165;
245;
240;
90;
185;
135;
255;
80;
240;
245;
70;
235;
70;];




%% linear fitting and median distance 
resprophf = nan(length(trange)+1,200);
resproplf = nan(length(trange)+1,50);

% create fit object with free constraints
fttpfree = fittype( @(a,b,x) a*x+b);
        
for i = 1: length(trange)
% for i = 1: 21
%     i=3;
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
    scatter(f.ax(1),relacont(:,1),relacont(:,2),20,'k','o','LineWidth',1.5);
    
%     [x,y] = circle_chao(relacont(3,1),relacont(3,2),5,1);
    [x,y] = circle_chao(-9,2,5,1);
%     plot(f.ax(1),x,y,'linewidth',1.5,'linestyle','-','color',[0.6 0.6 0.6]);
%     
%     [x,y] = circle_chao(-10,0,2,1);
%     plot(f.ax(1),x,y,'linewidth',1.5,'linestyle','-','color',[0.6 0.6 0.6]);
%     
%     [x,y] = circle_chao(-7.5,3,3,1);
%     plot(f.ax(1),x,y,'linewidth',1.5,'linestyle','-','color',[0.6 0.6 0.6]);

    [x,y] = sector_chao(-9,2,5,120,300,1);
    x = [x x(1)];
    y = [y y(1)];
    plot(f.ax(1),x,y,'linewidth',1.5,'linestyle','-','color',[0.6 0.6 0.6]);

    
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
    text(f.ax(1),0.85,0.15,'HF','FontSize',12,'unit','normalized');
    text(f.ax(1),0.04,0.93,'a','FontSize',11,'unit','normalized','EdgeColor','k','Margin',2);
    f.ax(1).Box = 'on';
    grid(f.ax(1), 'on');
    axis(f.ax(1), 'equal');
    f.ax(1).GridLineStyle = '--';
    f.ax(1).XAxisLocation = 'top';
    medx = median(mighf(:,1));
    medy = median(mighf(:,2));
    [rotx, roty] = complex_rot(0,5,-angbest(i));
    xvect = [medx-rotx medx+rotx];
    yvect = [medy-roty medy+roty];
    drawArrow(f.ax(1),xvect,yvect,xran,yran,'linewidth',1.5);
    [rotx, roty] = complex_rot(-5,0,-angbest(i));
    xvect = [medx-rotx medx+rotx];
    yvect = [medy-roty medy+roty];
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
    hold(f.ax(1),'off');

    % subplot 2 of figure i
    hold(f.ax(2),'on');
    plot(f.ax(2),[-100 100],[0 0],'k--');
    plot(f.ax(2),[0 0],[-100 100],'k--');
    f.ax(2).FontSize = 9;
    miglf = sortrows(miglf,-15);
    scatter(f.ax(2),miglf(:,1),miglf(:,2), 20, miglf(:,15)/3600, 'filled','o');  %, 'MarkerEdgeColor', 'w')
    scatter(f.ax(2),relacont(:,1),relacont(:,2),20,'k','o','LineWidth',1.5);
    
    plot(f.ax(2),x,y,'linewidth',1.5,'linestyle','-','color',[0.6 0.6 0.6]);

    colormap(f.ax(2),'jet');
    c=colorbar(f.ax(2),'SouthOutside');
    pos = f.ax(2).Position;
    c.Position = [pos(1), pos(2)-0.018, pos(3), 0.02];
%     c.TickLabels=[];
    c.Label.String = strcat({'Time (hr) on '}, day, mo, yr);
%     c.Label.String = strcat(num2str(trange(i,1)),' of LF',' (hr)');
    c.Label.FontSize = 11;
    caxis(f.ax(2),[trange(i,2)/3600 trange(i,3)/3600])
    text(f.ax(2),0.85,0.15,'LF','FontSize',12,'unit','normalized');
    text(f.ax(2),0.04,0.93,'b','FontSize',11,'unit','normalized','EdgeColor','k','Margin',2);
    f.ax(2).Box = 'on';
    grid(f.ax(2), 'on');
    axis(f.ax(2), 'equal');
    f.ax(2).GridLineStyle = '--';
    f.ax(2).XAxisLocation = 'top';
    medx = median(miglf(:,1));
    medy = median(miglf(:,2));
    [rotx, roty] = complex_rot(0,5,-angbest(i));
    xvect = [medx-rotx medx+rotx];
    yvect = [medy-roty medy+roty];
    drawArrow(f.ax(2),xvect,yvect,xran,yran,'linewidth',1.5);
    [rotx, roty] = complex_rot(-5,0,-angbest(i));
    xvect = [medx-rotx medx+rotx];
    yvect = [medy-roty medy+roty];
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
    hold(f.ax(2),'off');
    
    % subplot 3 of figure i
    hold(f.ax(3),'on');
    f.ax(3).FontSize = 9;
    scatter(f.ax(3),mighfdum(:,15)/3600,mighfdum(:,1),20,[1 0 0],'filled','o',...
            'MarkerEdgeColor','k');   
    scatter(f.ax(3),miglfdum(:,15)/3600,miglfdum(:,1),20,[0.6 1 1],'filled','o',...
            'MarkerEdgeColor','k');
    
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
    hold(f.ax(3),'off');
    
    % subplot 4 of figure i
    hold(f.ax(4),'on');
    f.ax(4).FontSize = 9;
    scatter(f.ax(4),mighfdum(:,15)/3600,mighfdum(:,2),20,[1 0 0],'filled','o',...
            'MarkerEdgeColor','k');   
    scatter(f.ax(4),miglfdum(:,15)/3600,miglfdum(:,2),20,[0.6 1 1],'filled','o',...
            'MarkerEdgeColor','k');
    
    f.ax(4).Box = 'on';
    grid(f.ax(4), 'on');
    f.ax(4).GridLineStyle = '--';
    xran1 = [trange(i,2)/3600 trange(i,3)/3600];
    yran1 = [round(min([mighfdum(:,2);miglfdum(:,2)]))-1 ...
             round(max([mighfdum(:,2);miglfdum(:,2)]))+1];
    xlim(f.ax(4),xran1);
    ylim(f.ax(4),yran1);
    text(f.ax(4),0.04,0.91,'d','FontSize',11,'unit','normalized','EdgeColor','k','Margin',2);
    xlabel(f.ax(4),'Time (hr)','fontsize',11);
    ylabel(f.ax(4),'Dist. along ort. (km)','fontsize',11); 
    hold(f.ax(4),'off');
    
    

    
    %%% save figure
    print(f.fig,'-dpdf',strcat(rstpath,'/LZB.mig.proj.ortho',num2str(trange(i,1)),...
          '_',num2str(trange(i,2)),'-',num2str(trange(i,3)),'.',num2str(winlenhf),'_',...
          num2str(winlenlf),'.',num2str(lofflf),'.',num2str(ccminlf),'.pdf'));
end


%% check the relationship of contemp. HF and LF in the target circle
%%% ONLY deal with the detection catalog composed by eligible migrations, rather than the entire
%%% catalog, cause the false positives or random detecitons may be misleading to some extent

%%% index of migrations that have considerable detections inside the target gray circle through 
%%% visual check
indvis = [8,9,10,11,14,15,17,18,20,28,29,30,31,33,34,35,37,38,40,41,42,43,44,47];
indvis = indvis';

%%% now use an automatic way to do it see if they are the same
indauto = [];

% use the gray circle as the boundary
[x,y] = circle_chao(-9,2,5,1);
bnd0 = [x' y'];

% use the gray half circle as the boundary
[x,y] = sector_chao(-9,2,5,300,480,1);
x = [x x(1)];
y = [y y(1)];
bnd1 = [x' y'];     % left half, down-dip half

[x,y] = sector_chao(-9,2,5,120,300,1);
x = [x x(1)];
y = [y y(1)];
bnd2 = [x' y'];

numhfbnd = zeros(length(trange),1);
for i = 1: length(trange)
%     i=33;    
%     disp(trange(i,:));
    
    % detections of that migration
    indhf = find(hftime(:,13)==trange(i,1) & hftime(:,15)>=trange(i,2) & ...
                 hftime(:,15)<=trange(i,3));
    mighf = hftime(indhf,:);
    
    % detections inside bound
    [is,ion] = inpolygon(mighf(:,1),mighf(:,2),bnd0(:,1),bnd0(:,2));
    isinbnd = is | ion;
    mighfbnd0 = mighf(isinbnd == 1, :);
    
    numhfbnd(i) = size(mighfbnd0,1);
    if numhfbnd(i) > 10
        indauto = [indauto; i]; 
    end
    
end

% index that can be splitted with spatial division
indspat = [];

% index that can be splitted with temporal divison
indtemp = [11; 34; 38; 40; 44];
timediv = [15.37; 6.58; 10.03; 0.79; 1.1 ];

% index that need no splitting
indnosp = [8; 9; 10; 14; 15; 17; 28; 29; 30; 31; 33; 41; 42; 43; 45];

% index that need to throw away some time, in fact applying a circle to 35 is misleading
indsave = [18; 35; 47];
timesave = [4.421; 11.45; 2.78];

% index that need to throw some farther lf detections
indrand = [20; 37];

%%% Now question where the contemp LF detections are for the HF detections inside the gray circle,
%%% and analyze their relationship
for i = 8: length(indauto)
%     i=1;    
%     disp(trange(i,:));
    
    % detections of that migration
    indhf = find(hftime(:,13)==trange(indauto(i),1) & hftime(:,15)>=trange(indauto(i),2) & ...
                 hftime(:,15)<=trange(indauto(i),3));
    mighf = hftime(indhf,:);
    
    indlf = find(lftime(:,13)==trange(indauto(i),1) & lftime(:,15)>=trange(indauto(i),2) & ...
                 lftime(:,15)<=trange(indauto(i),3));
    miglf = lftime(indlf,:); 

    % detections inside bound 0
    [is,ion] = inpolygon(mighf(:,1),mighf(:,2),bnd0(:,1),bnd0(:,2));
    isinbnd0 = is | ion;
    mighfbnd0 = mighf(isinbnd0 == 1, :);
    
    miglfbnd0 = miglf(miglf(:,15)>=mighfbnd0(1,15) & miglf(:,15)<=mighfbnd0(end,15), :);    
    
    % detections inside bound 1
    [is,ion] = inpolygon(mighf(:,1),mighf(:,2),bnd1(:,1),bnd1(:,2));
    isinbnd1 = is | ion;
    mighfbnd1 = mighf(isinbnd1 == 1, :);
    miglfbnd1 = miglf(miglf(:,15)>=mighfbnd1(1,15) & miglf(:,15)<=mighfbnd1(end,15), :);
    
    % detections inside bound 2
    [is,ion] = inpolygon(mighf(:,1),mighf(:,2),bnd2(:,1),bnd2(:,2));
    isinbnd2 = is | ion;
    mighfbnd2 = mighf(isinbnd2 == 1, :);
    miglfbnd2 = miglf(miglf(:,15)>=mighfbnd2(1,15) & miglf(:,15)<=mighfbnd2(end,15), :);
        
    %%% now it is case by case
    % deal with mig that split with time
    if ismember(indauto(i), indtemp)
        [~,~,ii] = intersect(indauto(i), indtemp);
        time = timediv(ii);
        mighfbnd2 = mighfbnd0(mighfbnd0(:,15)/3600<time, :);
        miglfbnd2 = miglfbnd0(miglfbnd0(:,15)/3600<time, :);        
        mighfbnd1 = mighfbnd0(mighfbnd0(:,15)/3600>=time, :);
        miglfbnd1 = miglfbnd0(miglfbnd0(:,15)/3600>=time, :);
    end
    
    % deal with mig that throw some detections away
    if ismember(indauto(i), indsave)
        [~,~,ii] = intersect(indauto(i), indsave);
        time = timesave(ii);
        mighfbnd0 = mighfbnd0(mighfbnd0(:,15)/3600<time, :);
        miglfbnd0 = miglfbnd0(miglfbnd0(:,15)/3600<time, :);
    end

    % deal with mig that has some related outliers in the lower right quadrant    
    if ismember(indauto(i), indrand)
        miglfbnd0(miglfbnd0(:,1)>0 & miglfbnd0(:,2)<-5, :) = [];
    end
    
    % although it might be not very useful, but let us get the median of all detections inside
    medhfbnd0 = median(mighfbnd0(:,1:2),1);
    medlfbnd0 = median(miglfbnd0(:,1:2),1);
    medhfbnd1 = median(mighfbnd1(:,1:2),1);
    medlfbnd1 = median(miglfbnd1(:,1:2),1);
    medhfbnd2 = median(mighfbnd2(:,1:2),1);
    medlfbnd2 = median(miglfbnd2(:,1:2),1);

    % relative coordinates
    mighfdum = mighfbnd0;
    for j = 1: size(mighfbnd0,1)
        x0 = mighfdum(j,1);
        y0 = mighfdum(j,2);
        [newx,newy] = coordinate_rot(x0,y0,-(angbest(indauto(i))-90),0,0);
        mighfdum(j,1) = newx;
        mighfdum(j,2) = newy;
    end
    
    miglfdum = miglfbnd0;
    for j = 1: size(miglfbnd0,1)
        x0 = miglfdum(j,1);
        y0 = miglfdum(j,2);
        [newx,newy] = coordinate_rot(x0,y0,-(angbest(indauto(i))-90),0,0);
        miglfdum(j,1) = newx;
        miglfdum(j,2) = newy;
    end
    
    %%% define and position the figure frame and axes of each plot
    f.fig=figure;
    widin = 8;  % maximum width allowed is 8.5 inches
    htin = 9;   % maximum height allowed is 11 inches
    set(f.fig,'Position',[scrsz(1)+1*scrsz(3)/10 scrsz(2)+scrsz(4)/20 widin*res htin*res]);
    nrow = 3;
    ncol = 2;    
    for isub = 1:nrow*ncol
        f.ax(isub) = subplot(nrow,ncol,isub);
    end

    %%% reposition
    set(f.ax(1), 'position', [ 0.1, 0.64, 0.32, 0.32]);
    set(f.ax(2), 'position', [ 0.5, 0.64, 0.32, 0.32]);
    set(f.ax(3), 'position', [ 0.1, 0.35, 0.32, 0.22]);
    set(f.ax(4), 'position', [ 0.5, 0.35, 0.32, 0.22]);
    set(f.ax(5), 'position', [ 0.1, 0.078, 0.32, 0.22]);
    set(f.ax(6), 'position', [ 0.5, 0.078, 0.32, 0.22]);
    
    % subplot 1 of figure i
    hold(f.ax(1),'on');
    plot(f.ax(1),[-100 100],[0 0],'k--');
    plot(f.ax(1),[0 0],[-100 100],'k--');
    f.ax(1).FontSize = 9;
    mighfbnd0 = sortrows(mighfbnd0,-15);
    scatter(f.ax(1),mighfbnd0(:,1),mighfbnd0(:,2), 20, mighfbnd0(:,15)/3600, 'filled','o');
%     mighfbnd1 = sortrows(mighfbnd1,-15);
%     scatter(f.ax(1),mighfbnd1(:,1),mighfbnd1(:,2), 20, mighfbnd1(:,15)/3600, 'filled','o');  %, 'MarkerEdgeColor', 'w')
%     scatter(f.ax(1),relacont(:,1),relacont(:,2),20,'k','o','LineWidth',1);
    
%     plot(f.ax(1),x,y,'linewidth',1,'linestyle','-','color',[0.6 0.6 0.6]);
    
    plot(f.ax(1),bnd0(:,1),bnd0(:,2),'linewidth',1,'linestyle','-','color',[0.6 0.6 0.6]);
    plot(f.ax(1),bnd1(:,1),bnd1(:,2),'linewidth',1,'linestyle','-','color',[0.6 0.6 0.6]);
    plot(f.ax(1),bnd2(:,1),bnd2(:,2),'linewidth',1,'linestyle','-','color',[0.6 0.6 0.6]);
    
    colormap(f.ax(1),'jet');
    c=colorbar(f.ax(1),'SouthOutside');
    pos = f.ax(1).Position;
    c.Position = [pos(1), pos(2)-0.018, pos(3), 0.02];
%     c.TickLabels=[];
    juldate = num2str(trange(indauto(i),1));
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
    caxis(f.ax(1),[trange(indauto(i),2)/3600 trange(indauto(i),3)/3600])
    text(f.ax(1),0.85,0.15,'HF','FontSize',12,'unit','normalized');
    text(f.ax(1),0.04,0.93,'a','FontSize',11,'unit','normalized','EdgeColor','k','Margin',2);
    f.ax(1).Box = 'on';
    grid(f.ax(1), 'on');
    axis(f.ax(1), 'equal');
    f.ax(1).GridLineStyle = '--';
    f.ax(1).XAxisLocation = 'top';
    medx = median(mighfbnd1(:,1));
    medy = median(mighfbnd1(:,2));
    [rotx, roty] = complex_rot(0,5,-angbest(indauto(i)));
    xvect = [medx-rotx medx+rotx];
    yvect = [medy-roty medy+roty];
    drawArrow(f.ax(1),xvect,yvect,xran,yran,'linewidth',1);
    [rotx, roty] = complex_rot(-5,0,-angbest(indauto(i)));
    xvect = [medx-rotx medx+rotx];
    yvect = [medy-roty medy+roty];
    drawArrow(f.ax(1),xvect,yvect,xran,yran,'linewidth',1,'linestyle','--','color',[0.4 0.4 0.4]);
    xticks(f.ax(1),xran(1):5:xran(2));
    yticks(f.ax(1),yran(1):5:yran(2));    
    xlabel(f.ax(1),'E (km)','fontsize',11);
    ylabel(f.ax(1),'N (km)','fontsize',11);
    text(f.ax(1),0.5,0.1,num2str(indauto(i)),'FontSize',12,'unit','normalized');
    text(f.ax(1),0.83,0.95,strcat(num2str(size(mighfbnd1,1)),{' detections'}),'FontSize',8,...
         'unit','normalized','horizontalalignment','center');
    text(f.ax(1),0.83,0.90,strcat({'in boundary'}),'FontSize',8,'unit','normalized',...
         'horizontalalignment','center');
    
    scatter(f.ax(1),medhfbnd0(1),medhfbnd0(2),40,'k','p','LineWidth',1);
    scatter(f.ax(1),medhfbnd1(1),medhfbnd1(2),40,'k','s','LineWidth',1);
    scatter(f.ax(1),medhfbnd2(1),medhfbnd2(2),40,'k','d','LineWidth',1);

    hold(f.ax(1),'off');

    % subplot 2 of figure indauto(i)
    hold(f.ax(2),'on');
    plot(f.ax(2),[-100 100],[0 0],'k--');
    plot(f.ax(2),[0 0],[-100 100],'k--');
    f.ax(2).FontSize = 9;
    miglfbnd0 = sortrows(miglfbnd0,-15);
    scatter(f.ax(2),miglfbnd0(:,1),miglfbnd0(:,2), 20, miglfbnd0(:,15)/3600, 'filled','o');  %, 'MarkerEdgeColor', 'w')
%     scatter(f.ax(2),relacont(:,1),relacont(:,2),20,'k','o','LineWidth',1);

%     plot(f.ax(2),x,y,'linewidth',1,'linestyle','-','color',[0.6 0.6 0.6]);
    
    plot(f.ax(2),bnd0(:,1),bnd0(:,2),'linewidth',1,'linestyle','-','color',[0.6 0.6 0.6]);
    plot(f.ax(2),bnd1(:,1),bnd1(:,2),'linewidth',1,'linestyle','-','color',[0.6 0.6 0.6]);
    plot(f.ax(2),bnd2(:,1),bnd2(:,2),'linewidth',1,'linestyle','-','color',[0.6 0.6 0.6]);

    colormap(f.ax(2),'jet');
    c=colorbar(f.ax(2),'SouthOutside');
    pos = f.ax(2).Position;
    c.Position = [pos(1), pos(2)-0.018, pos(3), 0.02];
%     c.TickLabels=[];
    c.Label.String = strcat({'Time (hr) on '}, day, mo, yr);
%     c.Label.String = strcat(num2str(trange(indauto(i),1)),' of LF',' (hr)');
    c.Label.FontSize = 11;
    caxis(f.ax(2),[trange(indauto(i),2)/3600 trange(indauto(i),3)/3600])
    text(f.ax(2),0.85,0.15,'LF','FontSize',12,'unit','normalized');
    text(f.ax(2),0.04,0.93,'b','FontSize',11,'unit','normalized','EdgeColor','k','Margin',2);
    f.ax(2).Box = 'on';
    grid(f.ax(2), 'on');
    axis(f.ax(2), 'equal');
    f.ax(2).GridLineStyle = '--';
    f.ax(2).XAxisLocation = 'top';
    medx = median(miglfbnd1(:,1));
    medy = median(miglfbnd1(:,2));
    [rotx, roty] = complex_rot(0,5,-angbest(indauto(i)));
    xvect = [medx-rotx medx+rotx];
    yvect = [medy-roty medy+roty];
    drawArrow(f.ax(2),xvect,yvect,xran,yran,'linewidth',1);
    [rotx, roty] = complex_rot(-5,0,-angbest(indauto(i)));
    xvect = [medx-rotx medx+rotx];
    yvect = [medy-roty medy+roty];
    drawArrow(f.ax(2),xvect,yvect,xran,yran,'linewidth',1,'linestyle','--','color',[0.4 0.4 0.4]);
    xticks(f.ax(2),xran(1):5:xran(2));
    yticks(f.ax(2),yran(1):5:yran(2));    
    xlabel(f.ax(2),'E (km)','fontsize',11);
    ylabel(f.ax(2),'N (km)','fontsize',11);
    text(f.ax(2),0.5,0.1,num2str(indauto(i)),'FontSize',12,'unit','normalized');
    text(f.ax(2),0.83,0.95,strcat(num2str(size(miglfbnd1,1)),{' detections'}),'FontSize',8,...
         'unit','normalized','horizontalalignment','center');
    text(f.ax(2),0.83,0.90,strcat({'contemp'}),'FontSize',8,'unit','normalized',...
         'horizontalalignment','center');
        
    scatter(f.ax(2),medlfbnd0(1),medlfbnd0(2),40,'k','p','LineWidth',1);
    scatter(f.ax(2),medlfbnd1(1),medlfbnd1(2),40,'k','s','LineWidth',1);
    scatter(f.ax(2),medlfbnd2(1),medlfbnd2(2),40,'k','d','LineWidth',1);
    
    xvect = [15 15+medlfbnd0(1)-medhfbnd0(1)];
    yvect = [5  5+medlfbnd0(2)-medhfbnd0(2)];
    scatter(f.ax(2),15,5,40,'k','.');
    plot(f.ax(2),xvect,yvect,'linewidth',1,'linestyle','-','color',[0.55 0.55 0.55]);
    
    xvect = [12 12+medlfbnd2(1)-medhfbnd2(1)];
    yvect = [5  5+medlfbnd2(2)-medhfbnd2(2)];
    scatter(f.ax(2),12,5,40,'k','.');
    plot(f.ax(2),xvect,yvect,'linewidth',1,'linestyle','-','color',[0.55 0.55 0.55]);
    
    xvect = [18 18+medlfbnd1(1)-medhfbnd1(1)];
    yvect = [5  5+medlfbnd1(2)-medhfbnd1(2)];
    scatter(f.ax(2),18,5,40,'k','.');
    plot(f.ax(2),xvect,yvect,'linewidth',1,'linestyle','-','color',[0.55 0.55 0.55]);
    
    hold(f.ax(2),'off');
    
    % subplot 3 of figure indauto(i)
    hold(f.ax(3),'on');
    f.ax(3).FontSize = 9;
    scatter(f.ax(3),mighfdum(:,15)/3600,mighfdum(:,1),20,[1 0 0],'filled','o',...
            'MarkerEdgeColor','k');   
    scatter(f.ax(3),miglfdum(:,15)/3600,miglfdum(:,1),20,[0.6 1 1],'filled','o',...
            'MarkerEdgeColor','k');
    if ~isempty(mighfbnd1)
        plot(f.ax(3),[mighfbnd1(1,15)/3600 mighfbnd1(1,15)/3600],[-100 100],'r--');
        plot(f.ax(3),[mighfbnd1(end,15)/3600 mighfbnd1(end,15)/3600],[-100 100],'r--');
    end
    if ~isempty(mighfbnd2)
        plot(f.ax(3),[mighfbnd2(1,15)/3600 mighfbnd2(1,15)/3600],[-100 100],'b--');
        plot(f.ax(3),[mighfbnd2(end,15)/3600 mighfbnd2(end,15)/3600],[-100 100],'b--');
    end
    f.ax(3).Box = 'on';
    grid(f.ax(3), 'on');
    f.ax(3).GridLineStyle = '--';
    xran1 = [trange(indauto(i),2)/3600 trange(indauto(i),3)/3600];
    yran1 = [round(min([mighfdum(:,1);miglfdum(:,1)]))-1 ...
             round(max([mighfdum(:,1);miglfdum(:,1)]))+1];
    xlim(f.ax(3),xran1);
    ylim(f.ax(3),yran1);
    text(f.ax(3),0.04,0.91,'c','FontSize',11,'unit','normalized','EdgeColor','k','Margin',2);
    xlabel(f.ax(3),'Time (hr)','fontsize',11);
    ylabel(f.ax(3),'Dist. along prop. (km)','fontsize',11); 
    hold(f.ax(3),'off');
    
    % subplot 4 of figure indauto(i)
    hold(f.ax(4),'on');
    f.ax(4).FontSize = 9;
    scatter(f.ax(4),mighfdum(:,15)/3600,mighfdum(:,2),20,[1 0 0],'filled','o',...
            'MarkerEdgeColor','k');   
    scatter(f.ax(4),miglfdum(:,15)/3600,miglfdum(:,2),20,[0.6 1 1],'filled','o',...
            'MarkerEdgeColor','k');
    if ~isempty(miglfbnd1)
        plot(f.ax(4),[miglfbnd1(1,15)/3600 miglfbnd1(1,15)/3600],[-100 100],'r--');
        plot(f.ax(4),[miglfbnd1(end,15)/3600 miglfbnd1(end,15)/3600],[-100 100],'r--');
    end
    if ~isempty(miglfbnd2)
        plot(f.ax(4),[miglfbnd2(1,15)/3600 miglfbnd2(1,15)/3600],[-100 100],'b--');
        plot(f.ax(4),[miglfbnd2(end,15)/3600 miglfbnd2(end,15)/3600],[-100 100],'b--');
    end
    f.ax(4).Box = 'on';
    grid(f.ax(4), 'on');
    f.ax(4).GridLineStyle = '--';
    xran1 = [trange(indauto(i),2)/3600 trange(indauto(i),3)/3600];
    yran1 = [round(min([mighfdum(:,2);miglfdum(:,2)]))-1 ...
             round(max([mighfdum(:,2);miglfdum(:,2)]))+1];
    xlim(f.ax(4),xran1);
    ylim(f.ax(4),yran1);
    text(f.ax(4),0.04,0.91,'d','FontSize',11,'unit','normalized','EdgeColor','k','Margin',2);
    xlabel(f.ax(4),'Time (hr)','fontsize',11);
    ylabel(f.ax(4),'Dist. along ort. (km)','fontsize',11); 
    hold(f.ax(4),'off');
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% It seems that median of the previous detections is a promising way to show offsets as well,
    %%% so try to deal with it the same way
    for j = 1: size(mighfdum,1)
        medhf(j,1) = median(mighfdum(1:j,1));
        medhf(j,2) = median(mighfdum(1:j,2));
        medhf(j,3) = mighfdum(j,15);
    end
    
    for j = 1: size(miglfdum,1)
        medlf(j,1) = median(miglfdum(1:j,1));
        medlf(j,2) = median(miglfdum(1:j,2));
        medlf(j,3) = miglfdum(j,15);
    end
       
    
    % subplot 5 of figure i
    hold(f.ax(5),'on');
    f.ax(5).FontSize = 9;
    scatter(f.ax(5),medhf(:,3)/3600,medhf(:,1),20,[1 0 0],'filled','o','MarkerEdgeColor',...
            'k');
    scatter(f.ax(5),medlf(:,3)/3600,medlf(:,1),20,[0.6 1 1],'filled','o','MarkerEdgeColor',...
            'k');            
    f.ax(5).Box = 'on';
    grid(f.ax(5), 'on');
    f.ax(5).GridLineStyle = '--';
    xran1 = [trange(indauto(i),2)/3600 trange(indauto(i),3)/3600];
    yran1 = [round(min([medhf(:,1);medlf(:,1)]))-1 ...
             round(max([medhf(:,1);medlf(:,1)]))+1];
    xlim(f.ax(5),xran1);
    ylim(f.ax(5),yran1);
    text(f.ax(5),0.97,0.1,'e','FontSize',11,'unit','normalized','EdgeColor','k','Margin',2,...
         'HorizontalAlignment','right');
    xlabel(f.ax(5),'Time (hr)','fontsize',11);
    ylabel(f.ax(5),'Med. dist. along prop. (km)','fontsize',11);
    hold(f.ax(5),'off');
    
    % subplot 6 of figure indauto(i)
    hold(f.ax(6),'on');
    f.ax(6).FontSize = 9;
    scatter(f.ax(6),medhf(:,3)/3600,medhf(:,2),20,[1 0 0],'filled','o','MarkerEdgeColor',...
            'k');
    scatter(f.ax(6),medlf(:,3)/3600,medlf(:,2),20,[0.6 1 1],'filled','o','MarkerEdgeColor',...
            'k');          
    f.ax(6).Box = 'on';
    grid(f.ax(6), 'on');
    f.ax(6).GridLineStyle = '--';
    xran1 = [trange(indauto(i),2)/3600 trange(indauto(i),3)/3600];
    yran1 = [round(min([medhf(:,2);medlf(:,2)]))-1 ...
             round(max([medhf(:,2);medlf(:,2)]))+1];
    xlim(f.ax(6),xran1);
    ylim(f.ax(6),yran1);
    text(f.ax(6),0.97,0.1,'f','FontSize',11,'unit','normalized','EdgeColor','k','Margin',2,...
         'HorizontalAlignment','right');
    xlabel(f.ax(6),'Time (hr)','fontsize',11);
    ylabel(f.ax(6),'Med. dist. along ort. (km)','fontsize',11);
    hold(f.ax(6),'off');

    
    %%% save figure
    print(f.fig,'-dpdf',strcat(rstpath,'/LZB.mig.proj.incircle',num2str(trange(indauto(i),1)),...
          '_',num2str(trange(indauto(i),2)),'-',num2str(trange(indauto(i),3)),'.',num2str(winlenhf),...
          '_',num2str(winlenlf),'.',num2str(lofflf),'.',num2str(ccminlf),'.pdf'));

end


%%% it seems that in most migrations, either along propagation or othogonal direction or both, shows
%%% the southern shift of LF detections relative to HF ones, with some exceptions;

% index of migrations that are almost centered 
indctr = [17,35,37,38,40,41];
indctr = indctr';

% index of migrations in which LF are more northly displaced
indnor = [43,44,45];
indnor = indnor';

















