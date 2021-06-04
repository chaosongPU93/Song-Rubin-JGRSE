% function plt_modtime_mig
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script is used to plot the migrations coming from the automated tremor
% bursts in 'identify_RTMs_v3', but need to be modified in time according to
% the spatial information.
%
%   
% NOTES:
%   2020/09/08, i am adding a new fam 006, so that i want to check if original
%               migrations are still reasonable
%   2020/09/10, i am adding a fam 001 back again for last try, 13 fam in total now
% 
% Chao Song, chaosong@princeton.edu
% First created date:   2021/02/01
% Last modified date:   2021/02/01
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Initialization
% format short e   % Set the format to 5-digit floating point
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
            '017';
            '006';
            '001';
            ];     % 2020/09/07, i am adding a new family 006 to the pool, final version

nfam = size(nfampool,1);
disp(nfam); 

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
           ];


relacont = loccont;
[dx, dy] = absloc2relaloc(loccont(:,1),loccont(:,2),loccont(2,1),loccont(2,2));
relacont(:,1) = dx;
relacont(:,2) = dy;


% sort according to day, sec
hftime = sortrows(hftime, [13, 15]);
lftime = sortrows(lftime, [13, 15]);

%%
randiv = [
    2004197   8.3402e+04   8.6659e+04
    2004197   8.3402e+04   8.4888e+04   % divided into 3
    2004197   8.4888e+04   8.5320e+04   % divided into 3
    2004197   8.5320e+04   8.6659e+04   % divided into 3
    ];

rancom = [
    2004199   4.6340e+04   4.7036e+04
    2004199   4.7227e+04   4.7633e+04
    2004199   4.7633e+04   4.9126e+04
    2004199   4.6340e+04   4.9126e+04   % combine, care speed direct, rmse, 120
    ];

ranmod = [
%     2004200   1.1706e+04   1.5799e+04
%     2004200   1.2600e+04   1.5799e+04   % modified time, y
%     2005256   8.0330e+03   1.0404e+04    
%     2005256   8.0330e+03   9.8280e+03   % change end time, y
%     2005256   2.0437e+04   2.2294e+04
%     2005256   2.1492e+04   2.2294e+04   % change start time, y
%     2005256   2.7719e+04   2.9172e+04
%     2005256   2.8080e+04   2.9172e+04   % change start time, y
%     2005256   3.4604e+04   3.5238e+04
%     2005256   3.4776e+04   3.5238e+04   % change start time, y
%     2005256   3.5597e+04   4.1112e+04
%     2005256   3.6900e+04   3.8520e+04   % change times, y
%     2005256   5.1394e+04   5.4053e+04
%     2005256   5.2020e+04   5.2596e+04   % change times, y
%     2005256   6.0723e+04   6.4865e+04
%     2005256   6.2640e+04   6.4865e+04   % change start time, y, check speed direc, y, but use rmse 130
%     2005256   7.6381e+04   7.8664e+04
%     2005256   7.6381e+04   7.7940e+04   % change end time, y
%     2005257   4.9440e+03   7.7460e+03
%     2005257   6.1200e+03   7.7460e+03   % change start time, y, but check speed direc, y
%     2005257   1.5973e+04   1.8655e+04
%     2005257   1.6920e+04   1.8655e+04   % change start time, y
%     2005257   3.8513e+04   4.5418e+04
%     2005257   3.9000e+04   4.5418e+04   % change start time, y
    2005257   7.3440e+04   8.2043e+04
    2005257   7.3440e+04   7.8840e+04   % change end time, y, but check speed direc, y, but use rmse 235
    ];



%%
trange = ranmod;

xran = [-15 25];
yran = [-20 20];

angbestl1 = zeros(size(trange,1),1);
angbestl2 = zeros(size(trange,1),1);
angbestl3 = zeros(size(trange,1),1);
angbestl4 = zeros(size(trange,1),1);
angbestl5 = zeros(size(trange,1),1);

resprophf = nan(size(trange,1)+1,200);
resproplf = nan(size(trange,1)+1,50);

% create fit object with free constraints
fttpfree = fittype( @(a,b,x) a*x+b);

indcheck = 1: size(trange,1);
indplt = indcheck;

angle = 0:5:355;

%     l1normhf = zeros(length(angle),1);
%     l2normhf = zeros(length(angle),1);
slopehf = zeros(size(trange,1),length(angle));
%     ssehf = zeros(length(angle),1);
rmsehf = zeros(size(trange,1),length(angle));
rsquarehf = zeros(size(trange,1),length(angle));

%     l1normlf = zeros(length(angle),1);
%     l2normlf = zeros(length(angle),1);
slopelf = zeros(size(trange,1),length(angle));
%     sself = zeros(length(angle),1);
rmself = zeros(size(trange,1),length(angle));
rsquarelf = zeros(size(trange,1),length(angle));
    
for i = 1: length(indplt)
% for i = 41: 63
    disp(trange(indplt(i),:));
    indhf = find(hftime(:,13)==trange(indplt(i),1) & hftime(:,15)>=trange(indplt(i),2) & ...
                 hftime(:,15)<=trange(indplt(i),3));
    mighf = hftime(indhf,:);

    indlf = find(lftime(:,13)==trange(indplt(i),1) & lftime(:,15)>=trange(indplt(i),2) & ...
                 lftime(:,15)<=trange(indplt(i),3));
    miglf = lftime(indlf,:);
    
    if size(mighf,1)< 15 || size(miglf,1) < 15
        disp(strcat(num2str(indplt(i)),' not enough detections'));
    end
% end

    
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
        slopehf(indplt(i),iang) = coef(1);
        fitprophf = feval(fitobj,mighfdum(:,15)/3600);
%         l1normhf(iang) = sum(abs(mighfdum(:,1)-fitprophf))/(length(mighfdum(:,1)));
%         l2normhf(iang) = sum((mighfdum(:,1)-fitprophf).^2)/(length(mighfdum(:,1)));
%         ssehf(iang) = gof.sse;
        rmsehf(indplt(i),iang) = gof.rmse;
        rsquarehf(indplt(i),iang) = gof.rsquare;    % r-square, defined as ratio of the sum of squares of the regression and the total sum of squares
    
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
        slopelf(indplt(i),iang) = coef(1);
        fitproplf = feval(fitobj,miglfdum(:,15)/3600);
%         l1normlf(iang) = sum(abs(miglfdum(:,1)-fitproplf))/(length(miglfdum(:,1)));
%         l2normlf(iang) = sum((miglfdum(:,1)-fitproplf).^2)/(length(miglfdum(:,1)));
%         sself(iang) = gof.sse;
        rmself(indplt(i),iang) = gof.rmse;
        rsquarelf(indplt(i),iang) = gof.rsquare;    % r-square, defined as ratio of the sum of squares of the regression and the total sum of squares
        
    end

    %%% best angle estimate from hf
    ind = find(slopehf(indplt(i),:)>0);
%     ind1 = find(l1normhf(ind)==min(l1normhf(ind)));
%     angbestl1(indplt(i)) = angle(ind(ind1(1)));
%     ind2 = find(l2normhf(ind)==min(l2normhf(ind)));
%     angbestl2(indplt(i)) = angle(ind(ind2(1)));
    ind3 = find(rmsehf(indplt(i),ind)==min(rmsehf(indplt(i),ind)));     % rmse, Root Mean Squared Error, the smaller, the better
    if length(ind3) > 1
        disp(strcat({'multiple angles have same rmse: '},num2str(angle(ind3))));
    end
    angrmsehf(indplt(i)) = angle(ind(ind3(1)));
%     ind4 = find(rsquarehf(ind)==max(rsquarehf(ind)));   % R-square, 0-1, the larger the better
%     angbestl4(indplt(i)) = angle(ind(ind4(1)));
%     ind5 = find(ssehf(ind)==min(ssehf(ind)));   % sum of square due to error, the smaller, the better
%     angbestl5(indplt(i)) = angle(ind(ind5(1)));
%     disp([angbestl3(indplt(i)),angbestl4(indplt(i)),angbestl5(indplt(i))])
%     angbesthf(indplt(i)) = median([angbestl3(indplt(i)),angbestl4(indplt(i)),angbestl5(indplt(i))]);

    ind6 = find(slopehf(indplt(i),:)==max(slopehf(indplt(i),:))); % one with the largest slope, i.e., migrating speed
    if length(ind6) > 1
        disp(strcat({'multiple angles have same slope:'},num2str(angle(ind6))));
    end
    angslopehf(indplt(i)) = angle(ind6(1));
    
    
    %%% best angle estimate from hf
    ind = find(slopelf(indplt(i),:)>0);
%     ind1 = find(l1normlf(ind)==min(l1normlf(ind)));
%     angbestl1(indplt(i)) = angle(ind(ind1(1)));
%     ind2 = find(l2normlf(ind)==min(l2normlf(ind)));
%     angbestl2(indplt(i)) = angle(ind(ind2(1)));
    ind3 = find(rmself(indplt(i),ind)==min(rmself(indplt(i),ind)));     % rmse, Root Mean Squared Error, the smaller, the better
    angrmself(indplt(i)) = angle(ind(ind3(1)));
%     ind4 = find(rsquarelf(ind)==max(rsquarelf(ind)));   % R-square, 0-1, the larger the better
%     angbestl4(indplt(i)) = angle(ind(ind4(1)));
%     ind5 = find(sself(ind)==min(sself(ind)));   % sum of square due to error, the smaller, the better
%     angbestl5(indplt(i)) = angle(ind(ind5(1)));
%     disp([angbestl3(indplt(i)),angbestl4(indplt(i)),angbestl5(indplt(i))])
%     angbestlf(indplt(i)) = median([angbestl3(indplt(i)),angbestl4(indplt(i)),angbestl5(indplt(i))]);

    ind6 = find(slopelf(indplt(i),:)==max(slopelf(indplt(i),:))); % one with the largest slope, i.e., migrating speed
    if length(ind6) > 1
        disp(strcat({'multiple angles have same slope:'},num2str(angle(ind6))));
    end
    angslopelf(indplt(i)) = angle(ind6(1));
    
    % determine if more inspection is needed
    if abs(angrmsehf(indplt(i))-angrmself(indplt(i))) > 20
        disp('Difference in propagation direction between hf & lf is too large!');
    end
    
end

% angbesthf = angbesthf';
% angbestlf = angbestlf';
% angbest = [angbesthf angbestlf];

angbest1 = angrmsehf';
angbest2 = angslopehf';

% angbest1 = angrmself';
% angbest2 = angslopelf';
% angbest1 = [245];

angcandihf = [angbest1 angbest2];


%%%
indspeed = [
     2
    27
    34
    37
    47
    ];

indspeedck = [];

% create fit object with free constraints
fttpfree = fittype( @(a,b,x) a*x+b);

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
        [newx,newy] = coordinate_rot(x0,y0,-(angbest1(i)-90),0,0);
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
        [newx,newy] = coordinate_rot(x0,y0,-(angbest1(i)-90),0,0);
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
    text(f.ax(1),0.04,0.1,'TWKB','FontSize',11,'unit','normalized','EdgeColor','k','Margin',2);
    f.ax(1).Box = 'on';
    grid(f.ax(1), 'on');
    axis(f.ax(1), 'equal');
    f.ax(1).GridLineStyle = '--';
    f.ax(1).XAxisLocation = 'top';
    medxhf = median(mighf(:,1));
    medyhf = median(mighf(:,2));
    [rotx, roty] = complex_rot(0,5,-angbest1(i));
    xvect = [medxhf-rotx medxhf+rotx];
    yvect = [medyhf-roty medyhf+roty];
    drawArrow(f.ax(1),xvect,yvect,xran,yran,'linewidth',1.5);
    [rotx, roty] = complex_rot(0,5,-angbest2(i));
    xvect = [medxhf-rotx medxhf+rotx];
    yvect = [medyhf-roty medyhf+roty];
    drawArrow(f.ax(1),xvect,yvect,xran,yran,'linewidth',1.5,'color','b');
    [rotx, roty] = complex_rot(-5,0,-angbest1(i));
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
    text(f.ax(1),0.83,0.75,strcat(num2str(angbest1(i)),{' ^o'}),'FontSize',10,...
         'unit','normalized','horizontalalignment','center');
    text(f.ax(1),0.83,0.65,strcat(num2str(angbest2(i)),{' ^o'}),'FontSize',10,...
         'unit','normalized','horizontalalignment','center','color','b'); 
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
    text(f.ax(2),0.04,0.1,'TWKB','FontSize',11,'unit','normalized','EdgeColor','k','Margin',2);
    f.ax(2).Box = 'on';
    grid(f.ax(2), 'on');
    axis(f.ax(2), 'equal');
    f.ax(2).GridLineStyle = '--';
    f.ax(2).XAxisLocation = 'top';
    medxlf = median(miglf(:,1));
    medylf = median(miglf(:,2));
    [rotx, roty] = complex_rot(0,5,-angbest1(i));
    xvect = [medxlf-rotx medxlf+rotx];
    yvect = [medylf-roty medylf+roty];
    drawArrow(f.ax(2),xvect,yvect,xran,yran,'linewidth',1.5);
    [rotx, roty] = complex_rot(0,5,-angbest2(i));
    xvect = [medxlf-rotx medxlf+rotx];
    yvect = [medylf-roty medylf+roty];
    drawArrow(f.ax(2),xvect,yvect,xran,yran,'linewidth',1.5,'color','b'); 
    [rotx, roty] = complex_rot(-5,0,-angbest1(i));
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
    text(f.ax(2),0.83,0.75,strcat(num2str(angbest1(i)),{' ^o'}),'FontSize',10,...
         'unit','normalized','horizontalalignment','center');
    text(f.ax(2),0.83,0.65,strcat(num2str(angbest2(i)),{' ^o'}),'FontSize',10,...
         'unit','normalized','horizontalalignment','center','color','b'); 
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
    offset1(i) = intcptproplf(i)-intcptprophf(i);
    plot(f.ax(3),miglfdum2(:,15)/3600,fitproplf,'-','linewidth',2,'color',[0.6 1 1]);
    f.ax(3).Box = 'on';
    grid(f.ax(3), 'on');
    f.ax(3).GridLineStyle = '--';
    xran1 = [trange(i,2)/3600 trange(i,3)/3600];
%     yran1 = [round(min([mighfdum(:,1);miglfdum(:,1)]))-1 ...
%              round(max([mighfdum(:,1);miglfdum(:,1)]))+1];
    yran1 = [round(prctile([mighfdum(:,1);miglfdum(:,1)], 2)) ...
             round(prctile([mighfdum(:,1);miglfdum(:,1)], 98))+5];     
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
    slopesehf = gofhfprop.rmse./sqrt(sum((x-mean(x)).^2));
    est = slopeprophf(i);
    slopeprophfCI(i,:) = confidence_interval_general(est,slopesehf,length(x)-2,95);
    interceptsehf = slopesehf.*sqrt(sum(x.^2)./length(x));
    est = intcptprophf(i);
    intcptprophfCI(i,:) = confidence_interval_general(est,interceptsehf,length(x)-2,95);
    sehf(i, :) = [slopesehf interceptsehf];
    
    x = miglfdum2(:,15)/3600;
    slopeself = goflfprop.rmse./sqrt(sum((x-mean(x)).^2));
    interceptself = slopeself.*sqrt(sum(x.^2)./length(x));
    self(i, :) = [goflfprop.rmse slopeself interceptself];
    
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
    
    text(f.ax(3),0.4,0.95,sprintf('Slope: %.1f km/h',slopeprophf(i)),'FontSize',8,...
         'unit','normalized','horizontalalignment','left');
    text(f.ax(3),0.98,0.95,sprintf('SE: %.2f',slopesehf),'FontSize',8,...
         'unit','normalized','horizontalalignment','right'); 
    text(f.ax(3),0.4,0.88,sprintf('Pearson: %.3f',pearwthf(i)),'FontSize',8,...
         'unit','normalized','horizontalalignment','left');
    text(f.ax(3),0.98,0.87,sprintf('SE_{LF}: %.2f',slopeself),'FontSize',8,...
         'unit','normalized','horizontalalignment','right'); 
    text(f.ax(3),0.4,0.8,sprintf('LF-HF: %.2f km',offset1(i)),'FontSize',10,'unit','normalized'); 
    hold(f.ax(3),'off');
    
    % subplot 4 of figure i
    hold(f.ax(4),'on');
    f.ax(4).FontSize = 9;
    f.ax(4).Box = 'on';
    yyaxis(f.ax(4),'left');
    plot(f.ax(4),angle,rmsehf(i,:),'o-','linew',1,'color','k','markers',4);
    if isempty(intersect(indspeed,i))
        plot(f.ax(4),[angbest1(i) angbest1(i)],f.ax(4).YLim,'--','linew',1.5,'color','k');
    else
        plot(f.ax(4),[angbest1(i) angbest1(i)],f.ax(4).YLim,'--','linew',0.5,'color','k');
    end
    f.ax(4).YColor = 'k';
    ylabel(f.ax(4),'RMSE of HF','fontsize',11);
    yyaxis(f.ax(4),'right');
    plot(f.ax(4),angle,slopehf(i,:),'^:','linew',1,'color',[0.12 0.56 1],'markers',4);
    if ~isempty(intersect(indspeed,i))
        plot(f.ax(4),[angbest2(i) angbest2(i)],f.ax(4).YLim,'--','linew',1.5,'color',[0.12 0.56 1]);
    else
        plot(f.ax(4),[angbest2(i) angbest2(i)],f.ax(4).YLim,'--','linew',0.5,'color',[0.12 0.56 1]);
    end
    f.ax(4).YColor = [0.12 0.56 1];
    ylabel(f.ax(4),'Slope of HF (km/h)','fontsize',11);
    
    xlim(f.ax(4),[0 360]);
%     ymax = f.ax(4).YLim(2)+0.1;
%     ylim(f.ax(4),[0 ymax]);
%     ylim(f.ax(4),[0 1]);
    xlabel(f.ax(4),'Trial propagation direction (^o)','fontsize',11);
    hold(f.ax(4),'off');
    
    % subplot 5 of figure i    
    mighfdum = mighf;
    for j = 1: size(mighf,1)
        x0 = mighfdum(j,1);
        y0 = mighfdum(j,2);
        [newx,newy] = coordinate_rot(x0,y0,-(angbest2(i)-90),0,0);
        mighfdum(j,1) = newx;
        mighfdum(j,2) = newy;
    end
        
    miglfdum = miglf;
    for j = 1: size(miglf,1)
        x0 = miglfdum(j,1);
        y0 = miglfdum(j,2);
        [newx,newy] = coordinate_rot(x0,y0,-(angbest2(i)-90),0,0);
        miglfdum(j,1) = newx;
        miglfdum(j,2) = newy;
    end
    
    hold(f.ax(5),'on');
    f.ax(5).FontSize = 9;
    scatter(f.ax(5),mighfdum(:,15)/3600,mighfdum(:,1),20,[0.6 0.6 0.6],'filled','o',...
            'MarkerEdgeColor','k');    
    scatter(f.ax(5),miglfdum(:,15)/3600,miglfdum(:,1),20,[0.6 1 1],'filled','o',...
            'MarkerEdgeColor','k');
        
        mintime = max(min(mighfdum(:,15)), min(miglfdum(:,15)));
        maxtime = min(max(mighfdum(:,15)), max(miglfdum(:,15)));
        mighfdum2 = mighfdum(mighfdum(:,15)>=mintime & mighfdum(:,15)<=maxtime, :);
        miglfdum2 = miglfdum(miglfdum(:,15)>=mintime & miglfdum(:,15)<=maxtime, :);

    [fitobjhfprop,gofhfprop,outphfprop] = fit(mighfdum2(:,15)/3600, mighfdum2(:,1),fttpfree,'Robust',...
                                              'Bisquare','StartPoint',[1 1]);
    % output fit parameters
    coefprop = coeffvalues(fitobjhfprop);
    slopeprophf(i) = coefprop(1);
    intcptprophf(i) = coefprop(2);
    fitprophf = feval(fitobjhfprop,mighfdum2(:,15)/3600);
    plot(f.ax(5),mighfdum2(:,15)/3600,fitprophf,'-','linewidth',2,'color',[0.6 0.6 0.6]);
        
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
    offset2(i) = intcptproplf(i)-intcptprophf(i);
    plot(f.ax(5),miglfdum2(:,15)/3600,fitproplf,'-','linewidth',2,'color',[0.6 1 1]);
    f.ax(5).Box = 'on';
    grid(f.ax(5), 'on');
    f.ax(5).GridLineStyle = '--';
    xran1 = [trange(i,2)/3600 trange(i,3)/3600];
%     yran1 = [round(min([mighfdum(:,1);miglfdum(:,1)]))-1 ...
%              round(max([mighfdum(:,1);miglfdum(:,1)]))+1];
    yran1 = [round(prctile([mighfdum(:,1);miglfdum(:,1)], 2)) ...
             round(prctile([mighfdum(:,1);miglfdum(:,1)], 98))+5];     
    xlim(f.ax(5),xran1);
    ylim(f.ax(5),yran1);
    text(f.ax(5),0.04,0.91,'e','FontSize',11,'unit','normalized','EdgeColor','k','Margin',2);
    xlabel(f.ax(5),'Time (hr)','fontsize',11);
    ylabel(f.ax(5),'Dist. along prop. (km)','fontsize',11); 
    
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
    slopesehf2 = gofhfprop.rmse./sqrt(sum((x-mean(x)).^2));
    est = slopeprophf(i);
    slopeprophfCI(i,:) = confidence_interval_general(est,slopesehf2,length(x)-2,95);
    interceptsehf2 = slopesehf2.*sqrt(sum(x.^2)./length(x));
    est = intcptprophf(i);
    intcptprophfCI(i,:) = confidence_interval_general(est,interceptsehf2,length(x)-2,95);
    sehf2(i, :) = [slopesehf2 interceptsehf2]; 
    
    x = miglfdum2(:,15)/3600;
    slopeself2 = goflfprop.rmse./sqrt(sum((x-mean(x)).^2));
    interceptself2 = slopeself2.*sqrt(sum(x.^2)./length(x));
    self2(i, :) = [goflfprop.rmse slopeself2 interceptself2];
    
    if ~isempty(intersect(indspeed,i))
        sehfbest(i, :) = sehf2(i, :);
        selfbest(i, :) = self2(i, :);
    else
        sehfbest(i, :) = sehf(i, :);
        selfbest(i, :) = self(i, :);  
    end
    
    x = mighfdum2(:,15)/3600;
    y = mighfdum2(:,1);
    x_bar = wt_mean(x,wthf);
    y_bar = wt_mean(y,wthf);
    x_var = sum(wthf.*(x-x_bar).^2) / sum(wthf);
    y_var = sum(wthf.*(y-y_bar).^2) / sum(wthf);
    xy_cov = sum(wthf.*(x-x_bar).*(y-y_bar)) / sum(wthf);
    pearwthf2(i) = xy_cov / sqrt(x_var*y_var);
    
    if i~= 49 && pearwthf2(i) - pearwthf(i) >= 0.05
        indspeedck = [indspeedck; i];
    end
    
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
    
    text(f.ax(5),0.4,0.95,sprintf('Slope: %.1f km/h',slopeprophf(i)),'FontSize',8,...
         'unit','normalized','horizontalalignment','left');
    text(f.ax(5),0.98,0.95,sprintf('SE: %.2f',slopesehf2),'FontSize',8,...
         'unit','normalized','horizontalalignment','right'); 
    text(f.ax(5),0.4,0.88,sprintf('Pearson: %.3f',pearwthf2(i)),'FontSize',8,...
         'unit','normalized','horizontalalignment','left');
    text(f.ax(5),0.98,0.87,sprintf('SE_{LF}: %.2f',slopeself2),'FontSize',8,...
         'unit','normalized','horizontalalignment','right');  
    text(f.ax(5),0.4,0.8,sprintf('LF-HF: %.2f km',offset2(i)),'FontSize',10,'unit','normalized'); 
    hold(f.ax(5),'off');
    



end


%% migration needs to divided into 3 parts
trange = randiv;

xran = [-20 25];
yran = [-20 20];

%%% define and position the figure frame and axes of each plot
f.fig=figure;
widin = 8;  % maximum width allowed is 8.5 inches
htin = 9;   % maximum height allowed is 11 inches
set(f.fig,'Position',[scrsz(1)+1*scrsz(3)/15 scrsz(2)+scrsz(4)/20 widin*res htin*res]);
nrow = 2;
ncol = 2;
for isub = 1:nrow*ncol
    f.ax(isub) = subplot(nrow,ncol,isub);
end

%%% reposition
set(f.ax(1), 'position', [ 0.08, 0.64, 0.36, 0.32]);
set(f.ax(2), 'position', [ 0.52, 0.64, 0.36, 0.32]);
set(f.ax(3), 'position', [ 0.08, 0.22, 0.36, 0.32]);
set(f.ax(4), 'position', [ 0.52, 0.22, 0.36, 0.32]);

    
for i = 1: size(trange,1)
% for i = 1: 211
%     i=5;
%     disp(trange(i,:));
    indhf = find(hftime(:,13)==trange(i,1) & hftime(:,15)>=trange(i,2) & ...
                 hftime(:,15)<=trange(i,3));
    mighf = hftime(indhf,:);
        
        
    % subplot 1 of figure i
    ax = f.ax(i);
    hold(ax,'on');
    plot(ax,[-100 100],[0 0],'k--');
    plot(ax,[0 0],[-100 100],'k--');
    ax.FontSize = 9;
    mighf = sortrows(mighf,-15);
    scatter(ax,mighf(:,1),mighf(:,2), 20, mighf(:,15)/3600, 'filled','o');  %, 'MarkerEdgeColor', 'w')
    colormap(ax,'jet');
    c=colorbar(ax,'SouthOutside');
    pos = ax.Position;
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
    caxis(ax,[trange(i,2)/3600 trange(i,3)/3600])
    text(ax,0.85,0.1,'HF','FontSize',12,'unit','normalized');
%     text(ax,0.04,0.93,'a','FontSize',11,'unit','normalized','EdgeColor','k','Margin',2);
    text(ax,0.04,0.1,'TWKB','FontSize',11,'unit','normalized','EdgeColor','k','Margin',2);
    ax.Box = 'on';
    grid(ax, 'on');
    axis(ax, 'equal');
    ax.GridLineStyle = '--';
    ax.XAxisLocation = 'top';
%     medxhf = median(mighf(:,1));
%     medyhf = median(mighf(:,2));
%     [rotx, roty] = complex_rot(0,5,-angbest1(i));
%     xvect = [medxhf-rotx medxhf+rotx];
%     yvect = [medyhf-roty medyhf+roty];
%     drawArrow(ax,xvect,yvect,xran,yran,'linewidth',1.5);
%     [rotx, roty] = complex_rot(0,5,-angbest2(i));
%     xvect = [medxhf-rotx medxhf+rotx];
%     yvect = [medyhf-roty medyhf+roty];
%     drawArrow(ax,xvect,yvect,xran,yran,'linewidth',1.5,'color','b');
%     [rotx, roty] = complex_rot(-5,0,-angbest1(i));
%     xvect = [medxhf-rotx medxhf+rotx];
%     yvect = [medyhf-roty medyhf+roty];
%     drawArrow(ax,xvect,yvect,xran,yran,'linewidth',1.5,'linestyle','--','color',[0.4 0.4 0.4]);
    xlim(ax, xran);
    ylim(ax, yran);
    xticks(ax,xran(1):5:xran(2));
    yticks(ax,yran(1):5:yran(2));    
    xlabel(ax,'E (km)','fontsize',11);
    ylabel(ax,'N (km)','fontsize',11);
%     text(ax,0.5,0.1,num2str(i),'FontSize',12,'unit','normalized');
    text(ax,0.83,0.95,strcat(num2str(size(mighf,1)),{' detections'}),'FontSize',8,...
         'unit','normalized','horizontalalignment','center');
    text(ax,0.83,0.90,strcat({'in '},num2str(trange(i,3)-trange(i,2)),{' s'}),'FontSize',8,...
         'unit','normalized','horizontalalignment','center');
    rate = sprintf('%.3f',size(mighf,1)/(trange(i,3)-trange(i,2)));
    text(ax,0.83,0.85,strcat({'rate: '},rate),'FontSize',8,...
         'unit','normalized','horizontalalignment','center');
%     text(ax,0.83,0.75,strcat(num2str(angbest1(i)),{' ^o'}),'FontSize',10,...
%          'unit','normalized','horizontalalignment','center');
%     text(ax,0.83,0.65,strcat(num2str(angbest2(i)),{' ^o'}),'FontSize',10,...
%          'unit','normalized','horizontalalignment','center','color','b'); 
    hold(ax,'off');


end

text(f.ax(1),0.04,0.93,'a','FontSize',11,'unit','normalized','EdgeColor','k','Margin',2);
text(f.ax(2),0.04,0.93,'b','FontSize',11,'unit','normalized','EdgeColor','k','Margin',2);
text(f.ax(3),0.04,0.93,'c','FontSize',11,'unit','normalized','EdgeColor','k','Margin',2);
text(f.ax(4),0.04,0.93,'d','FontSize',11,'unit','normalized','EdgeColor','k','Margin',2);

text(f.ax(1),0.5,0.1,'Original','FontSize',12,'unit','normalized','horizontalalignment','center');
text(f.ax(2),0.5,0.1,'Part 1','FontSize',12,'unit','normalized','horizontalalignment','center');
text(f.ax(3),0.5,0.1,'Part 2','FontSize',12,'unit','normalized','horizontalalignment','center');
text(f.ax(4),0.5,0.1,'Part 3','FontSize',12,'unit','normalized','horizontalalignment','center');

%%% save figure
print(f.fig,'-dpdf',strcat(rstpath,'/LZB.',num2str(nfam),'autortmv3.mig.divide.exp','.pdf'));


%% 3 migration needs to combined as 1
trange = rancom;

xran = [-20 25];
yran = [-20 20];

%%% define and position the figure frame and axes of each plot
f.fig=figure;
widin = 8;  % maximum width allowed is 8.5 inches
htin = 9;   % maximum height allowed is 11 inches
set(f.fig,'Position',[scrsz(1)+1*scrsz(3)/15 scrsz(2)+scrsz(4)/20 widin*res htin*res]);
nrow = 3;
ncol = 2;
for isub = 1:nrow*ncol-1
    f.ax(isub) = subplot(nrow,ncol,isub);
end

%%% reposition
set(f.ax(1), 'position', [ 0.08, 0.7, 0.315, 0.28]);
set(f.ax(2), 'position', [ 0.48, 0.7, 0.315, 0.28]);
set(f.ax(3), 'position', [ 0.08, 0.32, 0.315, 0.28]);
set(f.ax(4), 'position', [ 0.48, 0.32, 0.315, 0.28]); 
set(f.ax(5), 'position', [ 0.08, 0.05, 0.315, 0.2]); 
    
for i = 1: size(trange,1)-1
% for i = 1: 211
%     i=5;
%     disp(trange(i,:));
    indhf = find(hftime(:,13)==trange(i,1) & hftime(:,15)>=trange(i,2) & ...
                 hftime(:,15)<=trange(i,3));
    mighf = hftime(indhf,:);
        
    % subplot 1 of figure i
    ax = f.ax(i);
    hold(ax,'on');
    plot(ax,[-100 100],[0 0],'k--');
    plot(ax,[0 0],[-100 100],'k--');
    ax.FontSize = 9;
    mighf = sortrows(mighf,-15);
    scatter(ax,mighf(:,1),mighf(:,2), 20, mighf(:,15)/3600, 'filled','o');  %, 'MarkerEdgeColor', 'w')
    colormap(ax,'jet');
    c=colorbar(ax,'SouthOutside');
    pos = ax.Position;
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
    caxis(ax,[trange(i,2)/3600 trange(i,3)/3600])
    text(ax,0.85,0.1,'HF','FontSize',12,'unit','normalized');
%     text(ax,0.04,0.93,'a','FontSize',11,'unit','normalized','EdgeColor','k','Margin',2);
    text(ax,0.04,0.1,'TWKB','FontSize',11,'unit','normalized','EdgeColor','k','Margin',2);
    ax.Box = 'on';
    grid(ax, 'on');
    axis(ax, 'equal');
    ax.GridLineStyle = '--';
    ax.XAxisLocation = 'top';
%     medxhf = median(mighf(:,1));
%     medyhf = median(mighf(:,2));
%     [rotx, roty] = complex_rot(0,5,-angbest1(i));
%     xvect = [medxhf-rotx medxhf+rotx];
%     yvect = [medyhf-roty medyhf+roty];
%     drawArrow(ax,xvect,yvect,xran,yran,'linewidth',1.5);
%     [rotx, roty] = complex_rot(0,5,-angbest2(i));
%     xvect = [medxhf-rotx medxhf+rotx];
%     yvect = [medyhf-roty medyhf+roty];
%     drawArrow(ax,xvect,yvect,xran,yran,'linewidth',1.5,'color','b');
%     [rotx, roty] = complex_rot(-5,0,-angbest1(i));
%     xvect = [medxhf-rotx medxhf+rotx];
%     yvect = [medyhf-roty medyhf+roty];
%     drawArrow(ax,xvect,yvect,xran,yran,'linewidth',1.5,'linestyle','--','color',[0.4 0.4 0.4]);
    xlim(ax, xran);
    ylim(ax, yran);
    xticks(ax,xran(1):5:xran(2));
    yticks(ax,yran(1):5:yran(2));    
    xlabel(ax,'E (km)','fontsize',11);
    ylabel(ax,'N (km)','fontsize',11);
%     text(ax,0.5,0.1,num2str(i),'FontSize',12,'unit','normalized');
    text(ax,0.83,0.95,strcat(num2str(size(mighf,1)),{' detections'}),'FontSize',8,...
         'unit','normalized','horizontalalignment','center');
    text(ax,0.83,0.90,strcat({'in '},num2str(trange(i,3)-trange(i,2)),{' s'}),'FontSize',8,...
         'unit','normalized','horizontalalignment','center');
    rate = sprintf('%.3f',size(mighf,1)/(trange(i,3)-trange(i,2)));
    text(ax,0.83,0.85,strcat({'rate: '},rate),'FontSize',8,...
         'unit','normalized','horizontalalignment','center');
%     text(ax,0.83,0.75,strcat(num2str(angbest1(i)),{' ^o'}),'FontSize',10,...
%          'unit','normalized','horizontalalignment','center');
%     text(ax,0.83,0.65,strcat(num2str(angbest2(i)),{' ^o'}),'FontSize',10,...
%          'unit','normalized','horizontalalignment','center','color','b'); 
    hold(ax,'off');

end

i=size(trange,1);
angbest(i) = 120;

% create fit object with free constraints
fttpfree = fittype( @(a,b,x) a*x+b);

indhf = find(hftime(:,13)==trange(i,1) & hftime(:,15)>=trange(i,2) & ...
    hftime(:,15)<=trange(i,3));
mighf = hftime(indhf,:);

indlf = find(lftime(:,13)==trange(i,1) & lftime(:,15)>=trange(i,2) & ...
    lftime(:,15)<=trange(i,3));
miglf = lftime(indlf,:);
%%% Actually used is the pre-determined best prop direc to do the fitting
mighfdum = mighf;
for j = 1: size(mighf,1)
    x0 = mighfdum(j,1);
    y0 = mighfdum(j,2);
    [newx,newy] = coordinate_rot(x0,y0,-(angbest(i)-90),0,0);
    mighfdum(j,1) = newx;
    mighfdum(j,2) = newy;
end

miglfdum = miglf;
for j = 1: size(miglf,1)
    x0 = miglfdum(j,1);
    y0 = miglfdum(j,2);
    [newx,newy] = coordinate_rot(x0,y0,-(angbest(i)-90),0,0);
    miglfdum(j,1) = newx;
    miglfdum(j,2) = newy;
end

% detections involved into fitting
mintime = max(min(mighfdum(:,15)), min(miglfdum(:,15)));
maxtime = min(max(mighfdum(:,15)), max(miglfdum(:,15)));
mighfdum2 = mighfdum(mighfdum(:,15)>=mintime & mighfdum(:,15)<=maxtime, :);
miglfdum2 = miglfdum(miglfdum(:,15)>=mintime & miglfdum(:,15)<=maxtime, :);

% subplot 1 of figure i
ax = f.ax(i);
hold(ax,'on');
plot(ax,[-100 100],[0 0],'k--');
plot(ax,[0 0],[-100 100],'k--');
ax.FontSize = 9;
mighf = sortrows(mighf,-15);
scatter(ax,mighf(:,1),mighf(:,2), 20, mighf(:,15)/3600, 'filled','o');  %, 'MarkerEdgeColor', 'w')
colormap(ax,'jet');
c=colorbar(ax,'SouthOutside');
pos = ax.Position;
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
caxis(ax,[trange(i,2)/3600 trange(i,3)/3600])
text(ax,0.85,0.1,'HF','FontSize',12,'unit','normalized');
%     text(ax,0.04,0.93,'a','FontSize',11,'unit','normalized','EdgeColor','k','Margin',2);
text(ax,0.04,0.1,'TWKB','FontSize',11,'unit','normalized','EdgeColor','k','Margin',2);
ax.Box = 'on';
grid(ax, 'on');
axis(ax, 'equal');
ax.GridLineStyle = '--';
ax.XAxisLocation = 'top';
medxhf = median(mighf(:,1));
medyhf = median(mighf(:,2));
[rotx, roty] = complex_rot(0,5,-angbest(i));
xvect = [medxhf-rotx medxhf+rotx];
yvect = [medyhf-roty medyhf+roty];
drawArrow(ax,xvect,yvect,xran,yran,'linewidth',1.5);
[rotx, roty] = complex_rot(-5,0,-angbest(i));
xvect = [medxhf-rotx medxhf+rotx];
yvect = [medyhf-roty medyhf+roty];
drawArrow(ax,xvect,yvect,xran,yran,'linewidth',1.5,'linestyle','--','color',[0.4 0.4 0.4]);
xlim(ax, xran);
ylim(ax, yran);
xticks(ax,xran(1):5:xran(2));
yticks(ax,yran(1):5:yran(2));
xlabel(ax,'E (km)','fontsize',11);
ylabel(ax,'N (km)','fontsize',11);
%     text(ax,0.5,0.1,num2str(i),'FontSize',12,'unit','normalized');
text(ax,0.83,0.95,strcat(num2str(size(mighf,1)),{' detections'}),'FontSize',8,...
    'unit','normalized','horizontalalignment','center');
text(ax,0.83,0.90,strcat({'in '},num2str(trange(i,3)-trange(i,2)),{' s'}),'FontSize',8,...
    'unit','normalized','horizontalalignment','center');
rate = sprintf('%.3f',size(mighf,1)/(trange(i,3)-trange(i,2)));
text(ax,0.83,0.85,strcat({'rate: '},rate),'FontSize',8,...
    'unit','normalized','horizontalalignment','center');
text(ax,0.84,0.78,strcat(num2str(angbest(i)),{'{\circ}'}),'FontSize',10,...
         'unit','normalized','horizontalalignment','center');
hold(ax,'off');

% subplot 5 of figure ind(i)
ax = f.ax(5);
hold(ax,'on');
ax.FontSize = 9;
scatter(ax,mighfdum(:,15)/3600,mighfdum(:,1),20,[0.6 0.6 0.6],'filled','o',...
    'MarkerEdgeColor','k');

% create fit object

[fitobjhfprop,gofhfprop,outphfprop] = fit(mighfdum2(:,15)/3600, mighfdum2(:,1),fttpfree,'Robust',...
    'Bisquare','StartPoint',[1 1]);
% output fit parameters
coefprop = coeffvalues(fitobjhfprop);
slopeprophf(i) = coefprop(1);
intcptprophf(i) = coefprop(2);
fitprophf = feval(fitobjhfprop,mighfdum2(:,15)/3600);
plot(ax,mighfdum2(:,15)/3600,fitprophf,'-','linewidth',2,'color',[0.6 0.6 0.6]);

ax.Box = 'on';
grid(ax, 'on');
ax.GridLineStyle = '--';
xran1 = [trange(i,2)/3600 trange(i,3)/3600];
%     yran1 = [round(min([mighfdum(:,1);miglfdum(:,1)]))-1 ...
%              round(max([mighfdum(:,1);miglfdum(:,1)]))+1];
aa = round(prctile([mighfdum(:,1);miglfdum(:,1)], 98));
bb = round(prctile([mighfdum(:,1);miglfdum(:,1)], 2));
aap = aa + ceil((aa-bb)/6);
bbp = bb - ceil((aa-bb)/3);
yran1 = [bbp aap];
xlim(ax,xran1);
ylim(ax,yran1);
xlabel(ax,'Time (hr)','fontsize',11);
ylabel(ax,'Dist. along prop. (km)','fontsize',11);

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
slopesehf = gofhfprop.rmse./sqrt(sum((x-mean(x)).^2));
est = slopeprophf(i);
slopeprophfCI(i,:) = confidence_interval_general(est,slopesehf,length(x)-2,95);
interceptsehf = slopesehf.*sqrt(sum(x.^2)./length(x));
est = intcptprophf(i);
intcptprophfCI(i,:) = confidence_interval_general(est,interceptsehf,length(x)-2,95);
% sehf(i, :) = [gofhfprop.rmse slopesehf interceptsehf];


x = mighfdum2(:,15)/3600;
y = mighfdum2(:,1);
x_bar = wt_mean(x,wthf);
y_bar = wt_mean(y,wthf);
x_var = sum(wthf.*(x-x_bar).^2) / sum(wthf);
y_var = sum(wthf.*(y-y_bar).^2) / sum(wthf);
xy_cov = sum(wthf.*(x-x_bar).*(y-y_bar)) / sum(wthf);
pearwthf(i) = xy_cov / sqrt(x_var*y_var);

text(ax,0.4,0.2,sprintf('Slope: %.1f km/h',slopeprophf(i)),'FontSize',8,...
    'unit','normalized','horizontalalignment','left');
text(ax,0.97,0.2,sprintf('SE: %.2f',slopesehf),'FontSize',8,...
    'unit','normalized','horizontalalignment','right');
text(ax,0.4,0.13,sprintf('Pearson: %.3f',pearwthf(i)),'FontSize',8,...
    'unit','normalized','horizontalalignment','left');

hold(ax,'off');
    

text(f.ax(1),0.04,0.92,'a','FontSize',11,'unit','normalized','EdgeColor','k','Margin',2);
text(f.ax(2),0.04,0.92,'b','FontSize',11,'unit','normalized','EdgeColor','k','Margin',2);
text(f.ax(3),0.04,0.92,'c','FontSize',11,'unit','normalized','EdgeColor','k','Margin',2);
text(f.ax(4),0.04,0.92,'d','FontSize',11,'unit','normalized','EdgeColor','k','Margin',2);
text(f.ax(5),0.04,0.90,'e','FontSize',11,'unit','normalized','EdgeColor','k','Margin',2);

text(f.ax(1),0.5,0.1,'Original 1','FontSize',12,'unit','normalized','horizontalalignment','center');
text(f.ax(2),0.5,0.1,'Original 2','FontSize',12,'unit','normalized','horizontalalignment','center');
text(f.ax(3),0.5,0.1,'Original 3','FontSize',12,'unit','normalized','horizontalalignment','center');
text(f.ax(4),0.5,0.1,'Combined','FontSize',12,'unit','normalized','horizontalalignment','center');

%%% save figure
print(f.fig,'-dpdf',strcat(rstpath,'/LZB.',num2str(nfam),'autortmv3.mig.combine.exp','.pdf'));


%% 1 migration needs to be modified in time
trange = ranmod;
angbest = [230; 235];

xran = [-20 25];
yran = [-20 20];

%%% define and position the figure frame and axes of each plot
f.fig=figure;
widin = 8;  % maximum width allowed is 8.5 inches
htin = 9;   % maximum height allowed is 11 inches
set(f.fig,'Position',[scrsz(1)+1*scrsz(3)/15 scrsz(2)+scrsz(4)/20 widin*res htin*res]);
nrow = 2;
ncol = 2;
for isub = 1:nrow*ncol
    f.ax(isub) = subplot(nrow,ncol,isub);
end

%%% reposition
set(f.ax(1), 'position', [ 0.08, 0.64, 0.36, 0.32]);
set(f.ax(2), 'position', [ 0.52, 0.64, 0.36, 0.32]);
set(f.ax(3), 'position', [ 0.08, 0.35, 0.36, 0.22]);
set(f.ax(4), 'position', [ 0.52, 0.35, 0.36, 0.22]);

 % create fit object with free constraints
fttpfree = fittype( @(a,b,x) a*x+b);
   
for i = 1: size(trange,1)



indhf = find(hftime(:,13)==trange(i,1) & hftime(:,15)>=trange(i,2) & ...
    hftime(:,15)<=trange(i,3));
mighf = hftime(indhf,:);

indlf = find(lftime(:,13)==trange(i,1) & lftime(:,15)>=trange(i,2) & ...
    lftime(:,15)<=trange(i,3));
miglf = lftime(indlf,:);
%%% Actually used is the pre-determined best prop direc to do the fitting
mighfdum = mighf;
for j = 1: size(mighf,1)
    x0 = mighfdum(j,1);
    y0 = mighfdum(j,2);
    [newx,newy] = coordinate_rot(x0,y0,-(angbest(i)-90),0,0);
    mighfdum(j,1) = newx;
    mighfdum(j,2) = newy;
end

miglfdum = miglf;
for j = 1: size(miglf,1)
    x0 = miglfdum(j,1);
    y0 = miglfdum(j,2);
    [newx,newy] = coordinate_rot(x0,y0,-(angbest(i)-90),0,0);
    miglfdum(j,1) = newx;
    miglfdum(j,2) = newy;
end

% detections involved into fitting
mintime = max(min(mighfdum(:,15)), min(miglfdum(:,15)));
maxtime = min(max(mighfdum(:,15)), max(miglfdum(:,15)));
mighfdum2 = mighfdum(mighfdum(:,15)>=mintime & mighfdum(:,15)<=maxtime, :);
miglfdum2 = miglfdum(miglfdum(:,15)>=mintime & miglfdum(:,15)<=maxtime, :);

% subplot 1 of figure i
ax = f.ax(i);
hold(ax,'on');
plot(ax,[-100 100],[0 0],'k--');
plot(ax,[0 0],[-100 100],'k--');
ax.FontSize = 9;
mighf = sortrows(mighf,-15);
scatter(ax,mighf(:,1),mighf(:,2), 20, mighf(:,15)/3600, 'filled','o');  %, 'MarkerEdgeColor', 'w')
colormap(ax,'jet');
c=colorbar(ax,'SouthOutside');
pos = ax.Position;
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
caxis(ax,[trange(i,2)/3600 trange(i,3)/3600])
text(ax,0.85,0.1,'HF','FontSize',12,'unit','normalized');
%     text(ax,0.04,0.93,'a','FontSize',11,'unit','normalized','EdgeColor','k','Margin',2);
text(ax,0.04,0.1,'TWKB','FontSize',11,'unit','normalized','EdgeColor','k','Margin',2);
ax.Box = 'on';
grid(ax, 'on');
axis(ax, 'equal');
ax.GridLineStyle = '--';
ax.XAxisLocation = 'top';
medxhf = median(mighf(:,1));
medyhf = median(mighf(:,2));
[rotx, roty] = complex_rot(0,5,-angbest(i));
xvect = [medxhf-rotx medxhf+rotx];
yvect = [medyhf-roty medyhf+roty];
drawArrow(ax,xvect,yvect,xran,yran,'linewidth',1.5);
[rotx, roty] = complex_rot(-5,0,-angbest(i));
xvect = [medxhf-rotx medxhf+rotx];
yvect = [medyhf-roty medyhf+roty];
drawArrow(ax,xvect,yvect,xran,yran,'linewidth',1.5,'linestyle','--','color',[0.4 0.4 0.4]);
xlim(ax, xran);
ylim(ax, yran);
xticks(ax,xran(1):5:xran(2));
yticks(ax,yran(1):5:yran(2));
xlabel(ax,'E (km)','fontsize',11);
ylabel(ax,'N (km)','fontsize',11);
%     text(ax,0.5,0.1,num2str(i),'FontSize',12,'unit','normalized');
text(ax,0.83,0.95,strcat(num2str(size(mighf,1)),{' detections'}),'FontSize',8,...
    'unit','normalized','horizontalalignment','center');
text(ax,0.83,0.90,strcat({'in '},num2str(trange(i,3)-trange(i,2)),{' s'}),'FontSize',8,...
    'unit','normalized','horizontalalignment','center');
rate = sprintf('%.3f',size(mighf,1)/(trange(i,3)-trange(i,2)));
text(ax,0.83,0.85,strcat({'rate: '},rate),'FontSize',8,...
    'unit','normalized','horizontalalignment','center');
text(ax,0.84,0.78,strcat(num2str(angbest(i)),{'{\circ}'}),'FontSize',10,...
         'unit','normalized','horizontalalignment','center');
hold(ax,'off');

% subplot 5 of figure ind(i)
ax = f.ax(i+2);
hold(ax,'on');
ax.FontSize = 9;
scatter(ax,mighfdum(:,15)/3600,mighfdum(:,1),20,[0.6 0.6 0.6],'filled','o',...
    'MarkerEdgeColor','k');

% create fit object

[fitobjhfprop,gofhfprop,outphfprop] = fit(mighfdum2(:,15)/3600, mighfdum2(:,1),fttpfree,'Robust',...
    'Bisquare','StartPoint',[1 1]);
% output fit parameters
coefprop = coeffvalues(fitobjhfprop);
slopeprophf(i) = coefprop(1);
intcptprophf(i) = coefprop(2);
fitprophf = feval(fitobjhfprop,mighfdum2(:,15)/3600);
plot(ax,mighfdum2(:,15)/3600,fitprophf,'-','linewidth',2,'color',[0.6 0.6 0.6]);

ax.Box = 'on';
grid(ax, 'on');
ax.GridLineStyle = '--';
xran1 = [trange(i,2)/3600 trange(i,3)/3600];
%     yran1 = [round(min([mighfdum(:,1);miglfdum(:,1)]))-1 ...
%              round(max([mighfdum(:,1);miglfdum(:,1)]))+1];
aa = round(prctile([mighfdum(:,1);miglfdum(:,1)], 98));
bb = round(prctile([mighfdum(:,1);miglfdum(:,1)], 2));
aap = aa + ceil((aa-bb)/6);
bbp = bb - ceil((aa-bb)/3);
yran1 = [bbp aap];
xlim(ax,xran1);
ylim(ax,yran1);
xlabel(ax,'Time (hr)','fontsize',11);
ylabel(ax,'Dist. along prop. (km)','fontsize',11);

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
slopesehf = gofhfprop.rmse./sqrt(sum((x-mean(x)).^2));
est = slopeprophf(i);
slopeprophfCI(i,:) = confidence_interval_general(est,slopesehf,length(x)-2,95);
interceptsehf = slopesehf.*sqrt(sum(x.^2)./length(x));
est = intcptprophf(i);
intcptprophfCI(i,:) = confidence_interval_general(est,interceptsehf,length(x)-2,95);
% sehf(i, :) = [gofhfprop.rmse slopesehf interceptsehf];


x = mighfdum2(:,15)/3600;
y = mighfdum2(:,1);
x_bar = wt_mean(x,wthf);
y_bar = wt_mean(y,wthf);
x_var = sum(wthf.*(x-x_bar).^2) / sum(wthf);
y_var = sum(wthf.*(y-y_bar).^2) / sum(wthf);
xy_cov = sum(wthf.*(x-x_bar).*(y-y_bar)) / sum(wthf);
pearwthf(i) = xy_cov / sqrt(x_var*y_var);

text(ax,0.45,0.2,sprintf('Slope: %.1f km/h',slopeprophf(i)),'FontSize',8,...
    'unit','normalized','horizontalalignment','left');
text(ax,0.97,0.2,sprintf('SE: %.2f',slopesehf),'FontSize',8,...
    'unit','normalized','horizontalalignment','right');
text(ax,0.45,0.13,sprintf('Pearson: %.3f',pearwthf(i)),'FontSize',8,...
    'unit','normalized','horizontalalignment','left');

hold(ax,'off');

end


    

text(f.ax(1),0.04,0.92,'a','FontSize',11,'unit','normalized','EdgeColor','k','Margin',2);
text(f.ax(2),0.04,0.92,'b','FontSize',11,'unit','normalized','EdgeColor','k','Margin',2);
text(f.ax(3),0.04,0.91,'c','FontSize',11,'unit','normalized','EdgeColor','k','Margin',2);
text(f.ax(4),0.04,0.91,'d','FontSize',11,'unit','normalized','EdgeColor','k','Margin',2);

text(f.ax(1),0.5,0.1,'Original','FontSize',12,'unit','normalized','horizontalalignment','center');
text(f.ax(2),0.5,0.1,'Modified','FontSize',12,'unit','normalized','horizontalalignment','center');
% text(f.ax(3),0.5,0.1,'Original 3','FontSize',12,'unit','normalized','horizontalalignment','center');
% text(f.ax(4),0.5,0.1,'Combined','FontSize',12,'unit','normalized','horizontalalignment','center');

%%% save figure
print(f.fig,'-dpdf',strcat(rstpath,'/LZB.',num2str(nfam),'autortmv3.mig.modify.exp','.pdf'));










