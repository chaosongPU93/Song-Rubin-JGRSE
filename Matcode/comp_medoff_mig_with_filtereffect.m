% function comp_medoff_mig_with_filtereffect
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script is used to find out the difference in median offset between
% the HF and LF detections inside every migration that are tied to different
% fams, and compare this difference with the filtering effect from the same
% fams.
% This may tell us something about how good our filtering effects corrections
% based on the templates are? How well do they satisfy the tremor windows many
% of which never contributed to that template?
% Use the first part of parameters from 'mig_linear_fit_LZB_addfam_autortm_v3.m'
%
%
%
% Chao Song, chaosong@princeton.edu
% First created date:   2021/05/24
% Last modified date:   2021/05/24
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
            '017';
            '006';
            '001';
%             '158';      % 158, 20200916,testing purpose
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
%            -123.967000 48.458667 33.7800;       % 158, 20200916,testing purpose
           ];
       

relacont = loccont;
[dx, dy] = absloc2relaloc(loccont(:,1),loccont(:,2),loccont(2,1),loccont(2,2));
relacont(:,1) = dx;
relacont(:,2) = dy;

%%% load the location resolution points, +-1 & 2 samples
reslzb1 = load(fullfile(rstpath,'evtloc.13fam_locres_hf1spl'));
reslzb2 = load(fullfile(rstpath,'evtloc.13fam_locres_hf2spl'));

% sort according to day, sec
hftime = sortrows(hftime, [13, 15]);
lftime = sortrows(lftime, [13, 15]);

% this is new time ranges that come from 'identify_RTMs_v3.m'
trange = [
    2004195   4.1773e+04   4.3362e+04   % combine, y
    2004197   4.0974e+04   4.2205e+04   % speed direct 90, y, larger pear
    2004197   5.4453e+04   5.5867e+04   % speed direct 205, y, but use rmse 220
    2004197   5.6270e+04   5.7227e+04   % speed direct 175, y, but use rmse 185
%     2004197   7.9949e+04   8.1453e+04   % speed direct 260, y, but use rmse 195
    2004197   8.3402e+04   8.4888e+04   % divided into 3, 1 discarded, y
    2004197   8.5320e+04   8.6659e+04   % divided into 3, 1 discarded, y
    2004198   5.9530e+03   9.7940e+03   % speed direct 270, y, but use rmse 235
    2004198   1.9151e+04   2.1616e+04   % acceptted, y
    2004198   4.2632e+04   4.3795e+04   % acceptted, y
    2004198   5.3506e+04   5.8050e+04   % acceptted, y
    2004198   6.2966e+04   6.3885e+04   % acceptted, y
    2004198   7.3540e+04   7.6029e+04   % acceptted, y 
    2004198   8.4222e+04   8.6400e+04   % combine, y
    2004199   4.1670e+03   6.2660e+03   % speed direct 280, y, larger pear
    2004199   4.2845e+04   4.3359e+04   % speed direct 270, y, but use rmse 250
    2004199   4.6340e+04   4.9126e+04   % combine, care speed direct, y, larger pear
    2004199   4.9744e+04   5.0614e+04   % speed direct 180, y, but use rmse 195
    2004199   8.0861e+04   8.2008e+04   % divided into 2, y, care speed direct, y, but use rmse 90
    2004199   8.2008e+04   8.3515e+04   % divided into 2, y, care speed direct, RECHECK PEAR, slightly larger pear    
    2004200   1.2600e+04   1.5799e+04   % modified time, y
    2004200   1.5914e+04   1.7125e+04   % acceptted, y
    2004200   1.9104e+04   1.9900e+04   % check OLD param. to determine, y, use OLD end time
    2004200   4.8109e+04   4.8666e+04   % acceptted, y
    2004201   4.3700e+02   1.6290e+03   % speed direct 225, y, but use rmse 260
    2004203   1.6586e+04   2.0776e+04   % combine them, y
%     2005255   3.4400e+04   3.5612e+04   % speed direct 80, y, larger pear
    2005255   6.0381e+04   6.1557e+04   % speed direct 150, y, but use rmse 110
    2005255   6.7268e+04   6.8535e+04   % speed direct 80, y, larger pear
    2005255   8.4296e+04   8.5634e+04   % combine them, y
%     2005256   3.6200e+03   4.9570e+03   % speed direct 265, y, larger pear
    2005256   8.0330e+03   9.8280e+03   % change end time, y
    2005256   1.3412e+04   1.4346e+04   % speed direct 65, y, but use rmse 80   
%     2005256   1.9069e+04   2.0203e+04   % speed direct 195, y, larger pear
    2005256   2.1492e+04   2.2294e+04   % change start time, y
    2005256   2.8080e+04   2.9172e+04   % change start time, y
    2005256   3.4776e+04   3.5238e+04   % change start time, y
    2005256   3.6900e+04   3.8520e+04   % change times, y
%     2005256   4.6440e+04   4.8421e+04   % change start time, y, check speed direc, y, a bit awful
    2005256   5.2020e+04   5.2596e+04   % change times, y
    2005256   6.2640e+04   6.4865e+04   % change start time, y, check speed direc, y, but use rmse 130
    2005256   7.2720e+04   7.3914e+04   % divided into 2, y, but check speed direc, y
    2005256   7.6381e+04   7.7940e+04   % change end time, y
    2005257   6.1200e+03   7.7460e+03   % change start time, y, but check speed direc, y
%     2005257   1.2240e+04   1.3248e+04   % CHANGED time, st-3.68, y, use speed direc, larger pear, SUSPICOUS
    2005257   1.6920e+04   1.8655e+04   % change start time, y
    2005257   2.1070e+04   2.5365e+04   % accepted, y 
    2005257   3.9000e+04   4.5418e+04   % change start time, y
    2005257   6.1449e+04   6.4012e+04   % accepted, y
    2005257   7.3440e+04   7.8840e+04   % change end time, y, but check speed direc, y, but use rmse 235
    2005257   8.2159e+04   8.3380e+04   % speed direct 350, y, larger pear
    2005258   3.5000e+04   3.6540e+04   % divided into 3, y, but check speed direc, RECHECK PEAR, slightly larger pear, SUSPICOUS
    2005258   3.6540e+04   3.8160e+04   % divided into 3, y
    2005258   3.8160e+04   4.0021e+04   % divided into 3, y, but check speed direc y, slightly larger pear 
%     2005259   1.5400e+02   1.6340e+03   % speed direct 75, y, larger pear
%     2005259   1.9560e+03   4.0460e+03   % speed direct 185, y, but use rmse 140
    2005259   5.2320e+03   7.7140e+03   % check OLD param. to determine, y, but check speed direc, larger pear, use old direc 255
    2005259   3.7000e+04   4.1293e+04   % speed direct 85, RECHECK PEAR, y, but use rmse 115
%     2005260   2.6510e+03   3.4380e+03   % speed direct 80, y, but use rmse 115
    2005260   3.6040e+03   4.5000e+03   % divided into 2, y
%     2005260   4.6080e+03   6.6780e+03   % divided into 2, y, but check speed direc, y, larger pear
    2005260   5.6300e+04   5.8200e+04   % use old times, y
%     2005260   9.4090e+03   1.0357e+04   % speed direct 55, y, but use rmse 115
    ]; 
    
% hardcopy of angbest from either min. rmse or max. slope of accepted ones ntranacc
angbest = [
   250
    90
   220
   185
   225
   270
   235
   200
   230
    80
    65
   235
   245
   310
   250
   120
   195
    90
   140
   130
    70
   255
   250
   260
   120
   110
    80
   240
   220
    80
   115
   105
   260
    95
   260
   130
   230
   250
   195
   115
   250
   165
   245
   235
    10
    95
    85
   185
   255
   115
   250
   235
   ];

filthf = zeros(nfam, 2);
filtlf = zeros(nfam, 2);
filtcor = zeros(nfam, 2);
for ifam = 1: nfam
    %%% load the template filtering effect result 
    lolf = 0.5;
    hilf = 1.25;
    lohf = 1.25;
    hihf = 6.5;
    winsechf = 4;
    winseclf = 16;
    lofflf = 4;
    ccminlf = 0.35;
    fam = nfampool(ifam, :);
    PREFIX = strcat(fam,'.lf',num2str(lolf),'-',num2str(hilf),'.hf',num2str(lohf),'-',num2str(hihf),...
        '.hw',num2str(winsechf),'.lw',num2str(winseclf),'.lo',num2str(lofflf),'.ccm', ...
        num2str(ccminlf),'.','80sps');
    datapath = strcat(getenv('ALLAN'),'/data-no-resp','/LZBtrio');
    fname = strcat(datapath, '/MAPS/tempfeff_LZB_enc','_',PREFIX);
    filt = load(fname);
    spsratio = 20/80;
    filtlf(ifam, :) = filt(:,2)*spsratio;   % lf template shift due to filterig effect, sign is the same
    
    spsratio = 40/80;
    filthf(ifam, :) = filt(:,1)*spsratio;   % lf template shift due to filterig effect, sign is the same

    filtcor(ifam, :) = filt(:,3)*spsratio;   % filtering effect correction, in 40 sps
    
end

% CORFLAG = 'before';
CORFLAG = 'after';


%%
% percentage of detections from specific fam relative to all in each RTM
famperchf = zeros(length(trange), nfam);
famperclf = zeros(length(trange), nfam);

angle = 0:5:355;

slopehf = zeros(size(trange,1),length(angle));
rmsehf = zeros(size(trange,1),length(angle));
angrmsehf = zeros(size(trange,1),1);
angslopehf = zeros(size(trange,1),1);

% create fit object with free constraints
fttpfree = fittype( @(a,b,x) a*x+b);

fammig = cell(size(trange,1), 1);
diffmedmig = cell(size(trange,1), 1);
difftiedoffmig = cell(size(trange,1), 1);
for i = 1: size(trange,1)
% for i = 1: 211
%     i=46;
    disp(trange(i,:));
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
%     mighfdum2 = mighfdum;
%     miglfdum2 = miglfdum;

    %first you need to know which fams are involved in the particular migration 
    famhf = unique(mighfdum2(:,end));
    famlf = unique(miglfdum2(:,end));
    famtied = intersect(famlf, famhf);  % i need fams that have both LF and HF detections
    nfamtied = length(famtied);
%     famstr = [];
    diffmedtiedoff = [];
    difftiedoff = cell(nfamtied, 1);
    for j = 1: nfamtied
        %then identify the HF and LF detections tied to those fams
        hftied = mighfdum2(mighfdum2(:,end) == famtied(j), :);
        lftied = miglfdum2(miglfdum2(:,end) == famtied(j), :);
        hftiedoff = hftied(:, 9:10);
        lftiedoff = lftied(:, 9:10);

        %if we are questioning the median offset of detections before filtering effect, then we need
        %to add back the filtering effect to the corrected
        if famtied(j) < 10
            famstr(j, :) = strcat('00',num2str(famtied(j)));
        elseif famtied(j) < 100
            famstr(j, :) = strcat('0',num2str(famtied(j)));
        else
            famstr(j, :) = num2str(famtied(j));  
        end
        if strcmp(CORFLAG, 'before')    % using the offset before filtering effect correction            
            [~,ind] = ismember(famstr(j, :), nfampool, 'rows');
            hftiedoff = hftiedoff + filthf(ind,:);
            lftiedoff = lftiedoff + filtlf(ind,:);
        end
        
        %for these HF and LF detections, what are their median offset12 and offset13 in samples?
        %assume the HF locations are relaible enough
        medhftiedoff = median(hftiedoff, 1);
        
        %unify the sampling rate of LF, make it all 40 sps, same as HF
        lftiedoff = lftiedoff/20*40;
        
        %get the difference in offsets between HF median and each LF, and the corresponding median,
        %compare the difference in median offsets between HF and LF with the filtering effect of the
        %same fam
        %This makes sense as it recognizes that HF and LF have difference numbers, and it is
        %true for different migrations, we can use the number of LFs as a weight 
        tmp = lftiedoff - medhftiedoff;
        difftiedoff{j} = tmp;
        diffmedtiedoff(j,:) = median(tmp, 1);
%         filtcor(ind, :)

% %%%%%%%% NOT used anymore, as it didn't consider different number of detections   %%%%%%%%
%         medlftiedoff = median(lftiedoff, 1);
%         diffmedtiedoff2(j,:) = medlftiedoff - medhftiedoff;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
    fammig{i}  = famstr(1: nfamtied, :);
    diffmedmig{i} = diffmedtiedoff(1: nfamtied, :);
    difftiedoffmig{i} = difftiedoff;
end

%%
%re-structure the result, and make a cell array that is sorted by each migration
diffmedfam = cell(nfam, 1);
difftiedofffam = cell(nfam, 1);
for ifam = 1: nfam
    fam = nfampool(ifam, :);
    ii = 0;
    tmp3 = [];
    tmp5 = [];
    for i = 1: size(trange,1)
% i = 4;
        %for this migration, which fams were involved?
        tmp = char(fammig{i});
        if ismember(fam,tmp,'rows')     % is the current fam involved? 
            %if yes, what is the order?
            [~,ind] = ismember(fam,tmp,'rows');     
            %get the corresponding median lf-hf offset from all involved fams of this migration
            tmp2 = diffmedmig{i};
            %count
            ii = ii+1;
            %note the migration number
            tmp3(ii, 1) = i;
            %note the median lf-hf offset from this fam from this migration 
            tmp3(ii, 2:3) = tmp2(ind, :);
            
            %get the corresponding all lf-hf offsets (which were used to get the median above) from
            %all involved fams of this migration
            tmp4 = difftiedoffmig{i};
            tmp5{ii} = tmp4{ind};
        end
    end
    if ii == 0
        diffmedfam{ifam} = [];
        difftiedofffam{ifam} = [];
    else
        diffmedfam{ifam} = tmp3(1: ii, :);
        difftiedofffam{ifam} = tmp5;
    end
end

%%
%%% define and position the figure frame and axes of each plot
f.fig=figure;
f.fig.Renderer = 'painters';
widin = 8;  % maximum width allowed is 8.5 inches
htin = 9;   % maximum height allowed is 11 inches
set(f.fig,'Position',[scrsz(1)+1*scrsz(3)/10 scrsz(2)+scrsz(4)/20 widin*res htin*res]);
nrow = 3;
ncol = 3;
for isub = 1:nrow*ncol
    f.ax(isub) = subplot(nrow,ncol,isub);
    f.ax(isub).Box = 'on';
    grid(f.ax(isub), 'on');
    f.ax(isub).GridLineStyle = '--';
end

%%% reposition
set(f.ax(9), 'position', [ 0.7, 0.1, 0.25, 0.21]);
set(f.ax(8), 'position', [ 0.4, 0.1, 0.25, 0.21]);
set(f.ax(7), 'position', [ 0.1, 0.1, 0.25, 0.21]);
set(f.ax(6), 'position', [ 0.7, 0.4, 0.25, 0.21]);
set(f.ax(5), 'position', [ 0.4, 0.4, 0.25, 0.21]);
set(f.ax(4), 'position', [ 0.1, 0.4, 0.25, 0.21]);
set(f.ax(3), 'position', [ 0.7, 0.7, 0.25, 0.21]);
set(f.ax(2), 'position', [ 0.4, 0.7, 0.25, 0.21]);
set(f.ax(1), 'position', [ 0.1, 0.7, 0.25, 0.21]);

famcheck = ['043';
            '068';
            '147';
            '141';
            '099';
            '006';
            '125';
            '017';
            '144';];
    
for i = 1: size(famcheck,1)
    ax = f.ax(i);
    hold(ax, 'on');
    fam = famcheck(i, :);
    [~, ind] = ismember(fam,nfampool,'rows');
    tmp = diffmedfam{ind};
    if strcmp(CORFLAG, 'before')
        plot(ax,[-100 100],[filtcor(ind, 2) filtcor(ind, 2)],'--','color',[0.5 0.5 0.5],...
             'linew',1);
        plot(ax,[filtcor(ind, 1) filtcor(ind, 1)],[-100 100],'--','color',[0.5 0.5 0.5],...
             'linew',1);
        xm = max(abs(tmp(:,2)-filtcor(ind,1)))+1;
        ym = max(abs(tmp(:,3)-filtcor(ind,2)))+1;
        axis(ax, [filtcor(ind,1)-xm filtcor(ind,1)+xm filtcor(ind,2)-ym filtcor(ind,2)+ym]);
    elseif strcmp(CORFLAG, 'after')
        plot(ax,[-100 100],[0 0],'--','color',[0.5 0.5 0.5],'linew',1);
        plot(ax,[0 0],[-100 100],'--','color',[0.5 0.5 0.5],'linew',1);
        xm = max(abs(tmp(:,2)))+1;
        ym = max(abs(tmp(:,3)))+1;
        axis(ax, [-xm xm -ym ym]);
    end
    tmp2 = difftiedofffam{ind};
    nlf = zeros(size(tmp2, 2), 1);
    for j = 1: size(tmp2, 2)
        nlf(j) = size(tmp2{j}, 1);
    end
    %plot median offsets LF-HF from each involved migration
%     scatter(ax,tmp(:,2), tmp(:,3),20,[0.6 0.6 0.6],'filled');
    scatter(ax,tmp(:,2), tmp(:,3),20,nlf,'filled','markeredgecolor','k');
    oldc = colormap(ax,'gray');
    newc = flipud(oldc);
    colormap(ax,newc);
    c=colorbar(ax,'NorthOutside');
    pos = ax.Position;
    c.Position = [pos(1), pos(2)+0.214, pos(3), 0.015];
    if max(nlf) > 20
        tkint = 10;
    else
        tkint = 5;
    end
    c.Ticks = 0: tkint: max(nlf);
    c.TickLength = 0.03;
    c.Label.String = strcat('Number of LF detections');
    c.Label.FontSize = 8;
    caxis(ax,[0 max(nlf)]);
    %plot filtering effect
    scatter(ax,filtcor(ind, 1), filtcor(ind, 2),50,'k','filled','s','markeredgecolor','k');
    %plot the median of (median offsets from migration), considering diff. number of LF detections 
    %in diff. migrations
    comb = [];
    for j = 1: size(tmp2, 2)
        comb = [comb; tmp2{j}];
    end
    scatter(ax,median(comb(:,1),1), median(comb(:,2),1),120,'k','filled','p','markeredgecolor','k');
%     scatter(ax,median(tmp(:,2)), median(tmp(:,3)),120,'b','filled','p');
    text(ax,0.95,0.95,fam,'fontsize',10,'unit','normalized','Horizontalalignment','right');
%     text(ax,0.05,0.05,strcat(num2str(median(comb(:,1))-filtcor(ind, 1)),{', '},...
%          num2str(median(comb(:,2))-filtcor(ind, 2))),'fontsize',10,'unit','normalized',...
%          'Horizontalalignment','left');
    
    hold(ax, 'off');
end
xlabel(f.ax(7),'Med. diff. (12) (samples at 40 Hz)');
ylabel(f.ax(7),'Med. diff. (13) (samples at 40 Hz)');
% aa = supertit(f.ax(1:3), 'All possible migrations', 12);
% aa.FontWeight = 'bold';
print(f.fig,'-dpdf',strcat(rstpath,'/medoff_allmig_',CORFLAG,'cor_vs_filtereffect.pdf'));


%%
%%% define and position the figure frame and axes of each plot
f.fig=figure;
f.fig.Renderer = 'painters';
widin = 8;  % maximum width allowed is 8.5 inches
htin = 9;   % maximum height allowed is 11 inches
set(f.fig,'Position',[scrsz(1)+1*scrsz(3)/10 scrsz(2)+scrsz(4)/20 widin*res htin*res]);
nrow = 3;
ncol = 3;
for isub = 1:nrow*ncol
    f.ax(isub) = subplot(nrow,ncol,isub);
    f.ax(isub).Box = 'on';
    grid(f.ax(isub), 'on');
    f.ax(isub).GridLineStyle = '--';
end

%%% reposition
set(f.ax(9), 'position', [ 0.7, 0.1, 0.25, 0.21]);
set(f.ax(8), 'position', [ 0.4, 0.1, 0.25, 0.21]);
set(f.ax(7), 'position', [ 0.1, 0.1, 0.25, 0.21]);
set(f.ax(6), 'position', [ 0.7, 0.4, 0.25, 0.21]);
set(f.ax(5), 'position', [ 0.4, 0.4, 0.25, 0.21]);
set(f.ax(4), 'position', [ 0.1, 0.4, 0.25, 0.21]);
set(f.ax(3), 'position', [ 0.7, 0.7, 0.25, 0.21]);
set(f.ax(2), 'position', [ 0.4, 0.7, 0.25, 0.21]);
set(f.ax(1), 'position', [ 0.1, 0.7, 0.25, 0.21]);

famcheck = ['043';
            '068';
            '147';
            '141';
            '099';
            '006';
            '125';
            '017';
            '144';];
        
indne = [10,18,21,30,34,46,47];
indsw = [5,6,7,9,12,13,15,22,23,28,33,35,37,38,41,43,44,49,51,52];
indmig = reshape(union(indsw, indne), [], 1);
for i = 1: size(famcheck,1)
    ax = f.ax(i);
    hold(ax, 'on');
    fam = famcheck(i, :);
    [~, ind] = ismember(fam,nfampool,'rows');
    tmp = diffmedfam{ind};
    [~, ind2] = intersect(tmp(:,1), indmig);    
    if strcmp(CORFLAG, 'before')
        plot(ax,[-100 100],[filtcor(ind, 2) filtcor(ind, 2)],'k--');
        plot(ax,[filtcor(ind, 1) filtcor(ind, 1)],[-100 100],'k--');
%         xm = max(abs(tmp(ind2,2)-filtcor(ind,1)))+1;
%         ym = max(abs(tmp(ind2,3)-filtcor(ind,2)))+1;
        xm = max(abs(tmp(:,2)-filtcor(ind,1)))+1;
        ym = max(abs(tmp(:,3)-filtcor(ind,2)))+1;
        axis(ax, [filtcor(ind,1)-xm filtcor(ind,1)+xm filtcor(ind,2)-ym filtcor(ind,2)+ym]);
    elseif strcmp(CORFLAG, 'after')
        plot(ax,[-100 100],[0 0],'k--');
        plot(ax,[0 0],[-100 100],'k--');
%         xm = max(abs(tmp(ind2,2)))+1;
%         ym = max(abs(tmp(ind2,3)))+1;
        xm = max(abs(tmp(:,2)))+1;
        ym = max(abs(tmp(:,3)))+1;
        axis(ax, [-xm xm -ym ym]);
    end
    tmp2 = difftiedofffam{ind};
    nlf = zeros(size(tmp2, 2), 1);
    for j = 1: size(tmp2, 2)
        nlf(j) = size(tmp2{j}, 1);
    end
    %plot median offsets LF-HF from each involved SW migration
    [~, ind3] = intersect(tmp(:,1), indsw);
    scatter(ax,tmp(ind3,2), tmp(ind3,3),20,nlf(ind3),'filled','markeredgecolor','k');
    colormap(ax,'whitered');
    c=colorbar(ax,'NorthOutside');
    pos = ax.Position;
    c.Position = [pos(1), pos(2)+0.214, pos(3)/2-0.003, 0.015];
%     c.FontSize = 7;
    if max(nlf(ind3)) > 20
        tkint = 10;
    else
        tkint = 5;
    end
    c.Ticks = 0: tkint: max(nlf(ind3));
    c.TickLength = 0.05;
    c.Label.String = strcat('Num. of LF (WSW)');
    c.Label.FontSize = 8;
    caxis(ax,[0 max(nlf(ind3))]);
    
    %plot median offsets LF-HF from each involved NE migration
    [~, ind4] = intersect(tmp(:,1), indne);
    ax2 = axes;
    scatter(ax2,tmp(ind4,2), tmp(ind4,3),20,nlf(ind4),'filled','markeredgecolor','k');
    hold(ax2,'on');
    %link two axes
    linkaxes([ax,ax2])
    %Hide the top axes
    ax2.Visible = 'off';
    ax2.XTick = [];
    ax2.YTick = [];
    ax2.Position = ax.Position;
    colormap(ax2,'whiteblue');
    c2=colorbar(ax2,'NorthOutside');
    c2.Position = [pos(1)+pos(3)/2+0.003, pos(2)+0.214, pos(3)/2-0.003, 0.015];
    if max(nlf(ind4)) > 20
        tkint = 10;
    else
        tkint = 5;
    end
    c2.Ticks = 5: tkint: max(nlf(ind4));
    c2.TickLength = 0.05;
    c2.Label.String = strcat('Num. of LF (ENE)');
    c2.Label.FontSize = 8;
    caxis(ax2,[0 max(nlf(ind4))]);

    %plot the median of (median offsets from migration), considering diff. number of LF detections 
    %in diff. migrations
    combsw = [];
    for j = 1: size(ind3, 1)
        combsw = [combsw; tmp2{ind3(j)}];
    end
    scatter(ax2,median(combsw(:,1),1), median(combsw(:,2),1),100,[1 96/255 0],'filled','p',...
            'markeredgecolor','k');
    combne = [];
    for j = 1: size(ind4, 1)
        combne = [combne; tmp2{ind4(j)}];
    end
    scatter(ax2,median(combne(:,1),1), median(combne(:,2),1),100,[0 1 1],'filled','p',...
            'markeredgecolor','k');
    combuse = [];
    for j = 1: size(ind2, 1)
        combuse = [combuse; tmp2{ind2(j)}];
    end
    scatter(ax2,median(combuse(:,1),1), median(combuse(:,2),1),120,'k','filled','p',...
            'markeredgecolor','k');
%     scatter(ax,median(tmp(ind3,2)), median(tmp(ind3,3)),80,'r','filled','p','markeredgecolor','k');
%     scatter(ax,median(tmp(ind4,2)), median(tmp(ind4,3)),80,'b','filled','p','markeredgecolor','k');
%     scatter(ax,median(tmp(ind2,2)), median(tmp(ind2,3)),120,'k','filled','p','markeredgecolor','k');
    scatter(ax2,filtcor(ind, 1), filtcor(ind, 2),50,'k','filled','s','markeredgecolor','k');
    text(ax,0.95,0.95,fam,'fontsize',10,'unit','normalized','Horizontalalignment','right');
%     text(ax,0.05,0.05,strcat(num2str(median(combuse(:,1))),{', '},...
%          num2str(median(combuse(:,2)))),'fontsize',10,'unit','normalized',...
%          'Horizontalalignment','left');
    hold(ax, 'off');
    hold(ax2,'off');
end
xlabel(f.ax(7),'Med. diff. (12) (samples at 40 Hz)');
ylabel(f.ax(7),'Med. diff. (13) (samples at 40 Hz)');
% aa = supertit(f.ax(1:3), 'Migrations with opposing directions', 12);
% aa.FontWeight = 'bold';
print(f.fig,'-dpdf',strcat(rstpath,'/medoff_usedmig_',CORFLAG,'cor_vs_filtereffect.pdf'));







