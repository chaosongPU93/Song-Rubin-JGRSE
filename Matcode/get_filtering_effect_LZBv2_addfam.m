% function get_filtering_effect_LZBv2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This is the function to get the filtering effect under different frequency
% band for LZB trio. Use the Bostock's catalog in one file
% '/BOSTOCK/total_mag_detect_0000_cull_NEW.txt'. Use the similar way to do
% the cross correlatation in Allan's detection code with all 3 stations,
% under the constraint that the enclosed sum of the offset need to be zero
% so that it necessarily would reach the individual max cc.
%
%   stack first, then filter with different band. It turns out that order of
%   filtering/stack doesn't really affect the result too much under 40 or
%   80 sps base
%
%
% By Chao Song, chaosong@princeton.edu
% First created date:   2020/08/26
% Last modified date:   2020/08/26
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Initialization
format short e   % Set the format to 5-digit floating point
clear
clc
% close all

% set(0,'DefaultFigureVisible','on');
set(0,'DefaultFigureVisible','off');   % switch to show the plots or not

% [scrsz, res] = pixelperinch(1);
scrsz=get(0,'ScreenSize');   % for example, scrsz=[1,1,1680,1050]
wid=scrsz(3);       % width
hite=scrsz(4);       % height

workpath = getenv('ALLAN');
%%% choose the data path here
% datapath = workpath;
datapath = strcat(workpath,'/data-no-resp');
temppath = strcat(datapath, '/templates/LZBtrio/');
rstpath = strcat(datapath, '/LZBtrio');

% family pool
ofampool = ['002';
    '013';
    '025';
    '028';
    '056';
    '084';
    '099';
    '115';
    '125';
    '147';
    '149'];

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

% according to raw detections, 001 a new fam worthwhile to add
nfampool = ['001'];

% from the obtained templates, 231 is apparently not a good family
% 234 is bit weird on timing, but the shape of dipole seems good enough
% overall, 234 and 158 is comparable, 015 is less good,
% given the their locations, let's try 234 and 158 together!
nfampool = [
           '001';
           '006';
%            '015';
           '158';
%            '231';
           '234'];

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

stas=['TWKB '
    'LZB  '
    'MGCB '];     % determine the trio and order

% number of used stations
nsta=size(stas,1);         %  number of stations

% the summary array for filtering effect
filtall = [];   % filtering effect in hf&lf for all fams
decall = [];    % cc width estimate in hf&lf for all fams

%% make template or read template
remake = 0;  % re-make the template or not, 1/0
dstack = [];
ccstack = [];
sps = 80;
templensec = 60;
for ifam = 1: size(nfampool,1)
%     ifam=2;
    fam = nfampool(ifam, :)
    
    if remake
        ccmethod = 2; % only using lfes passed the cc thresholds (ccmethod=1); using all lfes (ccmethod=2)  
        plflag = 0;
        [dstack, ccstack] = mk_bbtemp_LZB(fam,sps,ccmethod,plflag);
        % write into files
        for ista = 1: nsta
            fid = fopen(strcat(temppath, fam, '_', strtrim(stas(ista, :)), '_', ...
                num2str(sps), 'sps_', num2str(templensec), 's_', ...
                'BBDS_', 'opt_Nof_Non_Chao'), 'w+');  % bandpassed direct stack, no filter, no norm
            
            fprintf(fid, '%f \n', dstack(ista, :)');
            fclose(fid);
            
            fid = fopen(strcat(temppath, fam, '_', strtrim(stas(ista, :)), '_', ...
                num2str(sps), 'sps_', num2str(templensec), 's_', ...
                'BBCCS_', 'opt_Nof_Non_Chao'), 'w+');  % bandpassed cc stack, no filter, no norm
            fprintf(fid, '%f \n', ccstack(ista, :)');
            fclose(fid);
        end
        ind = [5 4 6];
        %     stack = dstack(ind, :);
        stack = ccstack(ind, :);
        
    else
        for ista = 1: nsta
            fname = strcat(temppath, fam, '_', strtrim(stas(ista, :)), '_', num2str(sps), 'sps_', ...
                num2str(templensec), 's_', 'BBDS_', 'opt_Nof_Non_Chao');
            dstack(ista,:) = load(fname);
            fname = strcat(temppath, fam, '_', strtrim(stas(ista, :)), '_', num2str(sps), 'sps_', ...
                num2str(templensec), 's_', 'BBCCS_', 'opt_Nof_Non_Chao');
            ccstack(ista,:) = load(fname);
        end
        stack = ccstack;
    end

    
    %% normalization
    templen = size(stack,2);
    stackno = zeros(nsta, templen);
    mid = templen/2;

    for ista = 1: nsta

        ampmax = max(stack(ista, mid-sps/2+1: mid+sps/2));
        ampmin = min(stack(ista, mid-sps/2+1: mid+sps/2));
        if ampmax >= -ampmin
            fprintf('At station %s positive is bigger \n', strtrim(stas(ista, :)));
        else
            fprintf('At station %s negative is bigger \n', strtrim(stas(ista, :)));
        end
        norm = max(ampmax, -ampmin);
        stackno(ista, :) = stack(ista, :)/ norm;
    end


    %% plot
    %%% plot the 1-step stacked template
    %%% figure 1
    figure('Position',[0.1*wid 0.1*hite 0.3*wid 0.3*hite]);
    for ista=1: nsta
        plot(stack(ista, :) + 0.1* ista, 'linewidth', 2); hold on
        title('direct stack');
    end
    box on
    grid on

    %%% figure 2
    figure('Position',[0.1*wid 0.3*hite 0.3*wid 0.3*hite]);
    for ista = 1: nsta
        plot(stack(ista, mid-5*sps: mid+5*sps) + 0.1* ista, 'linewidth', 2); hold on
        title('direct stack-short window');
    end
    box on
    grid on

    %%% plot the normalized 1-step template
    %%% figure 3
    figure('Position',[0.1*wid 0.2*hite 0.3*wid 0.3*hite]);
    color=['r','b','k'];
    for ista = 1: nsta
        plot((1:10*sps+1)/sps,stackno(ista, mid-5*sps: mid+5*sps) + 1* ista, 'linewidth', 1,'color',...
            color(ista)); hold on
        %     plot(stackno(ista, mid-5*sps: mid+5*sps) + 1* ista, 'linewidth', 1,'color',...
        %          color(ista)); hold on
        text(1,1*ista,stas(ista,:),'fontsize',12);
    end
    title('normalized direct stack');
    xlabel('Time (s)');
    % xlabel('Samples (on 40 sps)');
    ylabel('Normalized Amp.');
    box on
    grid on
    
    
    %% HF
    wlensechf = 4;
    wlenhf = wlensechf*sps;
    lohf = 1.25;
    hihf = 6.5;
    npo = 2;
    npa = 2;
    mshift = 14;    % 29/14, samples based on sps in detection, not on 80 sps, but should not affect
    loopoffmax = 1.5;
    xcmaxAVEnmin = 0.4;
    
    for ista = 1: nsta
        stacktemp(ista,:) = Bandpass(stack(ista,:), sps, lohf, hihf, npo, npa, 'butter');
    end
    stackuse = stacktemp(1:3,:);
    stackauto = stackuse.*stackuse;
    lenx = templen-2*mshift;
    stack12x = zeros(lenx, 2*mshift+1);
    stack13x = zeros(lenx, 2*mshift+1);    % stas: 1->PGC, 2->SSIB, 3->SILB
    stack32x = zeros(lenx, 2*mshift+1);
    
    for n=-mshift:mshift
        % PGC corr SSIB, 1+mshift:tracelen-mshift == 1+mshift-n:tracelen-mshift-n == lenx
        stack12x(:,n+mshift+1)=stackuse(1,1+mshift:templen-mshift).* ...
            stackuse(2,1+mshift-n:templen-mshift-n);
        % PGC corr SILB
        stack13x(:,n+mshift+1)=stackuse(1,1+mshift:templen-mshift).* ...
            stackuse(3,1+mshift-n:templen-mshift-n);
        % SILB corr SSIB
        stack32x(:,n+mshift+1)=stackuse(3,1+mshift:templen-mshift).* ...
            stackuse(2,1+mshift-n:templen-mshift-n);
    end
    
    sumstack12=zeros(1,2*mshift+1);   % PGC-SSIB
    sumstack13=zeros(1,2*mshift+1);   % PGC-SILB
    sumstack32=zeros(1,2*mshift+1);   % SILB-SSIB
    sumstack1sq=zeros(1,2*mshift+1);  % "sq" are the running cumulative sum, to be used later by differencing the window edpoints
    sumstack2sq=zeros(1,2*mshift+1);
    sumstack3sq=zeros(1,2*mshift+1);
    sumstack3Bsq=zeros(1,2*mshift+1);
    
    istart = mid-wlenhf/2;
    iend = istart+wlenhf-1;
    sumstack12(1,:) = sum(stack12x(istart-mshift: iend-mshift, :));
    sumstack13(1,:) = sum(stack13x(istart-mshift: iend-mshift, :));
    sumstack32(1,:) = sum(stack32x(istart-mshift: iend-mshift, :));
    sumstack1sq(1,:) = sum(stackauto(1, istart: iend));
    sumstack3Bsq(1,:) = sum(stackauto(3, istart: iend));
    for m = -mshift:mshift
        sumstack2sq(1,m+mshift+1)=sum(stackauto(2, istart-m: iend-m)); %+m??? (yes).
        sumstack3sq(1,m+mshift+1)=sum(stackauto(3, istart-m: iend-m));
    end
    
    %An attempt to bypass glitches in data.  Min value of good data typically ~10^{-2}
    glitches=1.e-7;
    sumstack1sq=max(sumstack1sq,glitches);
    sumstack2sq=max(sumstack2sq,glitches);
    sumstack3sq=max(sumstack3sq,glitches);    % return maximum between A and B
    
    denomstack12n=realsqrt(sumstack1sq.*sumstack2sq);    % Real square root, An error is produced if X is negative
    denomstack13n=realsqrt(sumstack1sq.*sumstack3sq);
    denomstack32n=realsqrt(sumstack3Bsq.*sumstack2sq);
    
    sumstack12n=sumstack12./denomstack12n;   % suffix 'n' means normalized
    sumstack13n=sumstack13./denomstack13n;
    sumstack32n=sumstack32./denomstack32n;
    
    [xcmaxstack12n,imaxstack12]=max(sumstack12n,[],2);   %Integer-offset max cross-correlation
    [xcmaxstack13n,imaxstack13]=max(sumstack13n,[],2);   % along row, max cc val and index in each window
    [xcmaxstack32n,imaxstack32]=max(sumstack32n,[],2);
    
    %Parabolic fit:
    [xmaxstack12n,ymaxstack12n,aSTA12]=parabol(1,mshift,sumstack12n,imaxstack12); %Interpolated max cross-correlation
    [xmaxstack13n,ymaxstack13n,aSTA13]=parabol(1,mshift,sumstack13n,imaxstack13);
    [xmaxstack32n,ymaxstack32n,aSTA32]=parabol(1,mshift,sumstack32n,imaxstack32);
    
    %Center them
    imaxstack12cent=imaxstack12-mshift-1;  % "cent" is "centered"; imaxSTA12 is original 1: 2*mshift+1, corresponds to -mshift: mshift
    imaxstack13cent=imaxstack13-mshift-1;
    imaxstack32cent=imaxstack32-mshift-1;
    %%% NOTICE: the right order of a closed 3-sta pair is +13, -12, +32, where 13 means 1-->3
    iloopoff=imaxstack13cent-imaxstack12cent+imaxstack32cent; %How well does the integer loop close?
    %
    xmaxstack12n=xmaxstack12n-mshift-1;
    xmaxstack13n=xmaxstack13n-mshift-1;
    xmaxstack32n=xmaxstack32n-mshift-1;
    loopoff=xmaxstack13n-xmaxstack12n+xmaxstack32n; %How well does the interpolated loop close?
    loopoffhf(ifam) = loopoff;
    xcmaxAVEn=(xcmaxstack12n+xcmaxstack13n+xcmaxstack32n)/3;
    medxcmaxAVEn=median(xcmaxAVEn);
    xmaxstack12ntmphf=xmaxstack12n;    % tmp == temporary
    xmaxstack13ntmphf=xmaxstack13n;
    xmaxstack32ntmphf=xmaxstack32n;
    
    if xcmaxAVEn<xcmaxAVEnmin || abs(loopoff)>loopoffmax || isequal(abs(imaxstack12cent),mshift)...
       || isequal(abs(imaxstack13cent),mshift) || isequal(abs(imaxstack32cent),mshift)
        xmaxstack12ntmphf=mshift+1; xmaxstack13ntmphf=mshift+1; xmaxstack32ntmphf=mshift+1; %dummy them, if these criteria are met
        disp('WRONG! The basic criteria is not met');
    else
        xcmaxconprev=-99999.;  %used to be 0; not good with glitches
        %%% NOTE here are using offset imaxSTA12/13/32 without centering, i.e. [1, 2shift+1]
        % Usually, the loop is happened between imaxSTA12n +- floor(loopoffmax+1)
        % width of which is 2*floor(loopoffmax+1)
        %%% floor(2.5)=2; floor(-2.6)=-3
        for iSTA12 = max(1,imaxstack12-floor(loopoffmax+1)): ...
                     min(imaxstack12+floor(loopoffmax+1),2*mshift+1)
            for iSTA13 = max(1,imaxstack13-floor(loopoffmax+1)): ...
                         min(imaxstack13+floor(loopoffmax+1),2*mshift+1)
                ibangon = (mshift+1)-iSTA13+iSTA12;     %%% SEE NOTES #2019/03/17# page 66 to understand better
                %%% i.e., -mshift <= -iSTA13+iSTA12 <= mshift
                if ibangon >= 1 && ibangon <= 2*mshift+1
                    xcmaxcon=sumstack12n(iSTA12)+sumstack13n(iSTA13)+sumstack32n(ibangon);
                    if xcmaxcon > xcmaxconprev
                        xcmaxconprev=xcmaxcon;
                        iSTA12bang=iSTA12;
                        iSTA13bang=iSTA13;
                    end
                end
            end
        end
        %%% will result in the max xcmaxcon and corresponding iSTA12,
        %%% iSTA13, and save them into xcmaxconprev, iSTA12bang and iSTA13bang
        iSTA32bang=(mshift+1)-iSTA13bang+iSTA12bang;
        if abs(iSTA12bang-imaxstack12) <= loopoffmax && ...  %not sure if these 3 lines are satisfied automatically ...
           abs(iSTA13bang-imaxstack13) <= loopoffmax && ...  % SHOULD not be, i think, could be floor(loopoffmax+1) > loopoffmax
           abs(iSTA32bang-imaxstack32) <= loopoffmax && ...
           sumstack12n(iSTA12bang)+sumstack13n(iSTA13bang)+sumstack32n(iSTA32bang) ...
               >= 3*xcmaxAVEnmin   % xcmaxAVEnmin, predetermined
            %%% ALSO, sumsSTA12n(n,iSTA12bang) == sumsSTA12nn(iSTA12bang)
            
            xmaxstack12ntmphf=iSTA12bang-(mshift+1); %without interpolation this is just centering.
            xmaxstack13ntmphf=iSTA13bang-(mshift+1);
            xmaxstack32ntmphf=iSTA32bang-(mshift+1);
            
            %%% xcmaxAVEnbang is added by Chao, to distinguish from
            %%% xcmaxAVEn, because it is max average CC coef
            xcmaxAVEnbang=(sumstack12n(iSTA12bang)+sumstack13n(iSTA13bang)+sumstack32n(iSTA32bang))/3;
        end
    end
    
%     xmaxstack12ntmphf
%     xmaxstack13ntmphf
%     xmaxstack32ntmphf

    if xmaxstack13ntmphf-xmaxstack12ntmphf+xmaxstack32ntmphf ~=0
        disp('WRONG! Loopoff is not enclosed');
    else
        disp('Loopoff in hf is 0');
    end
    filthf = [xmaxstack12ntmphf;  xmaxstack13ntmphf];
    
    % examine the widness of cc function
    cc = [sumstack12n; sumstack13n; sumstack32n];    
    off = [imaxstack12; imaxstack13; imaxstack32];
    
    for i = 1: 3 
        dechf(1,i) = (2*cc(i,off(i))-cc(i,off(i)-1)-cc(i,off(i)+1))/2; % decline of cc for +-1 sample off
        dechf(2,i) = (2*cc(i,off(i))-cc(i,off(i)-3)-cc(i,off(i)+3))/2; % 3 samples
        dechf(3,i) = (2*cc(i,off(i))-cc(i,off(i)-6)-cc(i,off(i)+6))/2; % 6 samples
        dechf(4,i) = (2*cc(i,off(i))-cc(i,off(i)-10)-cc(i,off(i)+10))/2;  % 10 samples
    end
    avedechf = sum(dechf,2)/3;
    dechf(:,4) = avedechf;

    
    % sumstackxxn is actually the cross-correlation between stations allowing mshift
    figure
    hold on
    p1=plot(-mshift:mshift,sumstack12n,'b-','linewidth', 1);
    p2=plot(-mshift:mshift,sumstack13n,'r-','linewidth', 1);
    p3=plot(-mshift:mshift,sumstack32n,'k-','linewidth', 1);
    
    p4=plot(imaxstack12cent,sumstack12n(imaxstack12),'bs','markersize',6);
    p5=plot(imaxstack13cent,sumstack13n(imaxstack13),'ro','markersize',6);
    p6=plot(imaxstack32cent,sumstack32n(imaxstack32),'kd','markersize',6);
    % plot(xmaxstack12n,ymaxstack12n,'bs','markersize',10);
    % plot(xmaxstack13n,ymaxstack13n,'ro','markersize',10);
    % plot(xmaxstack32n,ymaxstack32n,'kd','markersize',10);
    p7=plot(xmaxstack12ntmphf,sumstack12n(xmaxstack12ntmphf+mshift+1),'bs','markersize',10);
    p8=plot(xmaxstack13ntmphf,sumstack13n(xmaxstack13ntmphf+mshift+1),'ro','markersize',10);
    p9=plot(xmaxstack32ntmphf,sumstack32n(xmaxstack32ntmphf+mshift+1),'kd','markersize',10);
    text(-5,-0.4,'1.25-6.5 Hz','fontsize',12);
    legend([p1,p2,p3,p4,p5,p6,p7,p8,p9],{'CC TWKB-LZB(12)','CC TWKB-KLNB(13)','CC KLNB-LZB(32)',...
        'individual \Delta{t}_{12}',...
        'individual \Delta{t}_{13}', 'individual \Delta{t}_{32}',...
        'constrained \Delta{t}_{12}','constrained \Delta{t}_{13}',...
        'constrained \Delta{t}_{32}'});
    xlabel(strcat('Samples (on',{' '},num2str(sps),{' '}, 'sps)'),'fontsize',12);
    ylabel('Normalized CC Coef.','fontsize',12);
    box on;
    grid on;
    % vertical_cursors
    
%     % save fig
%     fignm = strcat(rstpath, '/FIGS/tempfeff_LZB_enc_',fam,'.hf',num2str(lohf),'-',num2str(hihf),...
%         '.hw',num2str(wlensechf),'.lo',num2str(loopoffmax),'.ccm', ...
%         num2str(xcmaxAVEnmin),'.',num2str(sps),'sps.eps');
%     print('-depsc',fignm);
 
    %%%% end for hf
    
    
    %% LF
    wlenseclf = 16;
    wlenlf = wlenseclf*sps;
    lolf = 0.5;
    hilf = 1.25;
    npo = 2;
    npa = 2;
    cyclskip = 20;
    mshift=14+cyclskip;  % 19/14, samples based on sps in detection, not on 80 sps, but should not affect
    loopoffmax = 4;
    xcmaxAVEnmin = 0.35;  % 0.4/0.35
    
    for ista = 1: nsta
        stacktemp(ista,:) = Bandpass(stack(ista,:), sps, lolf, hilf, npo, npa, 'butter');
    end
    stackuse = stacktemp(1:3,:);
    stackauto = stackuse.*stackuse;
    lenx = templen-2*mshift;
    stack12x = zeros(lenx, 2*mshift+1);
    stack13x = zeros(lenx, 2*mshift+1);    % stas: 1->PGC, 2->SSIB, 3->SILB
    stack32x = zeros(lenx, 2*mshift+1);
    
    for n=-mshift:mshift
        % PGC corr SSIB, 1+mshift:tracelen-mshift == 1+mshift-n:tracelen-mshift-n == lenx
        stack12x(:,n+mshift+1)=stackuse(1,1+mshift:templen-mshift).* ...
            stackuse(2,1+mshift-n:templen-mshift-n);
        % PGC corr SILB
        stack13x(:,n+mshift+1)=stackuse(1,1+mshift:templen-mshift).* ...
            stackuse(3,1+mshift-n:templen-mshift-n);
        % SILB corr SSIB
        stack32x(:,n+mshift+1)=stackuse(3,1+mshift:templen-mshift).* ...
            stackuse(2,1+mshift-n:templen-mshift-n);
    end
    
    sumstack12=zeros(1,2*mshift+1);   % PGC-SSIB
    sumstack13=zeros(1,2*mshift+1);   % PGC-SILB
    sumstack32=zeros(1,2*mshift+1);   % SILB-SSIB
    sumstack1sq=zeros(1,2*mshift+1);  % "sq" are the running cumulative sum, to be used later by differencing the window edpoints
    sumstack2sq=zeros(1,2*mshift+1);
    sumstack3sq=zeros(1,2*mshift+1);
    sumstack3Bsq=zeros(1,2*mshift+1);
    
    istart = mid-wlenlf/2;
    iend = istart+wlenlf-1;
    sumstack12(1,:) = sum(stack12x(istart-mshift: iend-mshift, :));
    sumstack13(1,:) = sum(stack13x(istart-mshift: iend-mshift, :));
    sumstack32(1,:) = sum(stack32x(istart-mshift: iend-mshift, :));
    sumstack1sq(1,:) = sum(stackauto(1, istart: iend));
    sumstack3Bsq(1,:) = sum(stackauto(3, istart: iend));
    for m = -mshift:mshift
        sumstack2sq(1,m+mshift+1)=sum(stackauto(2, istart-m: iend-m)); %+m??? (yes).
        sumstack3sq(1,m+mshift+1)=sum(stackauto(3, istart-m: iend-m));
    end
    
    %An attempt to bypass glitches in data.  Min value of good data typically ~10^{-2}
    glitches=1.e-7;
    sumstack1sq=max(sumstack1sq,glitches);
    sumstack2sq=max(sumstack2sq,glitches);
    sumstack3sq=max(sumstack3sq,glitches);    % return maximum between A and B
    
    denomstack12n=realsqrt(sumstack1sq.*sumstack2sq);    % Real square root, An error is produced if X is negative
    denomstack13n=realsqrt(sumstack1sq.*sumstack3sq);
    denomstack32n=realsqrt(sumstack3Bsq.*sumstack2sq);
    
    sumstack12n=sumstack12./denomstack12n;   % suffix 'n' means normalized
    sumstack13n=sumstack13./denomstack13n;
    sumstack32n=sumstack32./denomstack32n;
    
    [xcmaxstack12n,imaxstack12]=max(sumstack12n,[],2);   %Integer-offset max cross-correlation
    [xcmaxstack13n,imaxstack13]=max(sumstack13n,[],2);   % along row, max cc val and index in each window
    [xcmaxstack32n,imaxstack32]=max(sumstack32n,[],2);
    
    %Parabolic fit:
    [xmaxstack12n,ymaxstack12n,aSTA12]=parabol(1,mshift,sumstack12n,imaxstack12); %Interpolated max cross-correlation
    [xmaxstack13n,ymaxstack13n,aSTA13]=parabol(1,mshift,sumstack13n,imaxstack13);
    [xmaxstack32n,ymaxstack32n,aSTA32]=parabol(1,mshift,sumstack32n,imaxstack32);
    
    %Center them
    imaxstack12cent=imaxstack12-mshift-1;  % "cent" is "centered"; imaxSTA12 is original 1: 2*mshift+1, corresponds to -mshift: mshift
    imaxstack13cent=imaxstack13-mshift-1;
    imaxstack32cent=imaxstack32-mshift-1;
    %%% NOTICE: the right order of a closed 3-sta pair is +13, -12, +32, where 13 means 1-->3
    iloopofflf=imaxstack13cent-imaxstack12cent+imaxstack32cent; %How well does the integer loop close?
    %
    xmaxstack12n=xmaxstack12n-mshift-1;
    xmaxstack13n=xmaxstack13n-mshift-1;
    xmaxstack32n=xmaxstack32n-mshift-1;
    loopoff=xmaxstack13n-xmaxstack12n+xmaxstack32n; %How well does the interpolated loop close?
    loopofflf(ifam) = loopoff;
    xcmaxAVEn=(xcmaxstack12n+xcmaxstack13n+xcmaxstack32n)/3;
    medxcmaxAVEn=median(xcmaxAVEn);
    xmaxstack12ntmplf=xmaxstack12n;    % tmp == temporary
    xmaxstack13ntmplf=xmaxstack13n;
    xmaxstack32ntmplf=xmaxstack32n;
    
    if xcmaxAVEn<xcmaxAVEnmin || abs(loopoff)>loopoffmax || isequal(abs(imaxstack12cent),mshift)...
       || isequal(abs(imaxstack13cent),mshift) || isequal(abs(imaxstack32cent),mshift)
        xmaxstack12ntmplf=mshift+1; xmaxstack13ntmplf=mshift+1; xmaxstack32ntmplf=mshift+1; %dummy them, if these criteria are met
        disp('WRONG! The basic criteria is not met');
    else
        xcmaxconprev=-99999.;  %used to be 0; not good with glitches
        %%% NOTE here are using offset imaxSTA12/13/32 without centering, i.e. [1, 2shift+1]
        % Usually, the loop is happened between imaxSTA12n +- floor(loopoffmax+1)
        % width of which is 2*floor(loopoffmax+1)
        %%% floor(2.5)=2; floor(-2.6)=-3
        for iSTA12 = max(1,imaxstack12-floor(loopoffmax+1)): ...
                     min(imaxstack12+floor(loopoffmax+1),2*mshift+1)
            for iSTA13 = max(1,imaxstack13-floor(loopoffmax+1)): ...
                         min(imaxstack13+floor(loopoffmax+1),2*mshift+1)
                ibangon = (mshift+1)-iSTA13+iSTA12;     %%% SEE NOTES #2019/03/17# page 66 to understand better
                %%% i.e., -mshift <= -iSTA13+iSTA12 <= mshift
                if ibangon >= 1 && ibangon <= 2*mshift+1
                    xcmaxcon=sumstack12n(iSTA12)+sumstack13n(iSTA13)+sumstack32n(ibangon);
                    if xcmaxcon > xcmaxconprev
                        xcmaxconprev=xcmaxcon;
                        iSTA12bang=iSTA12;
                        iSTA13bang=iSTA13;
                    end
                end
            end
        end
        %%% will result in the max xcmaxcon and corresponding iSTA12,
        %%% iSTA13, and save them into xcmaxconprev, iSTA12bang and iSTA13bang
        iSTA32bang=(mshift+1)-iSTA13bang+iSTA12bang;
        if abs(iSTA12bang-imaxstack12) <= loopoffmax && ...  %not sure if these 3 lines are satisfied automatically ...
           abs(iSTA13bang-imaxstack13) <= loopoffmax && ...  % SHOULD not be, i think, could be floor(loopoffmax+1) > loopoffmax
           abs(iSTA32bang-imaxstack32) <= loopoffmax && ...
           sumstack12n(iSTA12bang)+sumstack13n(iSTA13bang)+sumstack32n(iSTA32bang) ...
               >= 3*xcmaxAVEnmin   % xcmaxAVEnmin, predetermined
            %%% ALSO, sumsSTA12n(n,iSTA12bang) == sumsSTA12nn(iSTA12bang)
            
            xmaxstack12ntmplf=iSTA12bang-(mshift+1); %without interpolation this is just centering.
            xmaxstack13ntmplf=iSTA13bang-(mshift+1);
            xmaxstack32ntmplf=iSTA32bang-(mshift+1);
            
            %%% xcmaxAVEnbang is added by Chao, to distinguish from
            %%% xcmaxAVEn, because it is max average CC coef
            xcmaxAVEnbang=(sumstack12n(iSTA12bang)+sumstack13n(iSTA13bang)+sumstack32n(iSTA32bang))/3;
        end
    end
    
%     xmaxstack12ntmplf
%     xmaxstack13ntmplf
%     xmaxstack32ntmplf

    if xmaxstack13ntmplf-xmaxstack12ntmplf+xmaxstack32ntmplf ~=0
        disp('WRONG! Loopoff is not enclosed');
    else
        disp('Loopoff in lf is 0');
    end
    filtlf = [xmaxstack12ntmplf;  xmaxstack13ntmplf];
    
    % examine the widness of cc function
    cc = [sumstack12n; sumstack13n; sumstack32n];    
    off = [imaxstack12; imaxstack13; imaxstack32];
    
    for i = 1: 3 
        declf(1,i) = (2*cc(i,off(i))-cc(i,off(i)-1)-cc(i,off(i)+1))/2; % decline of cc for +-1 sample off
        declf(2,i) = (2*cc(i,off(i))-cc(i,off(i)-3)-cc(i,off(i)+3))/2; % 3 samples
        declf(3,i) = (2*cc(i,off(i))-cc(i,off(i)-6)-cc(i,off(i)+6))/2; % 6 samples
        declf(4,i) = (2*cc(i,off(i))-cc(i,off(i)-10)-cc(i,off(i)+10))/2;  % 10 samples
    end
    avedeclf = sum(declf,2)/3;
    declf(:,4) = avedeclf;
    
    % sumstackxxn is actually the cross-correlation between stations allowing mshift
    figure
    hold on
    p1=plot(-mshift:mshift,sumstack12n,'b-','linewidth', 1);
    p2=plot(-mshift:mshift,sumstack13n,'r-','linewidth', 1);
    p3=plot(-mshift:mshift,sumstack32n,'k-','linewidth', 1);
    
    p4=plot(imaxstack12cent,sumstack12n(imaxstack12),'bs','markersize',6);
    p5=plot(imaxstack13cent,sumstack13n(imaxstack13),'ro','markersize',6);
    p6=plot(imaxstack32cent,sumstack32n(imaxstack32),'kd','markersize',6);
    % plot(xmaxstack12n,ymaxstack12n,'bs','markersize',10);
    % plot(xmaxstack13n,ymaxstack13n,'ro','markersize',10);
    % plot(xmaxstack32n,ymaxstack32n,'kd','markersize',10);
    p7=plot(xmaxstack12ntmplf,sumstack12n(xmaxstack12ntmplf+mshift+1),'bs','markersize',10);
    p8=plot(xmaxstack13ntmplf,sumstack13n(xmaxstack13ntmplf+mshift+1),'ro','markersize',10);
    p9=plot(xmaxstack32ntmplf,sumstack32n(xmaxstack32ntmplf+mshift+1),'kd','markersize',10);
    legend([p1,p2,p3,p4,p5,p6,p7,p8,p9],{'CC TWKB-LZB(12)','CC TWKB-KLNB(13)','CC KLNB-LZB(32)',...
        'individual \Delta{t}_{12}',...
        'individual \Delta{t}_{13}', 'individual \Delta{t}_{32}',...
        'constrained \Delta{t}_{12}','constrained \Delta{t}_{13}',...
        'constrained \Delta{t}_{32}'});
    xlabel(strcat('Samples (on',{' '},num2str(sps),{' '}, 'sps)'),'fontsize',12);
    ylabel('Normalized CC Coef.','fontsize',12);
    text(-10,-0.4,'0.5-1.25 Hz','fontsize',12);
    box on;
    grid on;
    % vertical_cursors
    
    % display results
    disp(filtlf-filthf);
        
%     % save fig
%     fignm = strcat(rstpath, '/FIGS/tempfeff_LZB_enc_',fam,'.lf',num2str(lolf),'-',num2str(hilf),...
%         '.lw',num2str(wlenseclf),'.lo',num2str(loopoffmax),'.ccm', ...
%         num2str(xcmaxAVEnmin),'.',num2str(sps),'sps.eps');
%     print('-depsc',fignm);
    
    %%%% end for lf
    
        
    %%% save the filtering effect result
    savemat = [filthf filtlf filtlf-filthf];
    
%     PREFIX = strcat(fam,'.lf',num2str(lolf),'-',num2str(hilf),'.hf',num2str(lohf),'-',num2str(hihf),...
%         '.hw',num2str(wlensechf),'.lw',num2str(wlenseclf),'.lo',num2str(loopoffmax),'.ccm', ...
%         num2str(xcmaxAVEnmin),'.',num2str(sps),'sps');
%     fid = fopen(strcat(rstpath, '/MAPS/tempfeff_LZB_enc','_',PREFIX),'w+'); % template filtering effect for LZB trio enclosed
%     fprintf(fid,'%d %d %d \n',savemat');
%     fclose(fid);
    
    %%% save to summary array
    filtall = [filtall; savemat];
    decall = [decall; dechf declf];
    
end

% %%% save the width estimate
% PREFIX = strcat('allfam.lf',num2str(lolf),'-',num2str(hilf),'.hf',num2str(lohf),'-',num2str(hihf),...
%                 '.hw',num2str(wlensechf),'.lw',num2str(wlenseclf),'.lo',num2str(loopoffmax),'.ccm', ...
%                 num2str(xcmaxAVEnmin),'.',num2str(sps),'sps');
% fid = fopen(strcat(rstpath, '/MAPS/ccwidth_LZB','_',PREFIX),'w+'); % cc width estimate
% fprintf(fid,'%.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f \n',decall');
% fclose(fid);

%%
famintind = [2,3,6,7,8,9,10,11,12];
disp(nfampool(famintind, :));
medloffhf = median(abs(loopoffhf(famintind)));
medlofflf = median(abs(loopofflf(famintind)));

medfilt12 = median(filtall((famintind-1).*2+1, 3));
medfilt13 = median(filtall(famintind.*2, 3));





%%
figure
hold on
for ifam = 1: nfam
%     plot(decall((ifam-1)*4+1: ifam*4, 1), '--');
%     plot(decall((ifam-1)*4+1: ifam*4, 2), '-.');
%     plot(decall((ifam-1)*4+1: ifam*4, 3), ':');
    p(ifam)=plot(decall((ifam-1)*4+1: ifam*4, 4), 'o-');
end
plot(dec(1:4, 4), 'ko-','linew',1.5);
legend(p,nfampool,'location','northwest');
title('HF; decrease of CC wrt. samples from the center; 1,3,6,10');

%
figure
hold on
for ifam = 1: nfam
%     plot(decall((ifam-1)*4+1: ifam*4, 1), '--');
%     plot(decall((ifam-1)*4+1: ifam*4, 2), '-.');
%     plot(decall((ifam-1)*4+1: ifam*4, 3), ':');
    p(ifam)=plot(decall((ifam-1)*4+1: ifam*4, 8), 'o-');
end
plot(dec(1:4, 8), 'ko-','linew',1.5);
legend(p,nfampool,'location','northwest');
title('LF; decrease of CC wrt. samples from the center; 1,3,6,10');


%%
figure
hold on
p1=plot(abs(loopoffhf), 'ro-');
p2=plot(abs(loopofflf), 'bo-');
legend([p1,p2],{'loop HF','loop LF'},'location','best');
title('Abs. circuit of CC offsets');





