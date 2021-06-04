% function comp_medoff_mig_with_filtereffect
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script is used to predict the theoretical time offset between diff.
% frequencies at a station pair using some classic attenuation model, and
% compare them with the observed phase of cross-spectrum
%
%
%
% Chao Song, chaosong@princeton.edu
% First created date:   2021/05/21
% Last modified date:   2021/05/26
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% load some known parameters
clc
close all
clear

[scrsz, res] = pixelperinch(1);

% FLAG = 'PGC';
FLAG = 'TWKB';
disp(FLAG);

if strcmp(FLAG, 'TWKB')
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
    
    %frequency that has the largest amplitude according to the spectrum, NOT necessarily the best
    %measure, but can be served as a reference
    ftampmaxft = [
                  3.4375e+00   2.4219e+00   4.4531e+00
                  2.1875e+00   2.5000e+00   2.8906e+00
                  2.1875e+00   2.8125e+00   3.2031e+00
                  2.1094e+00   4.0625e+00   4.0625e+00
                  2.1094e+00   2.4219e+00   2.8906e+00
                  2.1094e+00   1.2500e+00   3.2031e+00
                  2.1875e+00   2.8125e+00   3.3594e+00
                  1.9531e+00   1.9531e+00   3.3594e+00
                  3.6719e+00   2.3438e+00   3.9062e+00
                  2.5000e+00   2.7344e+00   3.0078e+01
                  2.1094e+00   2.7344e+00   6.0156e+00
                  2.5000e+00   2.8125e+00   2.1094e+00
                  1.9531e+00   1.5625e+00   1.4844e+00
                  ];
    
    % lat and lon of the stations
    stas=['TWKB '
          'LZB  '
          'MGCB '];     % determine the trio and order
    nsta=size(stas,1);         %  number of stations
    
    dtdir = ('/home/data2/chaosong/matlab/allan/2004/JUL/');
    flist= ('datalist');
    outf = ('staloc.txt');
    stainfo = getstainfo(dtdir, flist, outf);
    ind=[11,3,4];
    stainfo = stainfo(ind,:);
    
    stalat = str2double(stainfo(:,2));
    stalon = str2double(stainfo(:,3));
    staloc = [stalon stalat zeros(length(stalon),1)];
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% VERSION 1 for the locations of lfe families, by inversion from hypoinverse
    %%% this is inverted from (0,0) of all fams, same order, location of control points
    %%% inverted from /home/data2/chaosong/Seisbasics/hypoinverse/lzbrst/chat00p.reffam
    cont = [
            -123.492667 48.451500 38.1400;
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
    
elseif strcmp(FLAG, 'PGC')
    nfampool = ['002';
                ];
    nfam = size(nfampool,1);
    
    stas=['PGC  '
          'SSIB '
          'SILB '];
    nsta=size(stas,1);         %  number of stations

    dtdir = ('/home/data2/chaosong/matlab/allan/2004/JUL/');
    flist= ('datalist');
    outf = ('staloc.txt');
    stainfo = getstainfo(dtdir, flist, outf);
    ind=[5,8,6];
    stainfo = stainfo(ind,:);
    
    stalat = str2double(stainfo(:,2));
    stalon = str2double(stainfo(:,3));
    staloc = [stalon stalat zeros(length(stalon),1)];
    
    % location of control fam, this is infered from 'mig_linear_fit_v3' and the detections
    cont = [
            -123.5850 48.4367 36.93;
            ];
end

%% Azimi's Law in my paper
for ifam = 1: nfam
    
    Qs = 200;

    fam = nfampool(ifam,:);
    [starel(:,1),starel(:,2)] = absloc2relaloc(stalon,stalat,cont(ifam,1),cont(ifam,2));
    starel(:,3) = ones(size(starel,1),1)*cont(ifam,3);
    Dist = sqrt(sum(starel.^2, 2));  % distance from source to any single station

    %let's say the velocity model gives us the vel. at 1 Hz
    c0 = 3.84;
    f0 = 1;
    om0 = 2*pi*f0;
    Fs = 80;

    fhf = (1.25+6.5)/2;
    flf = (0.5+1.25)/2;
    omhf = 2*pi*fhf;
    omlf = 2*pi*flf;
    
    chf = c0*(1+log(omhf/om0)/(pi*Qs));
%     chf = 3.84;
    dtfreq_12 = (Dist(2) - Dist(1))*log(omhf/omlf)/(pi*Qs*chf); 
    dtfreq_13 = (Dist(3) - Dist(1))*log(omhf/omlf)/(pi*Qs*chf);
    dtfreq_23 = (Dist(3) - Dist(2))*log(omhf/omlf)/(pi*Qs*chf);
    dnfreq_12 = dtfreq_12 * Fs;  
    dnfreq_13 = dtfreq_13 * Fs;
    dnfreq_23 = dtfreq_23 * Fs;
    
    %if fix a frequency, and analyse the tshift/nshift, or normalized phase as a function of 
    %frequency, meanning that you align the records at the reference frequency at two stations, and
    %then analyze the time offset at other frequency wrt the aligned frequency 
    falign = [(0.5+1.25)/2 2.5 (1.25+6.5)/2];
%     falign = 2.5;
    freq = 0.01: 0.02: 8;
    
    f.fig = figure;
    f.fig.Renderer = 'painters';
    widin = 6;  % maximum width allowed is 8.5 inches
    htin = 8;   % maximum height allowed is 11 inches
    set(f.fig,'Position',[1*scrsz(3)/20 scrsz(4)/10 widin*res htin*res]);
    nrow = 2;
    ncol = 1;
    for isub = 1:nrow*ncol
        f.ax(isub) = subplot(nrow,ncol,isub);
        grid(f.ax(isub), 'on');
        f.ax(isub).Box = 'on';
    end
    
    % reposition
    set(f.ax(2), 'position', [ 0.12, 0.08, 0.8, 0.4]);
    set(f.ax(1), 'position', [ 0.12, 0.55, 0.8, 0.4]);

    lines = {'--'; '-'; '-.'};
    linew = [0.8 1.5 0.8];
    
    for i = 1: length(falign)
        omalign = 2*pi*falign(i);
        omlf = 2*pi*freq;
        calign = c0*(1+log(omalign/om0)/(pi*Qs));
        dtfreq_12 = (Dist(2) - Dist(1))*log(omalign./omlf)./(pi*Qs*calign);
        dtfreq_13 = (Dist(3) - Dist(1))*log(omalign./omlf)./(pi*Qs*calign);
        dtfreq_23 = (Dist(3) - Dist(2))*log(omalign./omlf)./(pi*Qs*calign);
        
        dnfreq_12 = dtfreq_12 * Fs;
        dnfreq_13 = dtfreq_13 * Fs;
        dnfreq_23 = dtfreq_23 * Fs;
        
        phalagfreq_12 = dtfreq_12 .* omlf /pi;
        phalagfreq_13 = dtfreq_13 .* omlf /pi;
        phalagfreq_23 = dtfreq_23 .* omlf /pi;
        
        
        %if dtfreq_ij/dnfreq_ij/phalagfreq_ij is positive, according to the defination in the eqn and
        %physics behind, it means if you align the reference frequency at two stations, the requested
        %frequency would arrive later at station j relative to station i
        ax1 = f.ax(1);
        hold(ax1, 'on');
        p((i-1)*3+1) = plot(ax1,freq,phalagfreq_12,lines{i},'linewidth',linew(i),'color','r');
        p((i-1)*3+2) = plot(ax1,freq,phalagfreq_13,lines{i},'linewidth',linew(i),'color','b');
        p((i-1)*3+3) = plot(ax1,freq,phalagfreq_23,lines{i},'linewidth',linew(i),'color','k');
        if strcmp(FLAG, 'TWKB')
            lgdstr{(i-1)*3+1} = strcat({'TWKB-LZB, '},num2str(falign(i)),{' Hz'});
            lgdstr{(i-1)*3+2} = strcat({'TWKB-MGCB, '},num2str(falign(i)),{' Hz'});
            lgdstr{(i-1)*3+3} = strcat({'LZB-MGCB, '},num2str(falign(i)),{' Hz'});
        elseif strcmp(FLAG, 'PGC')
            lgdstr{(i-1)*3+1} = strcat({'PGC-SSIB, '},num2str(falign(i)),{' Hz'});
            lgdstr{(i-1)*3+2} = strcat({'PGC-SILB, '},num2str(falign(i)),{' Hz'});
            lgdstr{(i-1)*3+3} = strcat({'SSIB-SILB, '},num2str(falign(i)),{' Hz'});
        end
%         hold(ax1, 'off');
        
        ax2 = f.ax(2);
        hold(ax2,'on');
        p((i-1)*3+1) = plot(ax2,freq,dnfreq_12,lines{i},'linewidth',linew(i),'color','r');
        p((i-1)*3+2) = plot(ax2,freq,dnfreq_13,lines{i},'linewidth',linew(i),'color','b');
        p((i-1)*3+3) = plot(ax2,freq,dnfreq_23,lines{i},'linewidth',linew(i),'color','k');
        hold(ax2, 'off');
    end
    
    for isub = 1: nrow*ncol
        text(f.ax(isub),0.5, 0.9,'if positive, 2nd station of pair is lagging','fontsize',10,'unit',...
            'normalized','Horizontalalignment','center');
        text(f.ax(isub),0.5, 0.8,'if negative, 1st station of pair is lagging','fontsize',10,'unit',...
            'normalized','Horizontalalignment','center');
        lgd = legend(f.ax(isub),p,{string(lgdstr{1,1}),string(lgdstr{1,2}),string(lgdstr{1,3}),...
                     string(lgdstr{1,4}),string(lgdstr{1,5}),string(lgdstr{1,6}),string(lgdstr{1,7}),...
                     string(lgdstr{1,8}),string(lgdstr{1,9})},'location','south',...
                     'fontsize',6,'numcolumns',3);
        % make background transparent   
        set(lgd.BoxFace, 'ColorType','truecoloralpha', 'ColorData',uint8(255*[1;1;1;.8]));  
         
        xlabel(f.ax(isub),'Frequency (Hz)');
    end
    titstr = strcat(FLAG,{', '},fam,{', '},num2str(Fs),{' Hz, '},{'Qs: '},...
                    num2str(Qs),{', vel-1hz: '},num2str(c0),{' km/s'});
    title(ax1,titstr,'fontsize',12);
    ylabel(ax1,'Normalized Phase lag (phase/pi)');

    ylabel(ax2,'Lag samples');
    
end


