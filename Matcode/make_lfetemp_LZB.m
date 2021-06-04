%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This is the script to make the lfe templates for LZB trio of all fams
%
%
% By Chao Song, chaosong@princeton.edu
% First created date:   2019/11/13
% Last modified date:   2019/11/13

%% Initialization
format short e   % Set the format to 5-digit floating point
clear
close all

set(0,'DefaultFigureVisible','on');
% set(0,'DefaultFigureVisible','off');   % switch to show the plots or not

[scrsz, res] = pixelperinch(1);

% set path
workpath = getenv('ALLAN');
% datapath = workpath;
datapath = strcat(workpath,'/data-no-resp');
temppath = strcat(datapath, '/templates/LZBtrio/');

% old LFE family pool, from Yajun's selections
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

% the newest families from Bostock are different, so select the closest
% one as new family pool, also selected family have the most LFEs
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

filterc = [2 -3;
           -7 -10;
           -6 -12;
           1 -7;
           -6 -3;
           -5 -10;
           -3 -10;
           -6 -7;
           -5 -10;
           -5 -8;
           -3 -9];

stas=['PGC  '
    'SSIB '
    'SILB '
    'LZB  '
    'TWKB '
    'MGCB '
    'KLNB '];

nsta=size(stas,1);         %  number of stations
templensec = 60;

sps = 80;   % sampling rate you want
ccmethod = 2;  % only using lfes passed the cc thresholds (ccmethod=1); using all lfes (ccmethod=2)
plflag = 1;  % whether to plot (1/0)
remake = 0;  % whether to recalculate

%% caculate for each fam
for ifam = 1: length(nfampool)
    
    fam = nfampool(ifam,:);
    
    if remake
        % make broadband stacks
        [dstack, ccstack] = mk_bbtemp_LZB(fam,sps,ccmethod,plflag);
        % write into files
        for ista = 1: nsta
            fid = fopen(strcat(temppath, fam, '_', strtrim(stas(ista, :)), '_', ...
                num2str(sps), 'sps_', num2str(templensec), 's_', ...
                'BBDS_', 'opt_Nof_Non_Chao'), 'w+');  % broadband direct stack, no filter, no norm
            
            fprintf(fid, '%f \n', dstack(ista, :)');
            fclose(fid);
            
            fid = fopen(strcat(temppath, fam, '_', strtrim(stas(ista, :)), '_', ...
                num2str(sps), 'sps_', num2str(templensec), 's_', ...
                'BBCCS_', 'opt_Nof_Non_Chao'), 'w+');  % broadband cc stack, no filter, no norm
            fprintf(fid, '%f \n', ccstack(ista, :)');
            fclose(fid);
        end
        
    else
        for ista = 1: nsta
            fname = strcat(temppath, fam, '_', strtrim(stas(ista, :)), '_', num2str(sps), 'sps_', ...
                num2str(templensec), 's_', 'BBDS_', 'opt_Nof_Non_Chao');
            dstack(ista,:) = load(fname);
            fname = strcat(temppath, fam, '_', strtrim(stas(ista, :)), '_', num2str(sps), 'sps_', ...
                num2str(templensec), 's_', 'BBCCS_', 'opt_Nof_Non_Chao');
            ccstack(ista,:) = load(fname);
        end
        
        if plflag
            templensec = 60;
            templen = templensec * sps;
            mid = templen/2;
            %%% plot the 2 kinds of template together
            f.fig=figure;
            f.fig.Renderer='Painters';
            widin = 7;  % maximum width allowed is 8.5 inches
            htin = 9;   % maximum height allowed is 11 inches
            set(f.fig,'Position',[scrsz(1)+1*scrsz(3)/10 scrsz(2)+scrsz(4)/20 widin*res htin*res]);
                   
            for ista = 1: nsta
                f.ax(ista) = subplot(nsta,1,ista);
                ax = f.ax(ista);
                plot(ax,dstack(ista, mid-5*sps: mid+5*sps), 'k', 'linewidth', 1); hold on
                plot(ax,ccstack(ista, mid-5*sps: mid+5*sps), 'r', 'linewidth', 1);
                title(ax,strcat('Stacked template with/without CC at station', {' '}, strtrim(stas(ista, :))));
                legend(ax,'Direct stack', 'CC stack');
                xlabel(ax,'Samples');
                ylabel(ax,'Amplitude');
                ax.Box = 'on';
                grid(ax, 'on');
            end
            text(f.ax(1),0.05,0.8,fam,'FontSize',10,'unit','normalized');
            text(f.ax(4),0.05,0.2,num2str(filterc(ifam,1)),'FontSize',12,'unit','normalized');
            text(f.ax(6),0.05,0.2,num2str(filterc(ifam,2)),'FontSize',12,'unit','normalized');
            
            print(f.fig,'-depsc2',strcat(temppath,'/',fam,'.broadbandstack.eps'));
            
        end
        
    end
    
    
%     %%% make bandpassed stacks
%     bplo = 0.1;
%     bphi = 15;
%     if remake
%         [dstack, ccstack] = mk_bptemp_LZB(fam,sps,bplo,bphi,ccmethod,plflag);
%         % write into files
%         for ista = 1: nsta
%             fid = fopen(strcat(temppath, fam, '_', strtrim(stas(ista, :)), '_', ...
%                 num2str(sps), 'sps_', num2str(templensec), 's_', ...
%                 'BPDS_', 'opt_Nof_Non_Chao'), 'w+');  % bandpassed direct stack, no filter, no norm
%             
%             fprintf(fid, '%f \n', dstack(ista, :)');
%             fclose(fid);
%             
%             fid = fopen(strcat(temppath, fam, '_', strtrim(stas(ista, :)), '_', ...
%                 num2str(sps), 'sps_', num2str(templensec), 's_', ...
%                 'BPCCS_', 'opt_Nof_Non_Chao'), 'w+');  % bandpassed cc stack, no filter, no norm
%             fprintf(fid, '%f \n', ccstack(ista, :)');
%             fclose(fid);
%         end
%         
%     else
%         for ista = 1: nsta
%             fname = strcat(temppath, fam, '_', strtrim(stas(ista, :)), '_', num2str(sps), 'sps_', ...
%                 num2str(templensec), 's_', 'BPDS_', 'opt_Nof_Non_Chao');
%             dstack(ista,:) = load(fname);
%             fname = strcat(temppath, fam, '_', strtrim(stas(ista, :)), '_', num2str(sps), 'sps_', ...
%                 num2str(templensec), 's_', 'BPCCS_', 'opt_Nof_Non_Chao');
%             ccstack(ista,:) = load(fname);
%         end
%     end
%     
%     if plflag
%         templensec = 60;
%         templen = templensec * sps;
%         mid = templen/2;
%         %%% plot the 2 kinds of template together
%         figure('Position',[scrsz(3)/10 scrsz(4)/10 3*scrsz(3)/5 4*scrsz(4)/5]);
%         for ista = 1: nsta
%             subplot(nsta,2,2*ista-1);
%             plot(dstack(ista, mid-5*sps: mid+5*sps), 'k', 'linewidth', 1); hold on
%             plot(ccstack(ista, mid-5*sps: mid+5*sps), 'r', 'linewidth', 1);
%             title(strcat('Stacked template with/without CC at station', {' '}, strtrim(stas(ista, :))));
%             legend('Direct stack', 'CC stack');
%             xlabel('Samples');
%             ylabel('Amplitude');
%             box on
%             grid on
%         end
%         
%         for ista = 1: nsta
%             subplot(nsta,2,2*ista);
%             plot(dstackno(ista, mid-5*sps: mid+5*sps), 'k', 'linewidth', 1); hold on
%             plot(ccstackno(ista, mid-5*sps: mid+5*sps), 'r', 'linewidth', 1);
%             title(strcat('Normalized stacked template with/without CC at station', {' '}, ...
%                 strtrim(stas(ista, :))));
%             legend('Direct stack', 'CC stack');
%             xlabel('Samples');
%             ylabel('Normalized amplitude');
%             box on
%             grid on
%         end
%         
%         print('-depsc2',strcat(temppath,'/',fam,'.bandpassstack.eps'));
%     end

end    
    
    
    
    
    
    
    
    
    
