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
scrsz=get(0,'ScreenSize');   % for example, scrsz=[1,1,1680,1050]
wid=scrsz(3);       % width
hite=scrsz(4);       % height

% set path
workpath = getenv('ALLAN');
% datapath = workpath;
datapath = strcat(workpath,'/data-no-resp');
temppath = strcat(datapath, '/templates/LZBtrio/');



% the newest families from Bostock are different, so select the closest
% one as new family pool, also selected family have the most LFEs
fampool = ['001';
           '019';
           '021';
           '045';
           '076';
           '176'];

% from a preliminary detection result, 001 seems to have most HF detections       
fampool = ['001'];

% from the obtained templates, 231 is apparently not a good family
% 234 is bit weird on timing, but the shape of dipole seems good enough
% overall, 234 and 158 is comparable, 015 is less good,
% given the their locations, let's try 234 and 158 together!
fampool = ['006';
           '015';
           '158';
           '231';
           '234'];

fampool = ['006';
           '001'];
       
nfam = length(fampool);

stas=['PGC  '
    'SSIB '
    'SILB '
    'LZB  '
    'TWKB '
    'MGCB '
    'KLNB '];

nsta=size(stas,1);         %  number of stations
templensec = 60;

sps = 40;   % sampling rate you want
ccmethod = 2;  % only using lfes passed the cc thresholds (ccmethod=1); using all lfes (ccmethod=2)
plflag = 1;  % whether to plot (1/0)

%% caculate for each fam
for ifam = 1: size(fampool,1)
    
    fam = fampool(ifam,:);    
        
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
    
    
%     % make bandpassed stacks
%     bplo = 0.1;
%     bphi = 15;
%     [dstack, ccstack] = mk_bptemp_LZB(fam,sps,bplo,bphi,ccmethod,plflag);
%     % write into files
%     for ista = 1: nsta
%         fid = fopen(strcat(temppath, fam, '_', strtrim(stas(ista, :)), '_', ...
%             num2str(sps), 'sps_', num2str(templensec), 's_', ...
%             'BPDS_', 'opt_Nof_Non_Chao'), 'w+');  % bandpassed direct stack, no filter, no norm
%         
%         fprintf(fid, '%f \n', dstack(ista, :)');
%         fclose(fid);
%         
%         fid = fopen(strcat(temppath, fam, '_', strtrim(stas(ista, :)), '_', ...
%             num2str(sps), 'sps_', num2str(templensec), 's_', ...
%             'BPCCS_', 'opt_Nof_Non_Chao'), 'w+');  % bandpassed cc stack, no filter, no norm
%         fprintf(fid, '%f \n', ccstack(ista, :)');
%         fclose(fid);       
%     end
    
end












