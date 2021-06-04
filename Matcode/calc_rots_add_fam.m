%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This is the script to call the function calc_rots_noresp to compute the 
% rotation parameters for some additional fams located at the southwestern
% part of our research area, and which fam would give a best stacked template
% which is a indictor of rotation, splitting correction reliability, then use
% this selected family for detections, to examine whether HF and LF have a 
% distinct spatial detection density pattern 
%       
%      
%
% NOTES:
%   1. there are fams that the stacked templates are pretty low-amplitude,
%   then the maximum offset allowed in CC would be necessary, if the result
%   reaches the boundary, then it is necessary to decrease the threshold.
%   However, check the N and E component to see if the splitting should be
%   small or large. For example, fam 010 should allow splitthres <=8, but
%   others could works fine when splitthres <=12, and usually this offset
%   ranges from 0~9.
%
%
% By Chao Song, chaosong@princeton.edu
% First created date:   2020/05/25
% Last modified date:   2020/05/25
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
datapath = strcat(workpath,'/data-no-resp');
rstpath = strcat(datapath,'/split_chao/lzbtrio/');

% families that are located at the southwestern part of our research area
% fampool = ['001';
%            '019';
%            '021';
%            '045';
%            '076';
%            '176'];

% from the raw stacking of LFEs, 231 is apparently not a good family
fampool = ['006';
           '015';
           '158';
           '231';
           '234'];
       
fampool = ['002';
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
       
nfam = length(fampool);
        
%%% corresponds to PERMROTS
PERMSTA=['PGC'        % permanent station names
    'LZB'];
POLSTA=['SSIB '           % polaris station names
    'SILB '
    'KLNB '
    'MGCB '
    'TWKB '];

stas=['PGC  '
    'SSIB '
    'SILB '
    'LZB  '
    'TWKB '
    'MGCB '
    'KLNB '];
% number of used stations
nsta=size(stas,1);         %  number of stations


%% calculation
ROTSALL = [];

remake = 1;  % whether to calculate the rots again. 0/1

lo = 0.5;
hi = 6.5;
sps = 40;
bef = 25;
aft = 35;

% [ROTS]=calc_rots_noresp(fam,lo,hi,sps,splitoff,staoff,mainsta)

if remake   % re-calculate all results
%     %%% fam 001
%     fam = '001';
%     splitoff = 12;
%     staoff = 10;
%     mainsta = 4; 
%     [ROTS]=calc_rots_noresp(fam,lo,hi,sps,bef,aft,splitoff,staoff,mainsta);
%     % save all in one file
% %     ROTSALL((1-1)*nsta+1:1*nsta, :) = ROTS;
%     % save results for each family
%     PREFIX = strcat(fam,'rot_para',num2str(lo),'-',num2str(hi),'_',num2str(bef),'_',num2str(aft),...
%                     stas(mainsta,:));
%     fid = fopen(strcat(rstpath, PREFIX,'_',num2str(sps),'sps',num2str(nsta),'stas'),'w+');
%     fprintf(fid,'%d %d %d %d \n',ROTS');
%     fclose(fid);
%     
%     %%% fam 019
%     fam = '019';
%     splitoff = 12;
%     staoff = 8;
%     mainsta = 4; 
%     [ROTS]=calc_rots_noresp(fam,lo,hi,sps,bef,aft,splitoff,staoff,mainsta);
%     % save all in one file
% %     ROTSALL((1-1)*nsta+1:1*nsta, :) = ROTS;
%     % save results for each family
%     PREFIX = strcat(fam,'rot_para',num2str(lo),'-',num2str(hi),'_',num2str(bef),'_',num2str(aft),...
%                     stas(mainsta,:));
%     fid = fopen(strcat(rstpath, PREFIX,'_',num2str(sps),'sps',num2str(nsta),'stas'),'w+');
%     fprintf(fid,'%d %d %d %d \n',ROTS');
%     fclose(fid);
%     
%     %%% fam 021
%     fam = '021';
%     splitoff = 12;
%     staoff = 8;
%     mainsta = 4; 
%     [ROTS]=calc_rots_noresp(fam,lo,hi,sps,bef,aft,splitoff,staoff,mainsta);
%     % save all in one file
% %     ROTSALL((1-1)*nsta+1:1*nsta, :) = ROTS;
%     % save results for each family
%     PREFIX = strcat(fam,'rot_para',num2str(lo),'-',num2str(hi),'_',num2str(bef),'_',num2str(aft),...
%                     stas(mainsta,:));
%     fid = fopen(strcat(rstpath, PREFIX,'_',num2str(sps),'sps',num2str(nsta),'stas'),'w+');
%     fprintf(fid,'%d %d %d %d \n',ROTS');
%     fclose(fid);
% 
%     %%% fam 045
%     fam = '045';
%     splitoff = 12;
%     staoff = 8;
%     mainsta = 4; 
%     [ROTS]=calc_rots_noresp(fam,lo,hi,sps,bef,aft,splitoff,staoff,mainsta);
%     % save all in one file
% %     ROTSALL((1-1)*nsta+1:1*nsta, :) = ROTS;
%     % save results for each family
%     PREFIX = strcat(fam,'rot_para',num2str(lo),'-',num2str(hi),'_',num2str(bef),'_',num2str(aft),...
%                     stas(mainsta,:));
%     fid = fopen(strcat(rstpath, PREFIX,'_',num2str(sps),'sps',num2str(nsta),'stas'),'w+');
%     fprintf(fid,'%d %d %d %d \n',ROTS');
%     fclose(fid);
%     
%     %%% fam 076
%     fam = '076';
%     splitoff = 15;
%     staoff = 8;
%     mainsta = 4; 
%     [ROTS]=calc_rots_noresp(fam,lo,hi,sps,bef,aft,splitoff,staoff,mainsta);
%     % save all in one file
% %     ROTSALL((1-1)*nsta+1:1*nsta, :) = ROTS;
%     % save results for each family
%     PREFIX = strcat(fam,'rot_para',num2str(lo),'-',num2str(hi),'_',num2str(bef),'_',num2str(aft),...
%                     stas(mainsta,:));
%     fid = fopen(strcat(rstpath, PREFIX,'_',num2str(sps),'sps',num2str(nsta),'stas'),'w+');
%     fprintf(fid,'%d %d %d %d \n',ROTS');
%     fclose(fid);
%     
%     %%% fam 176
%     fam = '176';
%     splitoff = 9;
%     staoff = 8;
%     mainsta = 4; 
%     [ROTS]=calc_rots_noresp(fam,lo,hi,sps,bef,aft,splitoff,staoff,mainsta);
%     % save all in one file
% %     ROTSALL((1-1)*nsta+1:1*nsta, :) = ROTS;
%     % save results for each family
%     PREFIX = strcat(fam,'rot_para',num2str(lo),'-',num2str(hi),'_',num2str(bef),'_',num2str(aft),...
%                     stas(mainsta,:));
%     fid = fopen(strcat(rstpath, PREFIX,'_',num2str(sps),'sps',num2str(nsta),'stas'),'w+');
%     fprintf(fid,'%d %d %d %d \n',ROTS');
%     fclose(fid);

    %%% fam 006
    fam = '006';
    splitoff = 10;
    staoff = 8;
    mainsta = 4; 
    [ROTS]=calc_rots_noresp(fam,lo,hi,sps,bef,aft,splitoff,staoff,mainsta);
    % save all in one file
%     ROTSALL((1-1)*nsta+1:1*nsta, :) = ROTS;
    % save results for each family
    PREFIX = strcat(fam,'rot_para',num2str(lo),'-',num2str(hi),'_',num2str(bef),'_',num2str(aft),...
                    stas(mainsta,:));
    fid = fopen(strcat(rstpath, PREFIX,'_',num2str(sps),'sps',num2str(nsta),'stas'),'w+');
    fprintf(fid,'%d %d %d %d \n',ROTS');
    fclose(fid);

%     %%% fam 015
%     fam = '015';
%     splitoff = 12;
%     staoff = 8;
%     mainsta = 4; 
%     [ROTS]=calc_rots_noresp(fam,lo,hi,sps,bef,aft,splitoff,staoff,mainsta);
%     % save all in one file
% %     ROTSALL((1-1)*nsta+1:1*nsta, :) = ROTS;
%     % save results for each family
%     PREFIX = strcat(fam,'rot_para',num2str(lo),'-',num2str(hi),'_',num2str(bef),'_',num2str(aft),...
%                     stas(mainsta,:));
%     fid = fopen(strcat(rstpath, PREFIX,'_',num2str(sps),'sps',num2str(nsta),'stas'),'w+');
%     fprintf(fid,'%d %d %d %d \n',ROTS');
%     fclose(fid);
% 
%     %%% fam 158
%     fam = '158';
%     splitoff = 12;
%     staoff = 8;
%     mainsta = 4; 
%     [ROTS]=calc_rots_noresp(fam,lo,hi,sps,bef,aft,splitoff,staoff,mainsta);
%     % save all in one file
% %     ROTSALL((1-1)*nsta+1:1*nsta, :) = ROTS;
%     % save results for each family
%     PREFIX = strcat(fam,'rot_para',num2str(lo),'-',num2str(hi),'_',num2str(bef),'_',num2str(aft),...
%                     stas(mainsta,:));
%     fid = fopen(strcat(rstpath, PREFIX,'_',num2str(sps),'sps',num2str(nsta),'stas'),'w+');
%     fprintf(fid,'%d %d %d %d \n',ROTS');
%     fclose(fid);
%     
    %%% fam 234
    fam = '234';
    splitoff = 12;
    staoff = 8;
    mainsta = 4; 
    [ROTS]=calc_rots_noresp(fam,lo,hi,sps,bef,aft,splitoff,staoff,mainsta);
    % save all in one file
%     ROTSALL((1-1)*nsta+1:1*nsta, :) = ROTS;
    % save results for each family
    PREFIX = strcat(fam,'rot_para',num2str(lo),'-',num2str(hi),'_',num2str(bef),'_',num2str(aft),...
                    stas(mainsta,:));
    fid = fopen(strcat(rstpath, PREFIX,'_',num2str(sps),'sps',num2str(nsta),'stas'),'w+');
    fprintf(fid,'%d %d %d %d \n',ROTS');
    fclose(fid);
    
%     %%% fam 231
%     fam = '231';
%     splitoff = 12;
%     staoff = 10;
%     mainsta = 4; 
%     [ROTS]=calc_rots_noresp(fam,lo,hi,sps,bef,aft,splitoff,staoff,mainsta);
%     % save all in one file
% %     ROTSALL((1-1)*nsta+1:1*nsta, :) = ROTS;
%     % save results for each family
%     PREFIX = strcat(fam,'rot_para',num2str(lo),'-',num2str(hi),'_',num2str(bef),'_',num2str(aft),...
%                     stas(mainsta,:));
%     fid = fopen(strcat(rstpath, PREFIX,'_',num2str(sps),'sps',num2str(nsta),'stas'),'w+');
%     fprintf(fid,'%d %d %d %d \n',ROTS');
%     fclose(fid);
    
else    % directly read results
    
    %%% fam 001
    fam = '001';
    mainsta = 4; 
    PREFIX = strcat(fam,'rot_para',num2str(lo),'-',num2str(hi),'_',num2str(bef),'_',num2str(aft),...
                    stas(mainsta,:));
    [ROTS]=load(strcat(rstpath, PREFIX,'_',num2str(sps),'sps',num2str(nsta),'stas'));
    % save all in one file
    ROTSALL = [ROTSALL; ROTS];
    
    %%% fam 019
    fam = '019';
    mainsta = 4; 
    PREFIX = strcat(fam,'rot_para',num2str(lo),'-',num2str(hi),'_',num2str(bef),'_',num2str(aft),...
                    stas(mainsta,:));
    [ROTS]=load(strcat(rstpath, PREFIX,'_',num2str(sps),'sps',num2str(nsta),'stas'));
    ROTSALL = [ROTSALL; ROTS];
    
    %%% fam 021
    fam = '021';
    mainsta = 4; 
    PREFIX = strcat(fam,'rot_para',num2str(lo),'-',num2str(hi),'_',num2str(bef),'_',num2str(aft),...
                    stas(mainsta,:));
    [ROTS]=load(strcat(rstpath, PREFIX,'_',num2str(sps),'sps',num2str(nsta),'stas'));
    ROTSALL = [ROTSALL; ROTS];
    
    %%% fam 045
    fam = '045';
    mainsta = 4; 
    PREFIX = strcat(fam,'rot_para',num2str(lo),'-',num2str(hi),'_',num2str(bef),'_',num2str(aft),...
                    stas(mainsta,:));
    [ROTS]=load(strcat(rstpath, PREFIX,'_',num2str(sps),'sps',num2str(nsta),'stas'));
    ROTSALL = [ROTSALL; ROTS];
    
    %%% fam 076
    fam = '076';
    mainsta = 4; 
    PREFIX = strcat(fam,'rot_para',num2str(lo),'-',num2str(hi),'_',num2str(bef),'_',num2str(aft),...
                    stas(mainsta,:));
    [ROTS]=load(strcat(rstpath, PREFIX,'_',num2str(sps),'sps',num2str(nsta),'stas'));
    ROTSALL = [ROTSALL; ROTS];
    
    %%% fam 176
    fam = '176';
    mainsta = 4; 
    PREFIX = strcat(fam,'rot_para',num2str(lo),'-',num2str(hi),'_',num2str(bef),'_',num2str(aft),...
                    stas(mainsta,:));
    [ROTS]=load(strcat(rstpath, PREFIX,'_',num2str(sps),'sps',num2str(nsta),'stas'));
    ROTSALL = [ROTSALL; ROTS];
    
    %%% fam 176
    fam = '176';
    mainsta = 4; 
    PREFIX = strcat(fam,'rot_para',num2str(lo),'-',num2str(hi),'_',num2str(bef),'_',num2str(aft),...
                    stas(mainsta,:));
    [ROTS]=load(strcat(rstpath, PREFIX,'_',num2str(sps),'sps',num2str(nsta),'stas'));
    ROTSALL = [ROTSALL; ROTS];
    
    %%% fam 015
    fam = '015';
    mainsta = 4; 
    PREFIX = strcat(fam,'rot_para',num2str(lo),'-',num2str(hi),'_',num2str(bef),'_',num2str(aft),...
                    stas(mainsta,:));
    [ROTS]=load(strcat(rstpath, PREFIX,'_',num2str(sps),'sps',num2str(nsta),'stas'));
    ROTSALL = [ROTSALL; ROTS];
    
    %%% fam 158
    fam = '158';
    mainsta = 4; 
    PREFIX = strcat(fam,'rot_para',num2str(lo),'-',num2str(hi),'_',num2str(bef),'_',num2str(aft),...
                    stas(mainsta,:));
    [ROTS]=load(strcat(rstpath, PREFIX,'_',num2str(sps),'sps',num2str(nsta),'stas'));
    ROTSALL = [ROTSALL; ROTS];
    
    %%% fam 234
    fam = '234';
    mainsta = 4; 
    PREFIX = strcat(fam,'rot_para',num2str(lo),'-',num2str(hi),'_',num2str(bef),'_',num2str(aft),...
                    stas(mainsta,:));
    [ROTS]=load(strcat(rstpath, PREFIX,'_',num2str(sps),'sps',num2str(nsta),'stas'));
    ROTSALL = [ROTSALL; ROTS];
    
    %%% fam 231
    fam = '231';
    mainsta = 4; 
    PREFIX = strcat(fam,'rot_para',num2str(lo),'-',num2str(hi),'_',num2str(bef),'_',num2str(aft),...
                    stas(mainsta,:));
    [ROTS]=load(strcat(rstpath, PREFIX,'_',num2str(sps),'sps',num2str(nsta),'stas'));
    ROTSALL = [ROTSALL; ROTS];
    
end

%% save resuls of all fam
% PREFIX = strcat('rot_para_add6fam',num2str(lo),'-',num2str(hi),'_',num2str(bef),'_',num2str(aft),'LZB');
PREFIX = strcat('rot_para_SW4fam',num2str(lo),'-',num2str(hi),'_',num2str(bef),'_',num2str(aft),'LZB');
fid = fopen(strcat(rstpath,PREFIX,'_',num2str(sps),'sps',num2str(nsta),'stas'),'w+');
fprintf(fid,'%d %d %d %d \n',ROTSALL');
fclose(fid);











