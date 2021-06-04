%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% THIS is the script to merge the detections from all families IF you want 
% to set a distance cutoff first, then merge fams to remove duplicates base
% on time and CC coef, so run this AFTER hypoinverse.
% The input is already been removed off the double counting. Since for each
% fam, it wiped put the same data, so the detections should still have a lot of
% duplicates, preserve the ones within the one duration that have the highest 
% CC value. In the mean time, detections in each family have a different
% reference (0,0) point, so the merging can also convert all to the same
% frame (similar to changing centroid).
%
%   2020/12/07, try to exclude some fams as requested, for testing purpose,
%               bc the filtering effect from some fams is slightly different
%               from others, we need to make sure the result is not biased 
%               those 'abnormal' filtering effects
%
% NOTE :
%   1*. This should be ran AFTER hypoinverse, before do migration analyzing!!!
%       
%        
%   2. run with fcflag=0 (DEFAULT) would save the no-duplicate results for 
%       plotaddcheck_LZB_allfam and figures; run with fcflag=1 would correct 
%       the filtering effect for no-duplicate results to generate and save new
%        figures, but not save files
%   3. run with convflag=0 (DEFAULT) would keep the results in its own fam frame;
%       run with convflag=1 would convert results of all fams to the same frame
%       of the reference fam, which is now 043, can be changed easily.
%   4. IN DEFAULT, in order to write rst files for next-step plotaddcheck scripts, set
%       fcflag=0 and convflag=0; other options ONLY save plots to eyeball analysis
%   
%        
% Conversion algorithm:
%       (off12')_fam - (off12')_ref = (off12)_fam - (off12)_ref
%   where ' means the detection in that family, no ' means the original
%   reference point in each family, i.e. the location of that family 
%   (considered as 0,0 in each fam). Since we already have the travel time
%   difference between sta 12 & 13 in each family (4th col in rots para), 
%   what we need to do here is to convert them to a uniform frame.
%   Thus, relative to a selected fam, all detections in other fam can be
%   written in the new frame as:
%       (off12')_ref = (off12')_fam - [(off12)_fam - (off12)_ref]
%
% Features:
%   1. read detection offsets and trace files
%   2. option flag to correct filtering effect (1/0)
%   3. option flag to convert all fams to the same frame (1/0)
%   4. remove double counting, duplicates, NECESSARY
%   5. re-write those files (only when fcflag=0 and convflag=0)
%   6. re-plot the offsets variation w/ time after removing duplicates (always save plots)
%   
%
% Modified by Chao Song, chaosong@princeton.edu
% First created date:   2020/12/07
% Last modified date:   2020/12/07
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Initialization
format short e   % Set the format to 5-digit floating point
clear
% close all

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
%             '099';
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
        
famex = '099'

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
        
winlenhf = 4;

winlenlf = 16;
lofflf = 4;
ccminlf = 0.35;


%% load detections with or without a distance cutoff
SUFFIXhf = strcat('up.hf.time.',num2str(winlenhf),'_',num2str(winlenlf),'.',num2str(lofflf),...
                  '.',num2str(ccminlf));
SUFFIXlf = strcat('lf.time.',num2str(winlenhf),'_',num2str(winlenlf),'.',num2str(lofflf),...
                  '.',num2str(ccminlf));

% dcutflag = 1;
dcutflag = 0;

if dcutflag 
    hffname = strcat(rstpath, '/evtloc.ex',famex,'.dcut.',SUFFIXhf);
    lffname = strcat(rstpath, '/evtloc.ex',famex,'.dcut.',SUFFIXlf);
else
    hffname = strcat(rstpath, '/evtloc.ex',famex,'.nodcut.',SUFFIXhf);
    lffname = strcat(rstpath, '/evtloc.ex',famex,'.nodcut.',SUFFIXlf);
end

hfold = load(hffname);
lfold = load(lffname);
% 25 cols, format is:
%   E(043) N(043) E(own) N(own) dist(own) lon lat dep off12 off13 off12sec off13sec date 
%   main_arrival_time  cent_of_win  avecc ccadd1 ccadd2 ccadd3 ccadd4 ampmax ampmin eng normeng famnum
%



%% sort according to time, remove double counting
disp('Start removing double counts ...');

%%% for hf
%%% sort according to date and main arrival time, could add offset
hfsort = sortrows(hfold, [13,14]);
dateall = unique(hfsort(:,13));
hfnew = [];
for id = 1: length(dateall)
    hfday = hfsort(hfsort(:,13) == dateall(id), :);
    dtminhf = 0.5;      % min time during which only 1 detection is retained
    colnum = [14 16];
    indexthf = length(nfampool);    % the extreme case is that every fam contains a duplicate
    [hfdaynew,indsave,inddis] = RemoveDoubleCounting(hfday,dtminhf,indexthf,colnum);    
    hfnew = [hfnew; hfdaynew];    % hfdaynew is the file of all fam for one day  
    
end

%%% for lf
%%% sort according to date and main arrival time, could add offset
lfsort = sortrows(lfold, [13,14]);
dateall = unique(lfsort(:,13));
lfnew = [];
for id = 1: length(dateall)
    lfday = lfsort(lfsort(:,13) == dateall(id), :);
    dtminlf = 1.1;      % min time during which only 1 detection is retained
    colnum = [14 16];
    indextlf = length(nfampool);    % the extreme case is that every fam contains a duplicate
    [lfdaynew,indsave,inddis] = RemoveDoubleCounting(lfday,dtminlf,indextlf,colnum);    
    lfnew = [lfnew; lfdaynew];    % hfdaynew is the file of all fam for one day  
    
end


%%% save the results here
% 25 cols, format is:
%   E(043) N(043) E(own) N(own) dist(own) lon lat dep off12 off13 off12sec off13sec date 
%   main_arrival_time  cent_of_win  avecc ccadd1 ccadd2 ccadd3 ccadd4 ampmax ampmin eng normeng famnum
%
if dcutflag 
    fidhf = fopen(strcat(rstpath, '/evtloc.ex',famex,'.dcutnodou.',SUFFIXhf),'w+');
    fidlf = fopen(strcat(rstpath, '/evtloc.ex',famex,'.dcutnodou.',SUFFIXlf),'w+');
else
    fidhf = fopen(strcat(rstpath, '/evtloc.ex',famex,'.nodcutnodou.',SUFFIXhf),'w+');
    fidlf = fopen(strcat(rstpath, '/evtloc.ex',famex,'.nodcutnodou.',SUFFIXlf),'w+');
end

fprintf(fidhf,'%.4f %.4f %.4f %.4f %.4f %.4f %.4f %.4f %.4f %.4f %.4f %.4f %d %.4f %.1f %.4f %.4f %.4f %.4f %.4f %.4e %.4e %10.3e %10.3f %d \n',...
        hfnew');
fclose(fidhf);

fprintf(fidlf,'%.4f %.4f %.4f %.4f %.4f %.4f %.4f %.4f %.4f %.4f %.4f %.4f %d %.4f %.1f %.4f %.4f %.4f %.4f %.4f %.4e %.4e %10.3e %10.3f %d \n',...
        lfnew');
fclose(fidlf);









