% function loc_resolution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script is used to obtain the location resolution of the detection in
% HF and LF. In time domain, HF detection is precise to 1/4 sample in 40 sps,
% on the contrary, LF is precise to 1 sample in 20 sps. Thus the resolution in
% time is 1:8. Choose several control points, (off12, off13) = (+-1, 0) and
% (0, +-1) to obtain the resulting location could illustrate the location
% resolution.
% It seems that at least in fam 002 with PGC detector, the (-1/4,0) and (0, -1/4)
% in 40 sps are indistinguishable. But in LZB detector, it should n't have this
% issue. The reason comes from the grid size of the input slab. However, the
% slab itself is interpolated to a denser grid. Doesm't really make sense to
% use a denser interpolation.
%
%
% Chao Song, chaosong@princeton.edu
% First created date:   2020/08/13
% Last modified date:   2020/08/13
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Initialization
format short e   % Set the format to 5-digit floating point
clear
% close all

set(0,'DefaultFigureVisible','on');
% set(0,'DefaultFigureVisible','off');   % switch to show the plots or not

[scrsz, res] = pixelperinch(1);

% set path
workpath = getenv('MHOME');
pgcpath = strcat(workpath, '/Seisbasics/hypoinverse/forsummary');
lzbpath = strcat(workpath, '/Seisbasics/hypoinverse/lzbrst');

%% PGC trio
% relative arrival time to main station for each fam
contpgcoff = [
              86 20;
              ];

% +1/-1 samples in 40 sps
off1 = [1 0;
       -1 0;
       0 -1;
       0 1;
       ];
% +2/-2 samples in 40 sps, equivalent to +1/-1 samples in 20 sps, ie, LF resolution 
off2 = [2 0;
       -2 0;
       0 -2;
       0 2;
       ];   
   
%%% fc is inverted with (lf-hf_12, lf-hf_13) relative to fam 002
fcpgc = [-123.592000 48.436500 36.7900];

%%% this is inverted from (0,0) relative to fam 002, location of control points
contpgc = [-123.585000 48.436667 36.8800];

arr1 = [];
arr2 = [];
for i = 1:4
    tmp(i,1:2) = contpgcoff;
end
tmp(1:4,3:4) = off1;
arr1 = [arr1; tmp];

for i = 1:4
    tmp(i,1:2) = contpgcoff;
end
tmp(1:4,3:4) = off2;
arr2 = [arr2; tmp];
    

fid1 = fopen(strcat(pgcpath, '/offset_002_locres_hf1spl'),'w+');
fprintf(fid1,'%d %d %d %d \n',...
        arr1');
fclose(fid1);

fid2 = fopen(strcat(pgcpath, '/offset_002_locres_hf2spl'),'w+');
fprintf(fid2,'%d %d %d %d \n',...
        arr2');
fclose(fid2);


%% LZB trio
% fams used
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
%             '158';
%             '234';
            '006';
            '001';
            ];
nfam = length(nfampool);

% relative arrival time to main station for each fam
contlzboff = [
              14 -22;
              -14 -1;
              -23 7;
              2 -15;
              -25 -3;
              -26 10;
              -26 6;
              -8 -1;
              -4 -6;
              -15 4;
              -18 6;
              -32 9;
              -32 5;
              ];

% +1/-1 samples in 40 sps
off1 = [1 0;
       -1 0;
       0 -1;
       0 1;
       ];
% +2/-2 samples in 40 sps, equivalent to +1/-1 samples in 20 sps, ie, LF resolution 
off2 = [2 0;
       -2 0;
       0 -2;
       0 2;
       ];   

%%% fc is inverted with (lf-hf_12, lf-hf_13) relative to all fams
%%% the same content as /home/data2/chaosong/Seisbasics/hypoinverse/lzbrst/LZBallfamfc, which is
%%% inverted by hypoinverse from /home/data2/chaosong/Seisbasics/hypoinverse/lzbrst/chat00p.fc
fclzb = [-123.509500 48.457000 38.0200;
    -123.801000 48.561333 36.2500;
    -123.896333 48.594000 35.7800;
    -123.638833 48.474000 36.7200;
    -123.797167 48.440333 34.8100;
    -123.925000 48.599333 35.5600;
    -123.898667 48.555833 35.2500;
    -123.772167 48.575333 36.7300;
    -123.734833 48.562667 36.8900;
    -123.837500 48.587500 36.3200;
    -123.867833 48.590000 36.0100;
%     -123.984500 48.498500 34.0300;       % 158
%     -123.974167 48.474333 33.9100;       % 234
    -123.930667 48.545167 34.8600;       % 006
    -123.892500 48.467333 34.3600;       % 001
    ];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% VERSION 1 for the locations of lfe families, by inversion from hypoinverse
%%% this is inverted from (0,0) of all fams, same order, location of control points
%%% inverted from /home/data2/chaosong/Seisbasics/hypoinverse/lzbrst/chat00p.reffam
contlzb = [-123.492667 48.451500 38.1400;
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
%     -123.967000 48.458667 33.7800;       % 158
%     -123.960833 48.470000 33.9300;       % 234
    -123.908000 48.494167 34.5100;       % 006
    -123.879667 48.446167 34.2600;       % 001
    ];

arr1 = [];
arr2 = [];
for ifam = 1: nfam
    for i = 1:4
        tmp(i,1:2) = contlzboff(ifam,:);
    end
    tmp(1:4,3:4) = off1;
    arr1 = [arr1; tmp];
    
    for i = 1:4
        tmp(i,1:2) = contlzboff(ifam,:);
    end
    tmp(1:4,3:4) = off2;
    arr2 = [arr2; tmp];    
    
end

fid1 = fopen(strcat(lzbpath, '/offset_13fam_locres_hf1spl'),'w+');
fprintf(fid1,'%d %d %d %d \n',...
        arr1');
fclose(fid1);

fid2 = fopen(strcat(lzbpath, '/offset_13fam_locres_hf2spl'),'w+');
fprintf(fid2,'%d %d %d %d \n',...
        arr2');
fclose(fid2);













