% function identify_RTMs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script is used to find the best propagation direction for any possible
% periods of tremor bursts from 'identify_tremor_bursts_visually.m'
% Search over a series of trial directions, project HF directions onto this
% direciton, then apply a robust linear regression to projections vs. time,
% obtain the standard error of the slope and the Pearson correlation coefficient
% Try to find a threshold to discard some bad ones. The rest are viewed as
% RTMs that have a unidirectional propagation direction
%
%   
% NOTES:
%   2020/09/08, i am adding a new fam 006, so that i want to check if original
%               migrations are still reasonable
%   2020/09/10, i am adding a fam 001 back again for last try, 13 fam in total now
% 
% Chao Song, chaosong@princeton.edu
% First created date:   2020/12/24
% Last modified date:   2020/12/24
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
% bursts resulting from automated methods in 'identify_tremor_bursts_intertime.m' 
% with ttol=1e-3, ntol=15, no removing ;
ntran03 = [
    2003064    1.7589e+04   1.9399e+04 
    ];

ntran04 = [
    2004195   4.1773e+04   4.2636e+04
    2004195   4.2775e+04   4.3362e+04
    2004196   1.3686e+04   1.4040e+04
    2004196   3.5355e+04   3.7143e+04
    2004196   4.6730e+04   4.9360e+04
    2004196   6.9408e+04   7.1696e+04
    2004196   7.1932e+04   7.2764e+04
    2004196   7.6077e+04   7.6547e+04
    2004196   7.8718e+04   7.9172e+04
    2004196   8.1312e+04   8.1761e+04
    2004196   8.5935e+04   8.6386e+04
    2004197   5.4110e+03   5.7690e+03
    2004197   6.7100e+03   7.1720e+03
    2004197   8.3630e+03   9.2810e+03
    2004197   1.3192e+04   1.4362e+04
    2004197   1.7973e+04   2.0381e+04
    2004197   2.3784e+04   2.5900e+04
    2004197   2.6469e+04   2.7685e+04
    2004197   2.9108e+04   3.1903e+04
    2004197   3.3384e+04   3.5470e+04
    2004197   3.6589e+04   3.7168e+04
    2004197   3.7271e+04   3.7731e+04
    2004197   4.0974e+04   4.2205e+04
    2004197   4.2322e+04   4.2941e+04
    2004197   4.3031e+04   4.3632e+04
    2004197   4.4663e+04   4.5777e+04
    2004197   4.6257e+04   4.6786e+04
    2004197   5.1096e+04   5.2510e+04
    2004197   5.4453e+04   5.5867e+04
    2004197   5.6270e+04   5.7227e+04
    2004197   6.9950e+04   7.1112e+04
    2004197   7.1231e+04   7.3184e+04
    2004197   7.7283e+04   7.9811e+04
    2004197   7.9949e+04   8.1453e+04
    2004197   8.1547e+04   8.2210e+04
    2004197   8.2312e+04   8.3103e+04
    2004197   8.3402e+04   8.6659e+04
    2004198   5.9530e+03   9.7940e+03
    2004198   9.9200e+03   1.1441e+04
    2004198   1.9151e+04   2.1616e+04
    2004198   2.6626e+04   2.7261e+04
    2004198   3.2525e+04   3.5346e+04
    2004198   4.2632e+04   4.3795e+04
    2004198   5.3506e+04   5.8050e+04
    2004198   5.9632e+04   6.0437e+04
    2004198   6.2966e+04   6.3885e+04
    2004198   7.3540e+04   7.6029e+04
    2004198   8.4222e+04   8.6059e+04
    2004198   8.6148e+04   8.7248e+04
    2004199   1.3010e+03   3.1870e+03
    2004199   3.5090e+03   3.8000e+03
    2004199   4.1670e+03   6.2660e+03
    2004199   6.3670e+03   6.8580e+03
    2004199   7.0820e+03   7.5180e+03
    2004199   3.2442e+04   3.3827e+04
    2004199   3.5985e+04   3.6346e+04
    2004199   4.1721e+04   4.2114e+04
    2004199   4.2845e+04   4.3359e+04
    2004199   4.5786e+04   4.6012e+04
    2004199   4.6340e+04   4.7036e+04
    2004199   4.7227e+04   4.7633e+04
    2004199   4.7807e+04   4.9126e+04
    2004199   4.9744e+04   5.0614e+04
    2004199   5.9508e+04   6.0735e+04
    2004199   6.6269e+04   6.7368e+04
    2004199   8.0861e+04   8.3515e+04
    2004200   1.1706e+04   1.5799e+04
    2004200   1.5914e+04   1.7125e+04
    2004200   1.9104e+04   2.0110e+04
    2004200   2.8316e+04   2.9468e+04
    2004200   4.8109e+04   4.8666e+04
    2004200   4.9663e+04   5.0407e+04
    2004201   4.3700e+02   1.6290e+03
    2004201   6.9270e+03   7.7690e+03
    2004202   4.7231e+04   4.8033e+04
    2004203   1.6586e+04   1.7619e+04
    2004203   1.8155e+04   2.0776e+04
    2004203   6.1220e+04   6.1600e+04
    2004203   6.4578e+04   6.5555e+04
    2004203   6.6022e+04   6.6872e+04
    ];

ntran04mod1 = [
    2004195   4.1773e+04   4.3362e+04   % combine
    2004196   1.3686e+04   1.4040e+04
    2004196   3.5355e+04   3.7143e+04
    2004196   4.6730e+04   4.9360e+04
    2004196   6.9408e+04   7.1696e+04   % speed direct
    2004196   7.1932e+04   7.2764e+04
    2004196   7.6077e+04   7.6547e+04
    2004196   7.8718e+04   7.9172e+04
    2004196   8.1312e+04   8.1761e+04
    2004196   8.5935e+04   8.6386e+04   % speed direct
    2004197   5.4110e+03   5.7690e+03
    2004197   6.7100e+03   7.1720e+03
    2004197   8.3630e+03   9.2810e+03
    2004197   1.3192e+04   1.4362e+04   % speed direct
    2004197   1.7973e+04   1.8738e+04   % divided into 2 
    2004197   1.8756e+04   2.0381e+04   % divided into 2 
    2004197   2.3784e+04   2.5900e+04
    2004197   2.6469e+04   2.7180e+04   % modified time
    2004197   2.9772e+04   3.0618e+04   % modified time
    2004197   3.3384e+04   3.5470e+04
    2004197   3.6589e+04   3.7168e+04
    2004197   3.7271e+04   3.7731e+04
    2004197   4.0974e+04   4.2205e+04   % speed direct 90
    2004197   4.2322e+04   4.2941e+04   % speed direct 70
    2004197   4.3031e+04   4.3632e+04   % speed direct 200
    2004197   4.4663e+04   4.5777e+04
    2004197   4.6257e+04   4.6786e+04
    2004197   5.1096e+04   5.2510e+04
    2004197   5.4453e+04   5.5867e+04   % speed direct 205
    2004197   5.6270e+04   5.7227e+04   % speed direct 175
    2004197   6.9950e+04   7.1112e+04   % speed direct 250
    2004197   7.1231e+04   7.3184e+04   % speed direct 305
    2004197   7.7283e+04   7.9811e+04   % speed direct 270
    2004197   7.9949e+04   8.1453e+04   % speed direct 260
    2004197   8.1547e+04   8.2210e+04
    2004197   8.2312e+04   8.3103e+04   % speed direct 85
    2004197   8.3402e+04   8.4888e+04   % divided into 3
    2004197   8.4888e+04   8.5320e+04   % divided into 3
    2004197   8.5320e+04   8.6659e+04   % divided into 3
    2004198   5.9530e+03   9.7940e+03   % speed direct 270
    2004198   9.9200e+03   1.1441e+04
    2004198   1.9151e+04   2.1616e+04   % acceptted
    2004198   2.6626e+04   2.7261e+04
    2004198   3.2525e+04   3.3480e+04   % too complicated so modified time, but bad fitting due to spatial scatter
    2004198   4.2632e+04   4.3795e+04   % acceptted
    2004198   5.3506e+04   5.8050e+04   % acceptted
    2004198   5.9632e+04   6.0437e+04
    2004198   6.2966e+04   6.3885e+04   % acceptted
    2004198   7.3540e+04   7.6029e+04   % acceptted
    2004198   8.4222e+04   8.6400e+04   % combine
    2004199   1.3010e+03   3.1870e+03
    2004199   3.5090e+03   3.8000e+03
    2004199   4.1670e+03   6.2660e+03   % speed direct 280
    2004199   6.3670e+03   6.8580e+03
    2004199   7.0820e+03   7.5180e+03
    2004199   3.2442e+04   3.3827e+04   % speed direct 140
    2004199   3.5985e+04   3.6346e+04   % speed direct 50
    2004199   4.1721e+04   4.2114e+04
    2004199   4.2845e+04   4.3359e+04   % speed direct 270
    2004199   4.5786e+04   4.6012e+04
    2004199   4.6340e+04   4.9126e+04   % combine, care speed direct
    2004199   4.9744e+04   5.0614e+04   % speed direct 180
    2004199   5.9508e+04   6.0735e+04   % speed direct 65
    2004199   6.6269e+04   6.7368e+04   % speed direct 110
    2004199   8.0861e+04   8.2008e+04   % divided into 2
    2004199   8.2008e+04   8.3515e+04   % divided into 2
    2004200   1.2600e+04   1.5799e+04   % modified time
    2004200   1.5914e+04   1.7125e+04   % acceptted
    2004200   1.9104e+04   2.0110e+04   % check old param
    2004200   2.8316e+04   2.9468e+04   % speed direct 120
    2004200   4.8109e+04   4.8666e+04   % acceptted
    2004200   4.9663e+04   5.0407e+04
    2004201   4.3700e+02   1.6290e+03   % speed direct 225
    2004201   6.9270e+03   7.7690e+03
    2004202   4.7231e+04   4.8033e+04
    2004203   1.6586e+04   2.0776e+04   % combine them
    2004203   6.1220e+04   6.1600e+04
    2004203   6.4578e+04   6.5555e+04   % speed direct 250
    2004203   6.6022e+04   6.6872e+04   % speed direct 235
    ];

% bursts resulting from automated methods in 'identify_tremor_bursts_intertime.m' 
% with ttol=1e-3, ntol=15, no removing ;
ntran05 = [
    2005254   5.2247e+04   5.2795e+04
    2005254   5.6369e+04   5.6994e+04
    2005254   7.0860e+04   7.1612e+04
    2005255   1.8219e+04   1.9459e+04
    2005255   3.4400e+04   3.5612e+04
    2005255   4.6099e+04   4.6949e+04
    2005255   5.1129e+04   5.1945e+04
    2005255   5.3329e+04   5.3807e+04
    2005255   5.8050e+04   5.9394e+04
    2005255   6.0381e+04   6.1557e+04
    2005255   6.1678e+04   6.2487e+04
    2005255   6.7268e+04   6.8535e+04
    2005255   7.1472e+04   7.2019e+04
    2005255   7.3862e+04   7.4339e+04
    2005255   7.4790e+04   7.5708e+04
    2005255   8.4296e+04   8.4871e+04
    2005255   8.5130e+04   8.5634e+04
    2005256   3.6200e+03   4.9570e+03
    2005256   8.0330e+03   1.0404e+04
    2005256   1.2631e+04   1.3269e+04
    2005256   1.3412e+04   1.4346e+04
    2005256   1.4465e+04   1.4962e+04
    2005256   1.5237e+04   1.6010e+04
    2005256   1.6237e+04   1.6733e+04
    2005256   1.9069e+04   2.0203e+04
    2005256   2.0437e+04   2.2294e+04
    2005256   2.7719e+04   2.9172e+04
    2005256   2.9715e+04   3.1203e+04
    2005256   3.2636e+04   3.4444e+04
    2005256   3.4604e+04   3.5238e+04
    2005256   3.5597e+04   4.1112e+04
    2005256   4.1204e+04   4.5189e+04
    2005256   4.5316e+04   4.8421e+04
    2005256   4.8848e+04   5.0185e+04
    2005256   5.1394e+04   5.4053e+04
    2005256   5.4816e+04   5.7088e+04
    2005256   6.0723e+04   6.4865e+04
    2005256   7.1060e+04   7.3914e+04
    2005256   7.6381e+04   7.8664e+04
    2005256   7.8780e+04   8.0534e+04
    2005256   8.1490e+04   8.2062e+04
    2005256   8.4804e+04   8.5328e+04
    2005257   4.9440e+03   7.7460e+03
    2005257   9.0080e+03   1.3853e+04
    2005257   1.5973e+04   1.8655e+04
    2005257   2.1070e+04   2.5365e+04
    2005257   2.5495e+04   2.8055e+04
    2005257   3.5350e+04   3.7598e+04
    2005257   3.8513e+04   4.5418e+04
    2005257   4.6153e+04   4.8618e+04
    2005257   6.1449e+04   6.4012e+04
    2005257   6.9140e+04   6.9892e+04
    2005257   7.3440e+04   8.2043e+04
    2005257   8.2159e+04   8.3380e+04
    2005257   8.3514e+04   8.4440e+04
    2005258   2.4602e+04   2.8962e+04
    2005258   3.5000e+04   4.0021e+04
    2005258   7.2884e+04   7.4202e+04
    2005259   1.5400e+02   1.6340e+03
    2005259   1.9560e+03   4.0460e+03
    2005259   5.2320e+03   7.7140e+03
    2005259   3.7000e+04   4.1293e+04
    2005259   7.4255e+04   7.4723e+04
    2005260   1.4860e+03   2.4810e+03
    2005260   2.6510e+03   3.4380e+03
    2005260   3.6040e+03   6.6780e+03
    2005260   6.7690e+03   7.0310e+03
    2005260   7.1270e+03   7.4700e+03
    2005260   7.6740e+03   8.1000e+03
    2005260   9.4090e+03   1.0357e+04
    2005260   4.0871e+04   4.1640e+04
    2005260   5.1893e+04   5.2367e+04
    2005260   5.6422e+04   5.7180e+04
    2005261   8.5290e+03   9.6210e+03
    2005261   9.7470e+03   1.0313e+04
    2005261   1.0683e+04   1.1169e+04
    2005261   1.1318e+04   1.2103e+04
    ];

ntran05mod1 = [
    2005254   5.2247e+04   5.2795e+04
    2005254   5.6369e+04   5.6994e+04
    2005254   7.0860e+04   7.1612e+04
    2005255   1.8219e+04   1.9459e+04
    2005255   3.4400e+04   3.5612e+04   % speed direct 80
    2005255   4.6099e+04   4.6949e+04   % speed direct 65
    2005255   5.1129e+04   5.1945e+04   % speed direct 250
    2005255   5.3329e+04   5.3807e+04
    2005255   5.8050e+04   5.9394e+04   % speed direct 95
    2005255   6.0381e+04   6.1557e+04   % speed direct 150
    2005255   6.1678e+04   6.2487e+04
    2005255   6.7268e+04   6.8535e+04   % speed direct 80
    2005255   7.1472e+04   7.2019e+04
    2005255   7.3862e+04   7.4339e+04
    2005255   7.5060e+04   7.5708e+04   % care speed direc
    2005255   8.4296e+04   8.5634e+04   % combine them
    2005256   3.6200e+03   4.9570e+03   % speed direct 265
    2005256   8.0330e+03   9.8280e+03   % change end time
    2005256   1.2631e+04   1.3269e+04
    2005256   1.3412e+04   1.4346e+04   % speed direct 65   
    2005256   1.4465e+04   1.4962e+04
    2005256   1.5237e+04   1.6010e+04   % speed direct 80
    2005256   1.6237e+04   1.6733e+04
    2005256   1.9069e+04   2.0203e+04   % speed direct 195
    2005256   2.1492e+04   2.2294e+04   % change start time
    2005256   2.8080e+04   2.9172e+04   % change start time
    2005256   2.9715e+04   3.1203e+04   % speed direct 70
    2005256   3.2636e+04   3.4444e+04   % speed direct 105
    2005256   3.4776e+04   3.5238e+04   % change start time
    2005256   3.6900e+04   3.8520e+04   % change times
    2005256   4.1204e+04   4.5189e+04
    2005256   4.6440e+04   4.8421e+04   % change start time
    2005256   4.8848e+04   5.0185e+04
    2005256   5.2020e+04   5.2596e+04   % change times
    2005256   5.4816e+04   5.7088e+04
    2005256   6.2640e+04   6.4865e+04   % change start time
    2005256   7.1496e+04   7.2720e+04   % divided into 2,
    2005256   7.2720e+04   7.3914e+04   % divided into 2,   
    2005256   7.6381e+04   7.7940e+04   % change end time
    2005256   7.8780e+04   8.0534e+04
    2005256   8.1490e+04   8.2062e+04
    2005256   8.4804e+04   8.5328e+04   % both direc bad, try 245
    2005257   6.1200e+03   7.7460e+03   % change start time
    2005257   9.0080e+03   1.2240e+04   % divided into 2,
    2005257   1.2240e+04   1.3853e+04   % divided into 2,
    2005257   1.6920e+04   1.8655e+04   % change start time
    2005257   2.1070e+04   2.5365e+04   % accepted
    2005257   2.6280e+04   2.8055e+04   % change start time
    2005257   3.6000e+04   3.7598e+04   % change start time
    2005257   3.9000e+04   4.5418e+04   % change start time
    2005257   4.6153e+04   4.8618e+04
    2005257   6.1449e+04   6.4012e+04   % accepted
    2005257   6.9140e+04   6.9892e+04
    2005257   7.3440e+04   7.8840e+04   % change end time
    2005257   8.2159e+04   8.3380e+04   % speed direct 350
    2005257   8.3514e+04   8.4440e+04
    2005258   2.6640e+04   2.8962e+04   % change start time
    2005258   3.5000e+04   3.6540e+04   % divided into 3,
    2005258   3.6540e+04   3.8160e+04   % divided into 3,
    2005258   3.8160e+04   4.0021e+04   % divided into 3,
    2005258   7.2884e+04   7.4202e+04
    2005259   1.5400e+02   1.6340e+03   % speed direct 75
    2005259   1.9560e+03   4.0460e+03   % speed direct 185
    2005259   5.2320e+03   7.7140e+03   % check OLD param. to determine
    2005259   3.7000e+04   4.1293e+04   % speed direct 85
    2005259   7.4255e+04   7.4723e+04
    2005260   1.4860e+03   2.4810e+03
    2005260   2.6510e+03   3.4380e+03   % speed direct 80
    2005260   3.6040e+03   4.5000e+03   % divided into 2,
    2005260   4.6080e+03   6.6780e+03   % divided into 2,
    2005260   6.7690e+03   7.0310e+03
    2005260   7.1270e+03   7.4700e+03
    2005260   7.6740e+03   8.1000e+03
    2005260   9.4090e+03   1.0357e+04   % speed direct 55
    2005260   4.0871e+04   4.1640e+04
    2005260   5.1893e+04   5.2367e+04
    2005260   5.6300e+04   5.8200e+04   % use old times
    2005261   8.5290e+03   9.6210e+03
    2005261   9.7470e+03   1.0313e+04
    2005261   1.0683e+04   1.1169e+04
    2005261   1.1318e+04   1.2103e+04
    ];

% Combined accepted ones from all ETS
ntranacc = [
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

test = [
%     2005257   7.3440e+04   7.5460e+04
%     2005257   7.5810e+04   7.6790e+04
%     2005257   7.6860e+04   7.8840e+04 
%     2005257   7.3440e+04   7.8840e+04   % change end time, y, but check speed direc, y, but use rmse 235
%     2004198   8.4222e+04   8.6400e+04
    2005256   3.6900e+04   3.8520e+04   % change times, y
    ];

trange = test;

xran = [-15 25];
yran = [-20 20];


%%
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


%%
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

% a box around fam 002 region
reg002 = [ -5 -15;
            25 0
            35 -20;
            5 -35;
           -5 -15];
       
% index of migrations whose majority is NOT at 002 region but have some detections there  
indmig002 = [2,4,5,6,7,8,9,12,13,14,17,19,20,24,25,28,29,31,32,35,36,37,39,41,42,43,44,46,47,48,49,...
             50,52];

for i = 1: size(trange,1)
% for i = 1: 211
%     i=5;
    disp(trange(i,:));
    indhf = find(hftime(:,13)==trange(i,1) & hftime(:,15)>=trange(i,2) & ...
                 hftime(:,15)<=trange(i,3));
    mighf = hftime(indhf,:);
    
    indlf = find(lftime(:,13)==trange(i,1) & lftime(:,15)>=trange(i,2) & ...
                 lftime(:,15)<=trange(i,3));
    miglf = lftime(indlf,:); 
    
    %%% Actually used is the pre-determined best prop direc to do the fitting
    mighfdum = mighf;
    miglfdum = miglf;
    
%     % for the migrations whose majority is occuring not at 002 region, leave out 
%     % the detections at 002 region
%     if ~isempty(intersect(indmig002,i))
%         % only use the detections ouside the bound
%         [is,ion] = inpolygon(mighfdum(:,1),mighfdum(:,2),reg002(:,1),reg002(:,2));
%         isinreg = is | ion;
%         mighfdum = mighfdum(isinreg ~= 1, :);
%         
%         [is,ion] = inpolygon(miglfdum(:,1),miglfdum(:,2),reg002(:,1),reg002(:,2));
%         isinreg = is | ion;
%         miglfdum = miglfdum(isinreg ~= 1, :);
%     end
    
    for j = 1: size(mighfdum,1)
        x0 = mighfdum(j,1);
        y0 = mighfdum(j,2);
        [newx,newy] = coordinate_rot(x0,y0,-(angbest1(i)-90),0,0);
        mighfdum(j,1) = newx;
        mighfdum(j,2) = newy;
    end

    for j = 1: size(miglfdum,1)
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
    text(f.ax(1),0.04,0.1,'LZB','FontSize',11,'unit','normalized','EdgeColor','k','Margin',2);
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
    text(f.ax(2),0.04,0.1,'LZB','FontSize',11,'unit','normalized','EdgeColor','k','Margin',2);
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
    %%% Actually used is the pre-determined best prop direc to do the fitting
    mighfdum = mighf;
    miglfdum = miglf;
    
%     % for the migrations whose majority is occuring not at 002 region, leave out 
%     % the detections at 002 region
%     if ~isempty(intersect(indmig002,i))
%         % only use the detections ouside the bound
%         [is,ion] = inpolygon(mighfdum(:,1),mighfdum(:,2),reg002(:,1),reg002(:,2));
%         isinreg = is | ion;
%         mighfdum = mighfdum(isinreg ~= 1, :);
%         
%         [is,ion] = inpolygon(miglfdum(:,1),miglfdum(:,2),reg002(:,1),reg002(:,2));
%         isinreg = is | ion;
%         miglfdum = miglfdum(isinreg ~= 1, :);
%     end
    
    for j = 1: size(mighfdum,1)
        x0 = mighfdum(j,1);
        y0 = mighfdum(j,2);
        [newx,newy] = coordinate_rot(x0,y0,-(angbest2(i)-90),0,0);
        mighfdum(j,1) = newx;
        mighfdum(j,2) = newy;
    end

    for j = 1: size(miglfdum,1)
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
    
%     % subplot 4 of figure i
%     hold(f.ax(4),'on');
%     f.ax(4).FontSize = 9;
%     numhf(i) = length(mighfdum2(:,1));
%     numlf(i) = length(miglfdum2(:,1));
%     resprophf(i,1:numhf(i)) = outphfprop.residuals;
%     resproplf(i,1:numlf(i)) = outplfprop.residuals+(intcptproplf(i)-intcptprophf(i));
%     [xlochf, a, b, pdfvalhf] = weighted_bar_pdf(resprophf(i,1:numhf(i)), wthf, 0.5, 'int');
%     hfHdl = bar(f.ax(4),xlochf, pdfvalhf,1,'stacked');
%     hfHdl(1).FaceColor = [0.6 0.6 0.6];
%     [xloclf, ~, ~, pdfvallf] = weighted_bar_pdf(resproplf(i,1:numlf(i)), wtlf, 0.5, 'int');
%     lfHdl = bar(f.ax(4),xloclf, pdfvallf,1,'stacked','facea',0.6);
%     lfHdl(1).FaceColor = [0.6 1 1];
%     text(f.ax(4),0.04,0.91,'d','FontSize',11,'unit','normalized','EdgeColor','k','Margin',2);
%     text(f.ax(4),0.4,0.92,sprintf('LF-HF = %.2f km',offset(i)),'FontSize',11,'unit','normalized');
%     f.ax(4).Box = 'on';
%     grid(f.ax(4), 'on');
%     f.ax(4).GridLineStyle = '--';
%     ymax = f.ax(4).YLim(2)+0.1;
%     ylim(f.ax(4),[0 ymax]);
% %     ylim(f.ax(4),[0 1]);
%     xlabel(f.ax(4),'Residual in prop. (km)','fontsize',11);
%     ylabel(f.ax(4),'PDF estimate','fontsize',11);
%     hold(f.ax(4),'off');

%     if pearwthf(i) >= 0.5 && sehf(i,1)<3
%         print(f.fig,'-dpdf',strcat('/home/data2/chaosong/CurrentResearch/Song_Rubin_2020/Figures/0105nrtm05_',...
%             num2str(i),'.pdf'));
%     end

%     if ~isempty(intersect(largediffind,i))
%         print(f.fig,'-dpdf',strcat(rstpath,'/bigoffdiff_',num2str(i),'.pdf'));
%     end


end

[tmp,tmpind] = sort(sehf(:,1),'descend');
sortsehf = [tmpind tmp];
[tmp,tmpind] = sort(pearwthf(:),'ascend');
sortpearwthf = [tmpind tmp];
[tmp,tmpind] = sort(self(:,2),'descend');
sortself = [tmpind tmp];
[tmp,tmpind] = sort(sehf2(:,1),'descend');
sortsehf2 = [tmpind tmp];
[tmp,tmpind] = sort(pearwthf2(:),'ascend');
sortpearwthf2 = [tmpind tmp];
[tmp,tmpind] = sort(self2(:,2),'descend');
sortself2 = [tmpind tmp];

[tmp,tmpind] = sort(sehfbest(:,1),'descend');
sortsehfbest = [tmpind tmp];
[tmp,tmpind] = sort(selfbest(:,2),'descend');
sortselfbest = [tmpind tmp];



%%
tmp1 = sortpearwthf(sortpearwthf(:,2)>=0.5, 1);
tmp2 = sortsehf(sortsehf(:,2)<3, 1);
tmp3 = sortself(sortself(:,2)<3, 1);

tmp21 = sortpearwthf2(sortpearwthf2(:,2)>=0.5, 1);
tmp22 = sortsehf2(sortsehf2(:,2)<3, 1);
tmp23 = sortself2(sortself2(:,2)<3, 1);


if issame(trange,ntran04mod1)
    ind1 = intersect(tmp1,tmp2);
    ind2 = intersect(tmp21,tmp22);
%     ind04acc = union(ind1,ind2);    
%     ind04dis = setdiff(1:length(trange),union(ind1,ind2))';    % neither direction works, throw way
%     ind04spd = setdiff(ind2,ind1)';
%     ind04rmse = setdiff(ind1,ind2)';
    ind04acc = intersect(ind1,ind2)';
    ind04dis = setdiff(1:length(trange),ind04acc)';
elseif issame(trange,ntran05mod1)
    ind1 = intersect(tmp1,tmp2);
    ind2 = intersect(tmp21,tmp22);
%     ind05acc = union(ind1,ind2);    
%     ind05dis = setdiff(1:length(trange),union(ind1,ind2))';    % neither direction works, throw way
%     ind04spd = setdiff(ind2,ind1)';
%     ind04rmse = setdiff(ind1,ind2)';
    ind05acc = intersect(ind1,ind2)';
    ind05dis = setdiff(1:length(trange),ind05acc)';
else
    ind03 = intersect(tmp1,tmp2);
    acctran03 = ntran03(ind03,:);
end


%% hardcopy of the migrations that need further check
% copy from ntran04, discarded are commented out
ntran04raw1 = [
%     2004195   4.1773e+04   4.2636e+04
    2004195   4.1773e+04   4.3362e+04   % combine
%     2004196   1.3686e+04   1.4040e+04
%     2004196   3.5355e+04   3.7143e+04
%     2004196   4.6730e+04   4.9360e+04
    2004196   6.9408e+04   7.1696e+04   % speed direct
%     2004196   7.1932e+04   7.2764e+04
%     2004196   7.6077e+04   7.6547e+04
%     2004196   7.8718e+04   7.9172e+04
%     2004196   8.1312e+04   8.1761e+04
    2004196   8.5935e+04   8.6386e+04   % speed direct
%     2004197   5.4110e+03   5.7690e+03
%     2004197   6.7100e+03   7.1720e+03
%     2004197   8.3630e+03   9.2810e+03
    2004197   1.3192e+04   1.4362e+04   % speed direct
    2004197   1.7973e+04   1.8738e+04   1.8756e+04   2.0381e+04  % divided into 2
%     2004197   2.3784e+04   2.5900e+04
    2004197   2.6469e+04   2.7180e+04   % modified time
    2004197   2.9772e+04   3.0618e+04   % modified time
%     2004197   3.3384e+04   3.5470e+04
%     2004197   3.6589e+04   3.7168e+04
%     2004197   3.7271e+04   3.7731e+04
    2004197   4.0974e+04   4.2205e+04   % speed direct 90
    2004197   4.2322e+04   4.2941e+04   % speed direct 70
    2004197   4.3031e+04   4.3632e+04   % speed direct 200
%     2004197   4.4663e+04   4.5777e+04
%     2004197   4.6257e+04   4.6786e+04
%     2004197   5.1096e+04   5.2510e+04
    2004197   5.4453e+04   5.5867e+04   % speed direct 205
    2004197   5.6270e+04   5.7227e+04   % speed direct 175
    2004197   6.9950e+04   7.1112e+04   % speed direct 250
    2004197   7.1231e+04   7.3184e+04   % speed direct 305
    2004197   7.7283e+04   7.9811e+04   % speed direct 270
    2004197   7.9949e+04   8.1453e+04   % speed direct 260
%     2004197   8.1547e+04   8.2210e+04
    2004197   8.2312e+04   8.3103e+04   % speed direct 85
    2004197   8.3402e+04   8.4888e+04   8.5320e+04   8.6659e+04   % divided into 3, 1 discarded
    2004198   5.9530e+03   9.7940e+03   % speed direct 270
%     2004198   9.9200e+03   1.1441e+04
    2004198   1.9151e+04   2.1616e+04   % acceptted
%     2004198   2.6626e+04   2.7261e+04
    2004198   3.2525e+04   3.3480e+04   % too complicated so modified time, but bad fitting due to spatial scatter
    2004198   4.2632e+04   4.3795e+04   % acceptted
    2004198   5.3506e+04   5.8050e+04   % acceptted
%     2004198   5.9632e+04   6.0437e+04
    2004198   6.2966e+04   6.3885e+04   % acceptted
    2004198   7.3540e+04   7.6029e+04   % acceptted
%     2004198   8.4222e+04   8.6059e+04
    2004198   8.4222e+04   8.7248e+04   % combine
%     2004199   1.3010e+03   3.1870e+03
%     2004199   3.5090e+03   3.8000e+03
    2004199   4.1670e+03   6.2660e+03   % speed direct 280
%     2004199   6.3670e+03   6.8580e+03
%     2004199   7.0820e+03   7.5180e+03
    2004199   3.2442e+04   3.3827e+04   % speed direct 140
    2004199   3.5985e+04   3.6346e+04   % speed direct 50
%     2004199   4.1721e+04   4.2114e+04
    2004199   4.2845e+04   4.3359e+04   % speed direct 270
%     2004199   4.5786e+04   4.6012e+04
%     2004199   4.6340e+04   4.7036e+04
%     2004199   4.7227e+04   4.7633e+04
    2004199   4.6340e+04   4.9126e+04   % combine, care speed direct
    2004199   4.9744e+04   5.0614e+04   % speed direct 180
    2004199   5.9508e+04   6.0735e+04   % speed direct 65
    2004199   6.6269e+04   6.7368e+04   % speed direct 110
    2004199   8.0861e+04   8.2008e+04   8.2008e+04   8.3515e+04  % divided into 2
    2004200   1.2600e+04   1.5799e+04   % modified time
    2004200   1.5914e+04   1.7125e+04   % acceptted
    2004200   1.9104e+04   2.0110e+04   % check old param
    2004200   2.8316e+04   2.9468e+04   % speed direct 120
    2004200   4.8109e+04   4.8666e+04   % acceptted
%     2004200   4.9663e+04   5.0407e+04
    2004201   4.3700e+02   1.6290e+03   % speed direct 225
%     2004201   6.9270e+03   7.7690e+03
%     2004202   4.7231e+04   4.8033e+04
%     2004203   1.6586e+04   1.7619e+04
    2004203   1.6586e+04   2.0776e+04   % combine them
%     2004203   6.1220e+04   6.1600e+04
    2004203   6.4578e+04   6.5555e+04   % speed direct 250
    2004203   6.6022e+04   6.6872e+04   % speed direct 235
    ];

ntran04mod1 = [
    2004195   4.1773e+04   4.3362e+04   % combine
    2004196   1.3686e+04   1.4040e+04
    2004196   3.5355e+04   3.7143e+04
    2004196   4.6730e+04   4.9360e+04
    2004196   6.9408e+04   7.1696e+04   % speed direct
    2004196   7.1932e+04   7.2764e+04
    2004196   7.6077e+04   7.6547e+04
    2004196   7.8718e+04   7.9172e+04
    2004196   8.1312e+04   8.1761e+04
    2004196   8.5935e+04   8.6386e+04   % speed direct
    2004197   5.4110e+03   5.7690e+03
    2004197   6.7100e+03   7.1720e+03
    2004197   8.3630e+03   9.2810e+03
    2004197   1.3192e+04   1.4362e+04   % speed direct
    2004197   1.7973e+04   1.8738e+04   % divided into 2 
    2004197   1.8756e+04   2.0381e+04   % divided into 2 
    2004197   2.3784e+04   2.5900e+04
    2004197   2.6469e+04   2.7180e+04   % modified time
    2004197   2.9772e+04   3.0618e+04   % modified time
    2004197   3.3384e+04   3.5470e+04
    2004197   3.6589e+04   3.7168e+04
    2004197   3.7271e+04   3.7731e+04
    2004197   4.0974e+04   4.2205e+04   % speed direct 90
    2004197   4.2322e+04   4.2941e+04   % speed direct 70
    2004197   4.3031e+04   4.3632e+04   % speed direct 200
    2004197   4.4663e+04   4.5777e+04
    2004197   4.6257e+04   4.6786e+04
    2004197   5.1096e+04   5.2510e+04
    2004197   5.4453e+04   5.5867e+04   % speed direct 205
    2004197   5.6270e+04   5.7227e+04   % speed direct 175
    2004197   6.9950e+04   7.1112e+04   % speed direct 250
    2004197   7.1231e+04   7.3184e+04   % speed direct 305
    2004197   7.7283e+04   7.9811e+04   % speed direct 270
    2004197   7.9949e+04   8.1453e+04   % speed direct 260
    2004197   8.1547e+04   8.2210e+04
    2004197   8.2312e+04   8.3103e+04   % speed direct 85
    2004197   8.3402e+04   8.4888e+04   % divided into 3
    2004197   8.4888e+04   8.5320e+04   % divided into 3
    2004197   8.5320e+04   8.6659e+04   % divided into 3
    2004198   5.9530e+03   9.7940e+03   % speed direct 270
    2004198   9.9200e+03   1.1441e+04
    2004198   1.9151e+04   2.1616e+04   % acceptted
    2004198   2.6626e+04   2.7261e+04
    2004198   3.2525e+04   3.3480e+04   % too complicated so modified time, but bad fitting due to spatial scatter
    2004198   4.2632e+04   4.3795e+04   % acceptted
    2004198   5.3506e+04   5.8050e+04   % acceptted
    2004198   5.9632e+04   6.0437e+04
    2004198   6.2966e+04   6.3885e+04   % acceptted
    2004198   7.3540e+04   7.6029e+04   % acceptted
    2004198   8.4222e+04   8.6400e+04   % combine
    2004199   1.3010e+03   3.1870e+03
    2004199   3.5090e+03   3.8000e+03
    2004199   4.1670e+03   6.2660e+03   % speed direct 280
    2004199   6.3670e+03   6.8580e+03
    2004199   7.0820e+03   7.5180e+03
    2004199   3.2442e+04   3.3827e+04   % speed direct 140
    2004199   3.5985e+04   3.6346e+04   % speed direct 50
    2004199   4.1721e+04   4.2114e+04
    2004199   4.2845e+04   4.3359e+04   % speed direct 270
    2004199   4.5786e+04   4.6012e+04
    2004199   4.6340e+04   4.9126e+04   % combine, care speed direct
    2004199   4.9744e+04   5.0614e+04   % speed direct 180
    2004199   5.9508e+04   6.0735e+04   % speed direct 65
    2004199   6.6269e+04   6.7368e+04   % speed direct 110
    2004199   8.0861e+04   8.2008e+04   % divided into 2
    2004199   8.2008e+04   8.3515e+04   % divided into 2
    2004200   1.2600e+04   1.5799e+04   % modified time
    2004200   1.5914e+04   1.7125e+04   % acceptted
    2004200   1.9104e+04   2.0110e+04   % check old param
    2004200   2.8316e+04   2.9468e+04   % speed direct 120
    2004200   4.8109e+04   4.8666e+04   % acceptted
    2004200   4.9663e+04   5.0407e+04
    2004201   4.3700e+02   1.6290e+03   % speed direct 225
    2004201   6.9270e+03   7.7690e+03
    2004202   4.7231e+04   4.8033e+04
    2004203   1.6586e+04   2.0776e+04   % combine them
    2004203   6.1220e+04   6.1600e+04
    2004203   6.4578e+04   6.5555e+04   % speed direct 250
    2004203   6.6022e+04   6.6872e+04   % speed direct 235
    ];

% finally accepted ones combined by uncommented ones in ntran04redo and ntran04speed
ntran04acc = [
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
    ];

%%
% copy from ntran05, discarded are commented out
ntran05raw1 = [
%     2005254   5.2247e+04   5.2795e+04
%     2005254   5.6369e+04   5.6994e+04
%     2005254   7.0860e+04   7.1612e+04
%     2005255   1.8219e+04   1.9459e+04
    2005255   3.4400e+04   3.5612e+04   % speed direct 80
    2005255   4.6099e+04   4.6949e+04   % speed direct 65
    2005255   5.1129e+04   5.1945e+04   % speed direct 250
%     2005255   5.3329e+04   5.3807e+04
    2005255   5.8050e+04   5.9394e+04   % speed direct 95
    2005255   6.0381e+04   6.1557e+04   % speed direct 150
%     2005255   6.1678e+04   6.2487e+04
    2005255   6.7268e+04   6.8535e+04   % speed direct 80
%     2005255   7.1472e+04   7.2019e+04
%     2005255   7.3862e+04   7.4339e+04
    2005255   7.5060e+04   7.5708e+04   % care speed direc
%     2005255   8.4296e+04   8.4871e+04
    2005255   8.4296e+04   8.5634e+04   % combine them
    2005256   3.6200e+03   4.9570e+03   % speed direct 265
    2005256   8.0330e+03   9.8280e+03   % change end time
%     2005256   1.2631e+04   1.3269e+04
    2005256   1.3412e+04   1.4346e+04   % speed direct 65   
%     2005256   1.4465e+04   1.4962e+04
    2005256   1.5237e+04   1.6010e+04   % speed direct 80
%     2005256   1.6237e+04   1.6733e+04
    2005256   1.9069e+04   2.0203e+04   % speed direct 195
    2005256   2.1492e+04   2.2294e+04   % change start time
    2005256   2.8080e+04   2.9172e+04   % change start time
    2005256   2.9715e+04   3.1203e+04   % speed direct 70
    2005256   3.2636e+04   3.4444e+04   % speed direct 105
    2005256   3.4776e+04   3.5238e+04   % change start time
    2005256   3.6900e+04   3.8520e+04   % change times
%     2005256   4.1204e+04   4.5189e+04
    2005256   4.6440e+04   4.8421e+04   % change start time
%     2005256   4.8848e+04   5.0185e+04
    2005256   5.2020e+04   5.2596e+04   % change times
%     2005256   5.4816e+04   5.7088e+04
    2005256   6.2640e+04   6.4865e+04   % change start time
    2005256   7.1496e+04   7.2720e+04   7.2720e+04  7.3914e+04   % divided into 2,  
    2005256   7.6381e+04   7.7940e+04   % change end time
%     2005256   7.8780e+04   8.0534e+04
%     2005256   8.1490e+04   8.2062e+04
    2005256   8.4804e+04   8.5328e+04   % both direc bad, try 245
    2005257   6.1200e+03   7.7460e+03   % change start time
    2005257   9.0080e+03   1.2240e+04   1.2240e+04   1.3853e+04   % divided into 2,
    2005257   1.6920e+04   1.8655e+04   % change start time
    2005257   2.1070e+04   2.5365e+04   % accepted
    2005257   2.6280e+04   2.8055e+04   % change start time
    2005257   3.6000e+04   3.7598e+04   % change start time
    2005257   3.9000e+04   4.5418e+04   % change start time
%     2005257   4.6153e+04   4.8618e+04
    2005257   6.1449e+04   6.4012e+04   % accepted
%     2005257   6.9140e+04   6.9892e+04
    2005257   7.3440e+04   7.8840e+04   % change end time
    2005257   8.2159e+04   8.3380e+04   % speed direct 350
%     2005257   8.3514e+04   8.4440e+04
    2005258   2.6640e+04   2.8962e+04   % change start time
    2005258   3.5000e+04   3.6540e+04   3.6540e+04   3.8160e+04   3.8160e+04   4.0021e+04   % divided into 3
%     2005258   7.2884e+04   7.4202e+04
    2005259   1.5400e+02   1.6340e+03   % speed direct 75
    2005259   1.9560e+03   4.0460e+03   % speed direct 185
    2005259   5.2320e+03   7.7140e+03   % check OLD param. to determine
    2005259   3.7000e+04   4.1293e+04   % speed direct 85
%     2005259   7.4255e+04   7.4723e+04
%     2005260   1.4860e+03   2.4810e+03
    2005260   2.6510e+03   3.4380e+03   % speed direct 80
    2005260   3.6040e+03   4.5000e+03   4.6080e+03  6.6780e+03   % divided into 2,
%     2005260   6.7690e+03   7.0310e+03
%     2005260   7.1270e+03   7.4700e+03
%     2005260   7.6740e+03   8.1000e+03
    2005260   9.4090e+03   1.0357e+04   % speed direct 55
%     2005260   4.0871e+04   4.1640e+04
%     2005260   5.1893e+04   5.2367e+04
    2005260   5.6300e+04   5.8200e+04   % use old times
%     2005261   8.5290e+03   9.6210e+03
%     2005261   9.7470e+03   1.0313e+04
%     2005261   1.0683e+04   1.1169e+04
%     2005261   1.1318e+04   1.2103e+04
    ];

ntran05mod1 = [
    2005254   5.2247e+04   5.2795e+04
    2005254   5.6369e+04   5.6994e+04
    2005254   7.0860e+04   7.1612e+04
    2005255   1.8219e+04   1.9459e+04
    2005255   3.4400e+04   3.5612e+04   % speed direct 80
    2005255   4.6099e+04   4.6949e+04   % speed direct 65
    2005255   5.1129e+04   5.1945e+04   % speed direct 250
    2005255   5.3329e+04   5.3807e+04
    2005255   5.8050e+04   5.9394e+04   % speed direct 95
    2005255   6.0381e+04   6.1557e+04   % speed direct 150
    2005255   6.1678e+04   6.2487e+04
    2005255   6.7268e+04   6.8535e+04   % speed direct 80
    2005255   7.1472e+04   7.2019e+04
    2005255   7.3862e+04   7.4339e+04
    2005255   7.5060e+04   7.5708e+04   % care speed direc
    2005255   8.4296e+04   8.5634e+04   % combine them
    2005256   3.6200e+03   4.9570e+03   % speed direct 265
    2005256   8.0330e+03   9.8280e+03   % change end time
    2005256   1.2631e+04   1.3269e+04
    2005256   1.3412e+04   1.4346e+04   % speed direct 65   
    2005256   1.4465e+04   1.4962e+04
    2005256   1.5237e+04   1.6010e+04   % speed direct 80
    2005256   1.6237e+04   1.6733e+04
    2005256   1.9069e+04   2.0203e+04   % speed direct 195
    2005256   2.1492e+04   2.2294e+04   % change start time
    2005256   2.8080e+04   2.9172e+04   % change start time
    2005256   2.9715e+04   3.1203e+04   % speed direct 70
    2005256   3.2636e+04   3.4444e+04   % speed direct 105
    2005256   3.4776e+04   3.5238e+04   % change start time
    2005256   3.6900e+04   3.8520e+04   % change times
    2005256   4.1204e+04   4.5189e+04
    2005256   4.6440e+04   4.8421e+04   % change start time
    2005256   4.8848e+04   5.0185e+04
    2005256   5.2020e+04   5.2596e+04   % change times
    2005256   5.4816e+04   5.7088e+04
    2005256   6.2640e+04   6.4865e+04   % change start time
    2005256   7.1496e+04   7.2720e+04   % divided into 2,
    2005256   7.2720e+04   7.3914e+04   % divided into 2,   
    2005256   7.6381e+04   7.7940e+04   % change end time
    2005256   7.8780e+04   8.0534e+04
    2005256   8.1490e+04   8.2062e+04
    2005256   8.4804e+04   8.5328e+04   % both direc bad, try 245
    2005257   6.1200e+03   7.7460e+03   % change start time
    2005257   9.0080e+03   1.2240e+04   % divided into 2,
    2005257   1.2240e+04   1.3853e+04   % divided into 2,
    2005257   1.6920e+04   1.8655e+04   % change start time
    2005257   2.1070e+04   2.5365e+04   % accepted
    2005257   2.6280e+04   2.8055e+04   % change start time
    2005257   3.6000e+04   3.7598e+04   % change start time
    2005257   3.9000e+04   4.5418e+04   % change start time
    2005257   4.6153e+04   4.8618e+04
    2005257   6.1449e+04   6.4012e+04   % accepted
    2005257   6.9140e+04   6.9892e+04
    2005257   7.3440e+04   7.8840e+04   % change end time
    2005257   8.2159e+04   8.3380e+04   % speed direct 350
    2005257   8.3514e+04   8.4440e+04
    2005258   2.6640e+04   2.8962e+04   % change start time
    2005258   3.5000e+04   3.6540e+04   % divided into 3,
    2005258   3.6540e+04   3.8160e+04   % divided into 3,
    2005258   3.8160e+04   4.0021e+04   % divided into 3,
    2005258   7.2884e+04   7.4202e+04
    2005259   1.5400e+02   1.6340e+03   % speed direct 75
    2005259   1.9560e+03   4.0460e+03   % speed direct 185
    2005259   5.2320e+03   7.7140e+03   % check OLD param. to determine
    2005259   3.7000e+04   4.1293e+04   % speed direct 85
    2005259   7.4255e+04   7.4723e+04
    2005260   1.4860e+03   2.4810e+03
    2005260   2.6510e+03   3.4380e+03   % speed direct 80
    2005260   3.6040e+03   4.5000e+03   % divided into 2,
    2005260   4.6080e+03   6.6780e+03   % divided into 2,
    2005260   6.7690e+03   7.0310e+03
    2005260   7.1270e+03   7.4700e+03
    2005260   7.6740e+03   8.1000e+03
    2005260   9.4090e+03   1.0357e+04   % speed direct 55
    2005260   4.0871e+04   4.1640e+04
    2005260   5.1893e+04   5.2367e+04
    2005260   5.6300e+04   5.8200e+04   % use old times
    2005261   8.5290e+03   9.6210e+03
    2005261   9.7470e+03   1.0313e+04
    2005261   1.0683e+04   1.1169e+04
    2005261   1.1318e+04   1.2103e+04
    ];

% finally accepted ones combined by uncommented ones in ntran04redo and ntran04speed
ntran05acc = [
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

%%
% Combined accepted ones from all ETS
ntranacc = [
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

indspeed = [
    2
    27
    34
    37
    47
    ];

indrmse = [
     1
     3
     4
     5
     6
     7
     8
     9
    10
    11
    12
    13
    14
    15
    16
    17
    18
    19
    20
    21
    22
    23
    24
    25
    26
    28
    29
    30
    31
    32
    33
    35
    36
    38
    39
    40
    41
    42
    43
    44
    45
    46
    48
    50
    51
    52
    ];

indlffit = [
    49
    ];

%%
% hardcopy of angbest from rmse (col 1) and slope (col 2)
angcandihf = [
   250   260
   165    90
   220   205
   185   175
   225   275
   270   275
   235   270
   200   235
   230   285
    80    90
    65   105
   235   265
   245   285
   310   280
   250   270
   120    90
   195   180
    90    60
   140   195
   130   130
    70    70
   255   250
   250   260
   260   225
   120   105
   110   150
   155    80
   240   255
   220   255
    80    65
   115   140
   105    90
   260   270
   160    95
   260   215
   130    95
   175   230
   250   270
   195   235
   115   150
   250   255
   165   120
   245   285
   235   280
    10   350
    95    60
   155    85
   185   125
   210   280
   115    85
   250   235
   235   265
    ];

%%
angdiff = angcandihf(:,2)-angcandihf(:,1);
for i = 1: size(angdiff,1)
    if angdiff(i) > 180
        angdiff(i) = 360-angdiff(i);
    elseif angdiff(i) < -180
        angdiff(i) = angdiff(i)+360;
    end
end
angdiffs = sort(abs(angdiff));

figure
h=histogram(angdiffs,'binw',5);
h.BinEdges = h.BinEdges-2.5;
histogram(angdiffs,'binedge',h.BinEdges); hold on
plot([27.5 27.5], ylim, 'r--', 'linew',2);
perc = round(sum(angdiffs<=25)/length(angdiffs)*100);
text(0.5,0.8, sprintf('%d%% abs(angdiff)<= 25 deg', perc),...
     'unit','normalized','fontsize',12);
xlabel('abs diff. in angle (deg)');
ylabel('count'); 

%%
offset = [offset1' offset2'];
offdiff = offset(:,2)-offset(:,1);
offdiffs = sort(abs(offdiff));
figure 
histogram(offdiff,'binw',0.25); hold on
tol = 1;
plot([tol tol], ylim, 'r--', 'linew',2);
plot([-tol -tol], ylim, 'r--', 'linew',2);
perc = round(sum(abs(offdiff)<=tol)/length(offdiff)*100);
text(0.5,0.8, sprintf('%d%% abs(offdiff) <= %.1f km', perc,tol),...
     'unit','normalized','fontsize',12);
perc = round(sum(abs(offdiff)<=0.5)/length(offdiff)*100);
text(0.5,0.7, sprintf('%d%% abs(offdiff) <= %.1f km', perc,0.5),...
     'unit','normalized','fontsize',12); 
xlabel('diff. in offset (km)');
ylabel('count'); 

%%
f1.fig=figure;
f1.fig.Renderer='Painters';
widin = 4;  % maximum width allowed is 8.5 inches
htin = 4.5;   % maximum height allowed is 11 inches
set(f1.fig,'Position',[scrsz(1)+1*scrsz(3)/10 scrsz(2)+scrsz(4)/20 widin*res htin*res]);
ax = gca;
hold(ax,'on');
ax.Box = 'on';
histogram(ax,offdiffs,'binw',0.25,'facec',[.3 .3 .3]); hold on
tol = 1;
plot(ax,[tol tol], [0 50], 'b--', 'linew',2);
perc = round(sum(abs(offdiff)<=tol)/length(offdiff)*100);
text(ax,0.7,0.9, sprintf('%d%%', perc),...
     'unit','normalized','fontsize',12);
xlabel(ax,'Absolute difference in offset (km)');
ylabel(ax,'Count'); 
ylim(ax,[0 25]);
xlim(ax,[0 3.5]);
hold(ax,'off');
print(f1.fig,'-dpdf',strcat(rstpath,'/absdiffinoffset_all.pdf'));


%%
tol = 1;
indlargediff = find(abs(offdiff)>tol);

figure
ax=subplot(3,1,1);
scatter(ax, angdiff,offdiff, 10, 'b','filled'); hold on
text(0.5,0.8, 'angdiff vs. offdiff',...
     'unit','normalized','fontsize',12);
ax.Box = 'on';
grid(ax,'on');
xlabel(ax, 'diff. in angle');
ylabel(ax, 'diff. in offset');

ax=subplot(3,1,2);
scatter(ax,abs(angdiff),offdiff, 10, 'b','filled'); hold on
text(0.5,0.8, 'abs(angdiff) vs. offdiff',...
     'unit','normalized','fontsize',12);
ax.Box = 'on';
grid(ax,'on'); 
xlabel(ax, 'abs diff. in angle');
ylabel(ax, 'diff. in offset');

ax=subplot(3,1,3);
scatter(ax,abs(angdiff),abs(offdiff), 10, 'b','filled'); hold on
plot(ax,xlim, [tol tol], 'r--', 'linew',2);
text(0.5,0.8, 'abs(angdiff) vs. abs(offdiff)',...
     'unit','normalized','fontsize',12);
ax.Box = 'on';
grid(ax,'on'); 
xlabel(ax, 'abs diff. in angle');
ylabel(ax, 'abs diff. in offset');


%%
f2.fig=figure;
f2.fig.Renderer='Painters';
widin = 6;  % maximum width allowed is 8.5 inches
htin = 4;   % maximum height allowed is 11 inches
set(f2.fig,'Position',[scrsz(1)+1*scrsz(3)/10 scrsz(2)+scrsz(4)/20 widin*res htin*res]);
ax = gca;
hold(ax,'on');
ax.Box = 'on';
indne = [10,18,21,30,34,46,47];
indsw = [5,6,7,9,12,13,15,22,23,28,33,35,37,38,41,43,44,49,51,52];
indlzb = [fliplr(indne), fliplr(indsw)];
indpoorlf = [40,4,5,33,3,15,11,23,35,21,45,17,2,31,30,24,32,39,6,9,27,37,51,19,48,26,46,22,8,29,18,47];
indlargediff = [2,3,11,14,24,29,32,34,36,40,43,44,49,52];
indrmse = [1,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,28,29,30,31,32,33,35,36,38,39,40,41,42,43,44,45,46,48,50,51,52];
indspeed = [2,27,34,37,47];
indlffit = [49];
indsmalldiff = find(abs(offdiff)<=0.5);
indsmalldiff = [1,5,6,8,9,10,13,15,16,17,18,19,20,21,22,23,25,31,33,38,41,42,45,46,47,48,50,51];
indne2 = intersect(indne, indsmalldiff);
indsw2 = intersect(indsw, indsmalldiff);
for i = 1: size(trange,1)
    if ~isempty(intersect(indne,i)) && ~isempty(intersect(indrmse,i))
        s1 = scatter(ax,abs(angdiff(i)),abs(offdiff(i)), 20, 'bo','filled');
    elseif ~isempty(intersect(indne,i)) && ~isempty(intersect(indspeed,i)) 
        s2 = scatter(ax,abs(angdiff(i)),abs(offdiff(i)), 25, 'b^','filled');
    elseif ~isempty(intersect(indsw,i)) && ~isempty(intersect(indrmse,i))
        s3 = scatter(ax,abs(angdiff(i)),abs(offdiff(i)), 20, 'ro','filled');
    elseif ~isempty(intersect(indsw,i)) && ~isempty(intersect(indspeed,i))
        s4 = scatter(ax,abs(angdiff(i)),abs(offdiff(i)), 25, 'r^','filled');
    elseif ~isempty(intersect(indsw,i)) && ~isempty(intersect(indlffit,i))
        s5 = scatter(ax,abs(angdiff(i)),abs(offdiff(i)), 25, 'rs','filled');    
    end
    
end

% for i = 1: size(trange,1)
%     if ~isempty(intersect(indlzb,i)) && isempty(intersect(indlargediff,i))
%         scatter(ax,abs(angdiff(i)),abs(offdiff(i)), 20, 'ko','filled');
%     elseif ~isempty(intersect(indlzb,i)) && ~isempty(intersect(indlargediff,i)) && ...
%            isempty(intersect(indpoorlf,i)) 
%         scatter(ax,abs(angdiff(i)),abs(offdiff(i)), 25, 'k^','filled');
%     elseif isempty(intersect(indlzb,i))
%         scatter(ax,abs(angdiff(i)),abs(offdiff(i)), 20, 'ko','linew',1);
%     end
%     
% end

% scatter(ax,abs(angdiff),abs(offdiff), 10, 'b','filled');
tol = 1;
plot(ax,[0 180], [tol tol], 'k--', 'linew',1.5);
plot(ax,[0 180], [0.5 0.5], 'k-.', 'linew',1);
perc = round(sum(abs(offdiff)<=tol)/length(offdiff)*100);
% text(ax,0.7,0.9, sprintf('%d%%', perc),...
%      'unit','normalized','fontsize',12);
lgd = legend(ax,[s1,s2,s3,s4,s5],{'ENE propagation, min-RMSE','ENE propagation, max-speed',...
                'WSW propagation, min-RMSE','WSW propagation, max-speed',...
                'WSW propagation, min-RMSE(LF)'},'FontSize',7.5,...
       'Position',[0.22 0.75 0.2 0.15]); %,'Box','off'
set(lgd.BoxFace, 'ColorType','truecoloralpha', 'ColorData',uint8(255*[1;1;1;.8]));  

xlabel(ax,'Absolute difference in direction (deg)');
ylabel(ax,'Absolute difference in offset (km)');
xlim(ax,[0 80]);
ylim(ax,[0 2]);
yticks(ax,0:.5:2);
hold(ax,'off');
print(f2.fig,'-dpdf',strcat(rstpath,'/absdiffinoffsetvsangle_used.pdf'));


%%
f2.fig=figure;
f2.fig.Renderer='Painters';
widin = 6;  % maximum width allowed is 8.5 inches
htin = 6;   % maximum height allowed is 11 inches
set(f2.fig,'Position',[scrsz(1)+1*scrsz(3)/10 scrsz(2)+scrsz(4)/20 widin*res htin*res]);
ax = gca;
hold(ax,'on');
ax.Box = 'on';
indne = [10,18,21,30,34,46,47];
indsw = [5,6,7,9,12,13,15,22,23,28,33,35,37,38,41,43,44,49,51,52];
indlzb = [fliplr(indne), fliplr(indsw)];
indpoorlf = [40,4,5,33,3,15,11,23,35,21,45,17,2,31,30,24,32,39,6,9,27,37,51,19,48,26,46,22,8,29,18,47];
indlargediff = [2,3,11,14,24,29,32,34,36,40,43,44,49,52];
indrmse = [1,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,28,29,30,31,32,33,35,36,38,39,40,41,42,43,44,45,46,48,50,51,52];
indspeed = [2,27,34,37,47];
indlffit = [49];
indsmalldiff = find(abs(offdiff)<=0.5);
indsmalldiff = [1,5,6,8,9,10,13,15,16,17,18,19,20,21,22,23,25,31,33,38,41,42,45,46,47,48,50,51];
indne2 = intersect(indne, indsmalldiff);
indsw2 = intersect(indsw, indsmalldiff);
for i = 1: size(trange,1)
    if ~isempty(intersect(indne,i)) && ~isempty(intersect(indrmse,i))
        s1 = scatter(ax,abs(angdiff(i)),abs(offdiff(i)), 20, 'bo','filled');
    elseif ~isempty(intersect(indne,i)) && ~isempty(intersect(indspeed,i)) 
        s2 = scatter(ax,abs(angdiff(i)),abs(offdiff(i)), 25, 'b^','filled');
    elseif ~isempty(intersect(indsw,i)) && ~isempty(intersect(indrmse,i))
        s3 = scatter(ax,abs(angdiff(i)),abs(offdiff(i)), 20, 'ro','filled');
    elseif ~isempty(intersect(indsw,i)) && ~isempty(intersect(indspeed,i))
        s4 = scatter(ax,abs(angdiff(i)),abs(offdiff(i)), 25, 'r^','filled');
    elseif ~isempty(intersect(indsw,i)) && ~isempty(intersect(indlffit,i))
        s5 = scatter(ax,abs(angdiff(i)),abs(offdiff(i)), 25, 'rs','filled');    
    elseif isempty(intersect(indlzb,i)) && ~isempty(intersect(indrmse,i))
        s6 = scatter(ax,abs(angdiff(i)),abs(offdiff(i)), 20, 'ko','linew',1);
    elseif isempty(intersect(indlzb,i)) && ~isempty(intersect(indspeed,i))
        s7 = scatter(ax,abs(angdiff(i)),abs(offdiff(i)), 25, 'k^','linew',1);
    end
    
end
tol = 1;
plot(ax,[0 180], [tol tol], 'k--', 'linew',1.5);
perc = round(sum(abs(offdiff)<=tol)/length(offdiff)*100);
% text(ax,0.7,0.9, sprintf('%d%%', perc),...
%      'unit','normalized','fontsize',12);
lgd = legend(ax,[s1,s2,s3,s4,s5,s6,s7],{'ENE propagation, min-RMSE','ENE propagation, max-speed',...
                'WSW propagation, min-RMSE','WSW propagation, max-speed',...
                'WSW propagation, min-RMSE(LF)', 'Other RTMs, min-RMSE',...
                'Other RTMs, max-speed'},'FontSize',7.5,...
       'Position',[0.22 0.59 0.2 0.2]); %,'Box','off'
set(lgd.BoxFace, 'ColorType','truecoloralpha', 'ColorData',uint8(255*[1;1;1;.8]));  

xlabel(ax,'Absolute difference in direction (deg)');
ylabel(ax,'Absolute difference in offset (km)');
xlim(ax,[0 80]);
% ylim(ax,[0 2]);
% yticks(ax,0:.5:2);
hold(ax,'off');
print(f2.fig,'-dpdf',strcat(rstpath,'/absdiffinoffsetvsangle_all.pdf'));



%%
largediff = [    
    2004197   4.0974e+04   4.2205e+04   % speed direct 90, y, larger pear
    2004197   5.4453e+04   5.5867e+04   % speed direct 205, y, but use rmse 220, poor LF constraint?
    2004198   6.2966e+04   6.3885e+04   % acceptted, y, though it is not perfect in both passbands
    2004199   4.1670e+03   6.2660e+03   % speed direct 280, y, larger pear
    2004201   4.3700e+02   1.6290e+03   % speed direct 225, y, but use rmse 260, poor LF constraint?
    2005256   8.0330e+03   9.8280e+03   % change end time, y, poor LF constraint?
    2005256   2.8080e+04   2.9172e+04   % change start time, y, poor LF constraint?
    2005256   3.6900e+04   3.8520e+04   % change times, y, mainly due to HF?
    2005256   6.2640e+04   6.4865e+04   % change start time, y, check speed direc, y, but use rmse 130, mainly due to HF?
    2005257   1.6920e+04   1.8655e+04   % change start time, y, poor LF constraint?
    2005257   6.1449e+04   6.4012e+04   % accepted, y, mainly due to HF?
    2005257   7.3440e+04   7.8840e+04   % change end time, y, but check speed direc, y, but use rmse 235, mainly due to HF?
    2005259   5.2320e+03   7.7140e+03   % check OLD param. to determine, y, but check speed direc, larger pear, use old direc 255, mainly due to HF?
    2005260   5.6300e+04   5.8200e+04   % use old times, y, okay one, mainly due to HF?
    ];

indlargediff = [
     2
     3
    11
    14
    24
    29
    32
    34
    36
    40
    43
    44
    49
    52
    ];

% if exlucding detections at 002 region
indlargediffnew = [2,3,7,11,12,14,29,32,34,36,40,43,49];


%%
% ones within angle range and excluding those at distant fams, eg. 002, 010
indne = setdiff(find(angbest>=40 & angbest<=100),[27]); % [2,10,11,18,21,30,34,46,47]
% ones have poor lf constraints
indpoorlf = sortselfbest(sortselfbest(:,2)>=3, 1); %[40,4,5,33,3,15,11,23,35,21,45,17,2,31,30,24,32,39,6,9,27,37,51,19,48,26,46,22,8,29,18,47]
% ones have large difference in offsets due to a different selection of propagation direction
tol = 1;
indlargediff = find(abs(offdiff)>tol);  %[2,3,11,14,24,29,32,34,36,40,43,44,49,52]
% so you can have several options:
% 1. all within angles
indplt = indne;  % [2,10,11,18,21,30,34,46,47]
% 2. all excluding EITHER poor LF constraints OR there will be very different offsets
tmp = union(indpoorlf, indlargediff);
indplt = setdiff(indne, tmp);   % [10]
% 3. all excluding ones have large difference in offsets
indplt = setdiff(indne, indlargediff);  % [10,18,21,30,46,47]
% 4. all excluding ones have large difference in offsets due to LF constraint
tmp = intersect(indlargediff, indpoorlf);
indplt = setdiff(indne, tmp);   % [10,18,21,30,34,46,47]


% ones within angle range and excluding those at distant fams, eg. 002, 010
indsw = setdiff(find(angbest>=220 & angbest<=280),[1]); % [3,5,6,7,9,12,13,15,22,23,24,28,29,33,35,37,38,41,43,44,49,51,52]
% ones have poor lf constraints
indpoorlf = sortself(sortself(:,2)>=3, 1); %[40,4,5,33,3,15,11,23,35,21,45,17,2,31,30,24,32,39,6,9,27,37,51,19,48,26,46,22,8,29,18,47]

% ones have large difference in offsets due to a different selection of propagation direction
tol = 1;
indlargediff = find(abs(offdiff)>tol);  %[2,3,11,14,24,29,32,34,36,40,43,44,49,52]
% so you can have several options:
% 1. all within angles
indplt = indsw;  % [3,5,6,7,9,12,13,15,22,23,24,28,29,33,35,37,38,41,43,44,49,51,52]
% 2. all excluding EITHER poor LF constraints OR there will be very different offsets
tmp = union(indpoorlf, indlargediff);
indplt = setdiff(indsw, tmp); % [7,12,13,28,38,41]
% 3. all excluding ones have large difference in offsets
indplt = setdiff(indsw, indlargediff);  % [5,6,7,9,12,13,15,22,23,28,33,35,37,38,41,51]
% 4. all excluding ones have large difference in offsets due to LF constraint
tmp = intersect(indlargediff, indpoorlf);
indplt = setdiff(indsw, tmp);   % [5,6,7,9,12,13,15,22,23,28,33,35,37,38,41,43,44,49,51,52]





