function stainfo=getstainfo(dtdir, flist, outf)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% This is the function to get the location information of the stations in
% your data, and output the info file
% 
% Features:
%   
%
%
% By Chao Song, chaosong@princeton.edu
% First created date:   2019/04/24
% Last modified date:   2019/04/24
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% read data from a random directory.

defval('dtdir', '/home/data2/chaosong/matlab/allan/2004/JUL/')
defval('flist','datalist')
defval('outf','staloc.txt')

fid = fopen(fullfile(dtdir, flist), 'r');
datalist = textscan(fid, '%s \n');
nsta = size(datalist{1},1);
stla = zeros(nsta, 1);
stlo = zeros(nsta, 1);
stnm = strings(nsta, 1);
for i = 1: nsta
%     i=1;
    stafnm = fullfile(dtdir,datalist{1}{i});

    [~, sachdr, ~, ~, ~] = readsac(stafnm, 0, 'l');   % 'l' means linux
    stla(i) = sachdr.STLA;
    stlo(i) = sachdr.STLO;
    stnm(i) = sachdr.KSTNM;

end

stainfo = [stnm stla stlo];

fid = fopen(fullfile(dtdir, outf), 'w+');
fprintf(fid,'%s %s %s \n',stainfo');
fclose(fid);

