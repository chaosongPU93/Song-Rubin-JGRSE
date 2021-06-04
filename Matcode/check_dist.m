%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% THIS is the script to check the distance between lfe and each station in
% the trio to see if it is consistent with the 4th col of rotation parameters
% from calc_rots_allfam. Use the function 'distaz' written by Chao.
%
%   
%
% Modified by Chao Song, chaosong@princeton.edu
% First created date:   2019/11/11
% Last modified date:   2019/11/11
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

format short
% clear
% close all
set(0,'DefaultFigureVisible','on');
% set(0,'DefaultFigureVisible','off');   % switch to show the plots or not
scrsz=get(0,'ScreenSize');   % for example, scrsz=[1,1,1680,1050]

% set path
workpath = getenv('ALLAN');
temppath = strcat(workpath, '/templates/');
if ~exist(temppath, 'dir')
    mkdir(temppath);
end
% datapath = workpath;
datapath = strcat(workpath,'/data-no-resp');

% station loc file
dtdir = ('/home/data2/chaosong/matlab/allan/2004/JUL/');
flist= ('datalist');
outf = ('staloc.txt');
stainfo = getstainfo(dtdir, flist, outf);
ind=[5,8,6,3,11,4,2];   % PGC,SSIB,SILB,LZB,TWKB,MGCB,KLNB
stainfo = stainfo(ind,:);

% lfe fam loc file
bosdir = ('/home/data2/chaosong/matlab/allan/BOSTOCK');
lfefnm = ('newlfeloc');
lfeloc = load(fullfile(bosdir, lfefnm));    
%%% lfeloc is loctions of all fams, but not necessarily each fam has detections

% lfe fam pool
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
nfam = length(nfampool);

[~,~,ind] = intersect(str2num(nfampool(:,:)), lfeloc(:, 1),'stable');
lfelocuse = lfeloc(ind,:);

% rotation parameter file
famp1 = ['002';'043';'141';'047';'010';'099';'017'];
famp2 = ['144';'068';'125';'147';];

ifam = 11;
fam = nfampool(ifam,:);

lo = 0.5;
hi = 6.5;
sps = 40;
bef = 25;
aft = 35;
nsta = 7;
if ismember(fam,famp1,'rows')
    PREFIX = strcat(fam,'rot_para',num2str(lo),'-',num2str(hi),'_',num2str(bef),'_',...
                    num2str(aft),'LZB');
elseif ismember(fam,famp2,'rows')
    PREFIX = strcat(fam,'rot_para',num2str(lo),'-',num2str(hi),'_',num2str(bef),'_',...
                    num2str(aft),'TWKB');
end
rots = load(strcat(datapath, '/split_chao/',PREFIX,'_',num2str(sps),'sps',num2str(nsta),'stas'));

for i = 1: size(stainfo,1)
    stla = str2double(stainfo(i,2));
    stlo = str2double(stainfo(i,3));
    eqla = lfelocuse(ifam,2);
    eqlo = lfelocuse(ifam,3);
    [distk(i),distd(i),~,~] = distaz(stla,stlo,eqla,eqlo);
    [distd2(i),~] = distance(stla,stlo,eqla,eqlo);
end

stadist = (distk-distk(4))'
staoff = rots(:,4)-rots(4,4)

%%
figure
scatter(lfelocuse(:,3), lfelocuse(:,2), 30, 'r', 'filled', 'MarkerEdgeColor', 'k'); hold on
text(lfelocuse(:,3)-0.02, lfelocuse(:,2)+0.02, nfampool(:,:),'fontsize',10);
for i = 1: size(stainfo,1)
    plot(str2num(stainfo(i,3)),str2num(stainfo(i,2)),'^','MarkerFaceColor','k','MarkerSize',8,'MarkerEdgeColor', 'k');
    text(str2num(stainfo(i,3))+0.01,str2num(stainfo(i,2))+0.01,stainfo(i,1));
end
xlabel('Longitude');
ylabel('Latitude');
box on
grid on










    
    
    
    
    
    