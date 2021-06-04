% This script stores the unused code segments from 'mig_linear_fit_LZB_addfam.m'
%
% By Chao Song, chaosong@princeton.edu
% First created date:   2020/10/19
% Last modified date:   2020/10/19




%% shift the histogram of opposite migration group to their average of median
%%% account for the possible incorrect filtering correction which results in the discrepancy of the
%%% medians of the distance distribtuion of migration group that has opposite directions
%%% NOTE: using methods to obtain median of data with unequal weights 

%%% 1. for SW and NE group at main LZB region, i.e. the NW part of the study region
medsw = wt_median(vdistsw, wtsw);
medne = wt_median(vdistne, wtne);
medave = (medsw+medne)/2;
% meaning shifting of the raw data, if you use normfit, meaning the shift of mean, but preserve
% weights 
vdistsftsw = vdistsw - (medsw-medave);
vdistsftne = vdistne - (medne-medave);
% also meaning shifting of the xloc of the pdf, but perserve the pdf
pdfxlocsftsw = pdfxlocsw - (medsw-medave);
pdfxlocsftne = pdfxlocne - (medne-medave);
% combine them together
vdistsftswne = [vdistsftsw; vdistsftne];
wtsftswne = [wtsw; wtne];
% but there is a slight difference between medsw and median_of_wt_data(vdistsftswne, wtsftswne)

% get the pdf of the combined data
[pdfxlocsftswne, ~, ~, pdfvalsftswne] = weighted_bar_pdf(vdistsftswne, wtsftswne, [], 'int');
% shift the combined data to make the median to be zero, but preserve weights, test for Null
% hypothesis 
vdisttmpswne = vdistsftswne - (medave -0);
wttmpswne = wtsftswne;
% also meaning shifting of the xloc of the pdf, but perserve the pdf
pdfxloctmpswne = pdfxlocsftswne - (medave -0);
pdfvaltmpswne = pdfvalsftswne;
% we could do a test here, since the median of the combined data should be 0 now, but in fact it is
% not due to the reason above, but the difference should be very small
disp(wt_median(vdisttmpswne, wttmpswne));

%%% 2. for SW and NE group at fam 002 region
med002sw = wt_median(vdist002sw, wt002sw);
med002ne = wt_median(vdist002ne, wt002ne);
med002ave = (med002sw+med002ne)/2;
% meaning shifting of the raw data, if you use normfit, meaning the shift of mean, but preserve
% weights 
vdist002sftsw = vdist002sw - (med002sw-med002ave);
vdist002sftne = vdist002ne - (med002ne-med002ave);
% also meaning shifting of the xloc of the pdf, but perserve the pdf
% pdfxloc002sftsw = pdfxloc002sw - (med002sw-med002ave);
% pdfxloc002sftne = pdfxloc002ne - (med002ne-med002ave);
% combine them together
vdist002sftswne = [vdist002sftsw; vdist002sftne];
wt002sftswne = [wt002sw; wt002ne];
% but there is a slight difference between medsw and median_of_wt_data(vdistsftswne, wtsftswne)

% get the pdf of the combined data
[pdfxloc002sftswne, ~, ~, pdfval002sftswne] = weighted_bar_pdf(vdist002sftswne, wt002sftswne, [],...
                                                               'int');
% shift the combined data to make the median to be zero, but preserve weights, test for Null
% hypothesis 
vdist002tmpswne = vdist002sftswne - (med002ave -0);
wt002tmpswne = wt002sftswne;
% also meaning shifting of the xloc of the pdf, but perserve the pdf
pdfxloc002tmpswne = pdfxloc002sftswne - (med002ave -0);
pdfval002tmpswne = pdfval002sftswne;
% we could do a test here, since the median of the combined data should be 0 now, but in fact it is
% not due to the reason above, but the difference should be very small
disp(wt_median(vdist002tmpswne, wt002tmpswne));


%% plot the combined migrations
%%% 1. for SW and NE group at main LZB region, i.e. the NW part of the study region
f.fig=figure;
widin = 8;  % maximum width allowed is 8.5 inches
htin = 4.5;   % maximum height allowed is 11 inches
set(f.fig,'Position',[scrsz(1)+1*scrsz(3)/10 scrsz(2)+scrsz(4)/20 widin*res htin*res]);
nrow = 1;
ncol = 2;
for isub = 1:nrow*ncol
    f.ax(isub) = subplot(nrow,ncol,isub);
end

%%% reposition
set(f.ax(1), 'position', [ 0.1, 0.1, 0.35, 0.75]);

% subplot 1
hold(f.ax(1),'on');
xran = [-15 20];
yran = [-15 15];
indplt = union(indsw,indne);    % index to plot
[f.ax(1),c] = plt_mig_prop_onmap(f.ax(1), indplt, ranvechf98, medallhf, angbest, xran,yran);
hold(f.ax(1),'off');

% subplot 2 
ywid = f.ax(1).Position(4)+f.ax(1).Position(2)-(c.Position(2));
set(f.ax(2), 'position', [ 0.55, c.Position(2), 0.35, ywid]);
hold(f.ax(2),'on');
[f.ax(2), ~, ~, ~, ~, ~] = plt_weighted_dist(f.ax(2), vdistsftswne,wtsftswne, [],'dec');
text(f.ax(2),0.5,0.9,strcat({'med = '},sprintf('%.2f',median_of_wt_data(vdistsftswne,wtsftswne)),...
     {' km'}),'fontsize',12,'unit','normalized');
hold(f.ax(2),'off');

%%% save figure
print(f.fig,'-dpdf',strcat(rstpath,'/vertical_dist_shifted_mainSWNE_migs',num2str(nfam),'.pdf'));


%%% 2. for SW and NE group at fam 002 region
f.fig=figure;
widin = 8;  % maximum width allowed is 8.5 inches
htin = 4.5;   % maximum height allowed is 11 inches
set(f.fig,'Position',[scrsz(1)+1*scrsz(3)/10 scrsz(2)+scrsz(4)/20 widin*res htin*res]);
nrow = 1;
ncol = 2;
for isub = 1:nrow*ncol
    f.ax(isub) = subplot(nrow,ncol,isub);
end

%%% reposition
set(f.ax(1), 'position', [ 0.1, 0.1, 0.35, 0.75]);

% subplot 1
hold(f.ax(1),'on');
xran = [-5 20];
yran = [-10 10];
indplt = union(ind002sw,ind002ne);    % index to plot
[f.ax(1),c] = plt_mig_prop_onmap(f.ax(1), indplt, ranvechf98, medallhf, angbest, xran,yran);
hold(f.ax(1),'off');

% subplot 2 
ywid = f.ax(1).Position(4)+f.ax(1).Position(2)-(c.Position(2)+0.014);
set(f.ax(2), 'position', [ 0.55, c.Position(2)-0.00, 0.35, ywid]);
hold(f.ax(2),'on');
[f.ax(2), ~, ~, ~, ~, ~] = plt_weighted_dist(f.ax(2), vdist002sftswne, wt002sftswne, [],'dec');
text(f.ax(2),0.5,0.9,strcat({'med = '},sprintf('%.2f',median_of_wt_data(vdist002sftswne,...
     wt002sftswne)),{' km'}),'fontsize',12,'unit','normalized');
hold(f.ax(2),'off');

%%% save figure
print(f.fig,'-dpdf',strcat(rstpath,'/vertical_dist_shifted_SWNE002_migs',num2str(nfam),'.pdf'));


%%
%%% plot 1, SW migs at NW regions of the above plot
f12.fig=figure;
widin = 8;  % maximum width allowed is 8.5 inches
htin = 4.5;   % maximum height allowed is 11 inches
set(f12.fig,'Position',[scrsz(1)+1*scrsz(3)/10 scrsz(2)+scrsz(4)/20 widin*res htin*res]);
nrow = 1;
ncol = 2;
for isub = 1:nrow*ncol
    f12.ax(isub) = subplot(nrow,ncol,isub);
end

%%% reposition
set(f12.ax(1), 'position', [ 0.1, 0.1, 0.35, 0.75]);
% set(f.ax(2), 'position', [ 0.55, 0.1, 0.35, 0.75]);

% subplot 1
hold(f12.ax(1),'on');
xran = [-15 10];
yran = [-10 10];
% indplt = indsw(11);    % index to plot
% indplt = 42;
indplt = indsw;    % index to plot
[f12.ax(1),c] = plt_mig_prop_onmap(f12.ax(1), indplt, ranvechf98, medallhf, angbest, xran,yran);
hold(f12.ax(1),'off');

% subplot 2 
ywid = f12.ax(1).Position(4)+f12.ax(1).Position(2)-(c.Position(2)+0.014);
set(f12.ax(2), 'position', [ 0.55, c.Position(2)-0.00, 0.35, ywid]);
hold(f12.ax(2),'on');
[f12.ax(2), barsw, pdfxlocsw, ~, ~, pdfvalsw] = plt_weighted_dist(f12.ax(2), vdistsw, wtsw, [],'dec');
hold(f12.ax(2),'off');

%%% save figure
print(f12.fig,'-dpdf',strcat(rstpath,'/vertical_dist_SW_migs',num2str(nfam),'.pdf'));


%%
%%% plot 2, NE migs at NW regions of the above plot
f.fig=figure;
widin = 8;  % maximum width allowed is 8.5 inches
htin = 4.5;   % maximum height allowed is 11 inches
set(f.fig,'Position',[scrsz(1)+1*scrsz(3)/10 scrsz(2)+scrsz(4)/20 widin*res htin*res]);
nrow = 1;
ncol = 2;
for isub = 1:nrow*ncol
    f.ax(isub) = subplot(nrow,ncol,isub);
end

%%% reposition
set(f.ax(1), 'position', [ 0.1, 0.1, 0.35, 0.75]);

% subplot 1
hold(f.ax(1),'on');
xran = [-15 10];
yran = [-10 10];
indplt = indne;    % index to plot
[f.ax(1),c] = plt_mig_prop_onmap(f.ax(1), indplt, ranvechf98, medallhf, angbest, xran,yran);
hold(f.ax(1),'off');

% subplot 2 
ywid = f.ax(1).Position(4)+f.ax(1).Position(2)-(c.Position(2)+0.014);
set(f.ax(2), 'position', [ 0.55, c.Position(2)-0.00, 0.35, ywid]);
hold(f.ax(2),'on');
[f.ax(2), barne, pdfxlocne, ~, ~, pdfvalne] = plt_weighted_dist(f.ax(2), vdistne, wtne, [],'dec');
hold(f.ax(2),'off');

%%% save figure
print(f.fig,'-dpdf',strcat(rstpath,'/vertical_dist_NE_migs',num2str(nfam),'.pdf'));


%%
%%% plot 3, NW + SE at region of 002, in single plot
f.fig=figure;
widin = 8;  % maximum width allowed is 8.5 inches
htin = 4.5;   % maximum height allowed is 11 inches
set(f.fig,'Position',[scrsz(1)+1*scrsz(3)/10 scrsz(2)+scrsz(4)/20 widin*res htin*res]);
nrow = 1;
ncol = 2;
for isub = 1:nrow*ncol
    f.ax(isub) = subplot(nrow,ncol,isub);
end

%%% reposition
set(f.ax(1), 'position', [ 0.1, 0.1, 0.35, 0.75]);

% subplot 1
hold(f.ax(1),'on');
xran = [-5 20];
yran = [-10 10];
indplt = ind002;    % index to plot
[f.ax(1),c] = plt_mig_prop_onmap(f.ax(1), indplt, ranvechf98, medallhf, angbest, xran,yran);
hold(f.ax(1),'off');

% subplot 2 
ywid = f.ax(1).Position(4)+f.ax(1).Position(2)-(c.Position(2)+0.014);
set(f.ax(2), 'position', [ 0.55, c.Position(2)-0.00, 0.35, ywid]);
hold(f.ax(2),'on');
[f.ax(2), bar002, pdfxloc002, ~, ~, pdfval002] = plt_weighted_dist(f.ax(2), vdist002, wt002, ...
                                                                   [],'dec');
hold(f.ax(2),'off');

%%% save figure
print(f.fig,'-dpdf',strcat(rstpath,'/vertical_dist_SWNE002_migs',num2str(nfam),'.pdf'));


%%
%%% plot 3.5, NW regions + region of 002, in single plot
f.fig=figure;
widin = 8;  % maximum width allowed is 8.5 inches
htin = 4.5;   % maximum height allowed is 11 inches
set(f.fig,'Position',[scrsz(1)+1*scrsz(3)/10 scrsz(2)+scrsz(4)/20 widin*res htin*res]);
nrow = 1;
ncol = 2;
for isub = 1:nrow*ncol
    f.ax(isub) = subplot(nrow,ncol,isub);
end

%%% reposition
set(f.ax(1), 'position', [ 0.1, 0.1, 0.35, 0.75]);

% subplot 1
hold(f.ax(1),'on');
xran = [-15 20];
yran = [-15 15];
indplt = indswne;    % index to plot
[f.ax(1),c] = plt_mig_prop_onmap(f.ax(1), indplt, ranvechf98, medallhf, angbest, xran,yran);
hold(f.ax(1),'off');

% subplot 2 
ywid = f.ax(1).Position(4)+f.ax(1).Position(2)-(c.Position(2));
set(f.ax(2), 'position', [ 0.55, c.Position(2), 0.35, ywid]);
hold(f.ax(2),'on');
[f.ax(2), barswne, pdfxlocswne, ~, ~, pdfvalswne] = plt_weighted_dist(f.ax(2), vdistswne, wtswne,...
                                                                      [],'dec');
hold(f.ax(2),'off');

%%% save figure
print(f.fig,'-dpdf',strcat(rstpath,'/vertical_dist_SWNE_migs',num2str(nfam),'.pdf'));


%%
%%% plot 4, all others, in single plot
f.fig=figure;
widin = 8;  % maximum width allowed is 8.5 inches
htin = 4.5;   % maximum height allowed is 11 inches
set(f.fig,'Position',[scrsz(1)+1*scrsz(3)/10 scrsz(2)+scrsz(4)/20 widin*res htin*res]);
nrow = 1;
ncol = 2;
for isub = 1:nrow*ncol
    f.ax(isub) = subplot(nrow,ncol,isub);
end

%%% reposition
set(f.ax(1), 'position', [ 0.1, 0.1, 0.35, 0.75]);

% subplot 1
hold(f.ax(1),'on');
xran = [-15 20];
yran = [-15 15];
indplt = indelse;    % index to plot
[f.ax(1),c] = plt_mig_prop_onmap(f.ax(1), indplt, ranvechf98, medallhf, angbest, xran,yran);
hold(f.ax(1),'off');

% subplot 2 
ywid = f.ax(1).Position(4)+f.ax(1).Position(2)-(c.Position(2));
set(f.ax(2), 'position', [ 0.55, c.Position(2), 0.35, ywid]);
hold(f.ax(2),'on');
[f.ax(2), barelse, pdfxlocelse, ~, ~, pdfvalelse] = plt_weighted_dist(f.ax(2), vdistelse, wtelse,...
                                                                      [],'dec');
hold(f.ax(2),'off');

%%% save figure
print(f.fig,'-dpdf',strcat(rstpath,'/vertical_dist_other_migs',num2str(nfam),'.pdf'));


%% set and save random seed state for reproduction 
N = 1000;   % sampling times
seed = [];
for i = 1: N
%     seed{i} = rng(i,'twister');    % fix and perserve the random seed state fro reproductivity
    seed{i} = rng('shuffle','twister');
end

% %% Theoretical test for random number generation from the custom PDF
% %%% CHECK the notes 'Meetingnotes/2020-09-25' and 'Research/latest'for detailed explanation
% %%%%%%%%%%%%%%%%%%%%%%%% THIS portion is important for testing the theory %%%%%%%%%%%%%%%%%%%%%%%%%%
% i = 1;
% dimsize = [ceil(sum(wtsw)),1];  % size of desired randum numbers per sampling
% rdnum = [];
% %%% 1. assume pdf is normal-like, and try to fit it with a gaussian
% rng(seed{i});
% [randnum,muHat,sigmaHat,~,~] = randgen_normal(vdistsw, wtsw, dimsize);
% rdnum = [rdnum randnum];
% % some sense of the misfit
% [pdfmisfit1, pdfmisfit2, fighdl] = compare_distribution2(randnum,muHat,sigmaHat,vdistsw, wtsw);
% 
% 
% %%% 3. take pdf as customized, use function 'normpdf' NOT written by me
% rng(seed{i});
% randnum = randpdf(pdfval, pdfxloc, dimsize);
% rdnum = [rdnum randnum];
% 
% %%% 4. take pdf as customized, use built-in function 'ksdensity'
% rng(seed{i});
% randnum = ksdensity(vdistsw, rand(dimsize), 'function', 'icdf', 'weight', wtsw);
% rdnum = [rdnum randnum];
% 
% %%% check the distribution of generated random numbers to see if it accords with the pdf of data
% pdfmisfit = zeros(size(rdnum,2), 1);
% titlestr = {'rand gen from normal fitting';
%             'rand gen from custom pdf';
%             'rand gen from ksdensity';
%            };
% for ii = 1: size(rdnum,2)
%     binw = pdfxloc(2)-pdfxloc(1);
%     [pdfmisfit(ii), f55] = compare_distribution1(rdnum(:,ii), pdfxloc, pdfval, binw);
%     supertit(f55.ax, string(titlestr{ii}),12);
%     
% end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Application of random number generation for real cases, turn testing into functions
%%% 1. for SW and NE group at main LZB region, i.e. the NW part of the study region
% compute the median of the randoms generated and the misfits of PDFs using the above 3 ways
% in order: normpdf, ksdensity, normfit
[medrandswne,pdfmisfswne] = randgen_custom(seed,N,vdisttmpswne,wttmpswne,pdfxloctmpswne,pdfvaltmpswne);

% plot the distribution of generated random numbers
[probswne] = plt_random_median_simple(medrandswne,median_of_wt_data(vdistsftswne,wtsftswne));


%%% 2. for SW and NE group at fam 002 region
% compute the median of the randoms generated and the misfits of PDFs
[medrand002swne,pdfmisf002swne] = randgen_custom(seed,N,vdist002tmpswne,wt002tmpswne,...
                                                 pdfxloc002tmpswne,pdfval002tmpswne);
 
% plot the distribution of generated random numbers
[prob002swne,f] = plt_random_median_simple(medrand002swne,median_of_wt_data(vdist002sftswne,...
                                           wt002sftswne));



%% plot the combined migrations, unshifted all
%%% 1. for SW and NE group at main LZB region, i.e. the NW part of the study region
f.fig=figure;
widin = 8;  % maximum width allowed is 8.5 inches
htin = 4.5;   % maximum height allowed is 11 inches
set(f.fig,'Position',[scrsz(1)+1*scrsz(3)/10 scrsz(2)+scrsz(4)/20 widin*res htin*res]);
nrow = 1;
ncol = 2;
for isub = 1:nrow*ncol
    f.ax(isub) = subplot(nrow,ncol,isub);
end

%%% reposition
set(f.ax(1), 'position', [ 0.1, 0.1, 0.35, 0.75]);

% subplot 1
hold(f.ax(1),'on');
xran = [-15 20];
yran = [-15 15];
indplt = union(indsw,indne);    % index to plot
[f.ax(1),c] = plt_mig_prop_onmap(f.ax(1), indplt, ranvechf98, medallhf, angbest, xran,yran);
hold(f.ax(1),'off');

% subplot 2 
binw = 1;

ywid = f.ax(1).Position(4)+f.ax(1).Position(2)-(c.Position(2));
set(f.ax(2), 'position', [ 0.55, c.Position(2), 0.35, ywid]);

meansw = wt_mean(vdistsw,wtsw);
meanne = wt_mean(vdistne,wtne);
sigmasw = sqrt(wt_var(vdistsw,wtsw,2));
neff = sum(wtsw);
conf = 99;
CIsw = confidence_interval(meansw,sigmasw,neff,conf);
sigmane = sqrt(wt_var(vdistne,wtne,2));
neff = sum(wtne);
% conf = 99;
CIne = confidence_interval(meanne,sigmane,neff,conf);

hold(f.ax(2),'on');
[f.ax(2), barsw, pdfxlocsw, ~, ~, pdfvalsw] = plt_weighted_dist(f.ax(2), vdistsw, wtsw, binw,'dec');
barsw(1).FaceAlpha = 1;
[f.ax(2), barne, pdfxlocne, ~, ~, pdfvalne] = plt_weighted_dist(f.ax(2), vdistne, wtne, binw,'dec');
barne(1).FaceColor = [0 1 1];
plot(f.ax(2),[meansw meansw],[-100 100],'r--','linew',1.5);    % median of wt dist
plot(f.ax(2),[meanne meanne],[-100 100],'b--','linew',1.5);    % median of wt dist
errorbar(f.ax(2),meansw,0.25,CIsw(1)-meansw,CIsw(2)-meansw,'horizontal','o','markersize',3,'color',...
         'r','linew',1,'MarkerEdgeColor','r','MarkerFaceColor','r');
errorbar(f.ax(2),meanne,0.25,CIne(1)-meanne,CIne(2)-meanne,'horizontal','o','markersize',3,'color',...
         'b','linew',1,'MarkerEdgeColor','b','MarkerFaceColor','b');
xvect = [-7 -12];
yvect = [0.015 0.015];
drawArrow(f.ax(2),xvect,yvect,f.ax(2).XLim,f.ax(2).YLim,'linewidth',1.5,'linestyle','-',...
          'color',[0.6 0.6 0.6]);
xvect = [7 12];
yvect = [0.015 0.015];
drawArrow(f.ax(2),xvect,yvect,f.ax(2).XLim,f.ax(2).YLim,'linewidth',1.5,'linestyle','-',...
          'color',[0.6 0.6 0.6]);
text(f.ax(2),0.05,0.1,strcat({'WSW'}),'fontsize',11,'unit','normalized');
text(f.ax(2),0.85,0.1,strcat({'ENE'}),'fontsize',11,'unit','normalized');
legend(f.ax(2),[barsw,barne],{'WSW propagating','ENE propagating'},'fontsize',7);
hold(f.ax(2),'off');

%%% save figure
print(f.fig,'-dpdf',strcat(rstpath,'/vertical_dist_unshifted_mainSW+NE_migs',num2str(nfam),'.pdf'));

ff.fig=figure;
widin = 3.5;  % maximum width allowed is 8.5 inches
htin = 3.8;   % maximum height allowed is 11 inches
set(ff.fig,'Position',[scrsz(1)+1*scrsz(3)/10 scrsz(2)+scrsz(4)/20 widin*res htin*res]);
ax = gca;
hold(ax,'on');
patarea = [0 0;
           -100 0;
           -100 1;
            0 1;
            0 0];
patch(ax,patarea(:,1),patarea(:,2),'k','Facealpha',0.05,'edgecolor','none');
plot(ax,[meansw meansw],[-100 100],'r--','linew',2);    % median of wt dist
plot(ax,[meanne meanne],[-100 100],'b--','linew',2);    % median of wt dist
errorbar(ax,meansw,0.25,CIsw(1)-meansw,CIsw(2)-meansw,'horizontal','o','markersize',6,'color',...
         'r','linewidth',1.5,'MarkerEdgeColor','r','MarkerFaceColor','r','CapSize',10);
errorbar(ax,meanne,0.25,CIne(1)-meanne,CIne(2)-meanne,'horizontal','o','markersize',6,'color',...
         'b','linewidth',1.5,'MarkerEdgeColor','b','MarkerFaceColor','b','CapSize',10);
text(ax,0.58,0.94,strcat({'10X zoom-in'}),'fontsize',12,'unit','normalized','EdgeColor','k',...
     'Margin',2);
xlimit = f.ax(2).XLim/10;
xlim(ax,xlimit); 
ylim(ax,[0.2 0.3]);
% xticks(xlimit(1): 0.2: xlimit(2));
yticks(ax,[0.2 0.25 0.3]);
xlabel(ax,'WSW-ENE (km)','fontsize',11);
ylabel(ax,'PDF estimate','fontsize',11);
ax.Box = 'on';
grid(ax, 'on');
ax.GridLineStyle = '--';
hold(ax,'off');


%%% 2. for SW and NE group at fam 002 region
f.fig=figure;
widin = 8;  % maximum width allowed is 8.5 inches
htin = 4.5;   % maximum height allowed is 11 inches
set(f.fig,'Position',[scrsz(1)+1*scrsz(3)/10 scrsz(2)+scrsz(4)/20 widin*res htin*res]);
nrow = 1;
ncol = 2;
for isub = 1:nrow*ncol
    f.ax(isub) = subplot(nrow,ncol,isub);
end

%%% reposition
set(f.ax(1), 'position', [ 0.1, 0.1, 0.35, 0.75]);

% subplot 1
hold(f.ax(1),'on');
xran = [-5 20];
yran = [-10 10];
indplt = union(ind002sw,ind002ne);    % index to plot
[f.ax(1),c] = plt_mig_prop_onmap(f.ax(1), indplt, ranvechf98, medallhf, angbest, xran,yran);
hold(f.ax(1),'off');

% subplot 2 
ywid = f.ax(1).Position(4)+f.ax(1).Position(2)-(c.Position(2)+0.014);
set(f.ax(2), 'position', [ 0.55, c.Position(2)-0.00, 0.35, ywid]);

mean002sw = wt_mean(vdist002sw,wt002sw);
mean002ne = wt_mean(vdist002ne,wt002ne);
sigma002sw = sqrt(wt_var(vdist002sw,wt002sw,2));
neff = sum(wt002sw);
% conf = 99;
CI002sw = confidence_interval(mean002sw,sigma002sw,neff,conf);
sigma002ne = sqrt(wt_var(vdist002ne,wt002ne,2));
neff = sum(wt002ne);
% conf = 99;
CI002ne = confidence_interval(mean002ne,sigma002ne,neff,conf);

hold(f.ax(2),'on');
[f.ax(2), barsw002, pdfxloc002sw, ~, ~, pdfval002sw] = plt_weighted_dist(f.ax(2), vdist002sw, wt002sw, ...
                                                                      binw,'dec');
barsw002(1).FaceAlpha = 1;
[f.ax(2), barne002, pdfxloc002ne, ~, ~, pdfval002ne] = plt_weighted_dist(f.ax(2), vdist002ne, wt002ne, ...
                                                                      binw,'dec');
barne002(1).FaceColor = [0 1 1];
plot(f.ax(2),[mean002sw mean002sw],[-100 100],'r--','linew',1.5);    % median of wt dist
plot(f.ax(2),[mean002ne mean002ne],[-100 100],'b--','linew',1.5);    % median of wt dist
errorbar(f.ax(2),mean002sw,0.25,CI002sw(1)-mean002sw,CI002sw(2)-mean002sw,'horizontal','o',...
         'markersize',3,'color','r','linewidth',1,'MarkerEdgeColor','r','MarkerFaceColor','r');
errorbar(f.ax(2),mean002ne,0.25,CI002ne(1)-mean002ne,CI002ne(2)-mean002ne,'horizontal','o',...
         'markersize',3,'color','b','linewidth',1,'MarkerEdgeColor','b','MarkerFaceColor','b');
% text(f.ax(2),0.5,0.9,strcat({'med = '},sprintf('%.2f',median_of_wt_data(vdist002sftswne,...
%      wt002sftswne)),{' km'}),'fontsize',12,'unit','normalized');
xlabel(f.ax(2),'SW-NE (km)','fontsize',11);
legend(f.ax(2),[barsw002,barne002],{'SW group','NE group'},'fontsize',8);
hold(f.ax(2),'off');

%%% save figure
print(f.fig,'-dpdf',strcat(rstpath,'/vertical_dist_unshifted_SW+NE002_migs',num2str(nfam),'.pdf'));

ff.fig=figure;
widin = 3.5;  % maximum width allowed is 8.5 inches
htin = 3.8;   % maximum height allowed is 11 inches
set(ff.fig,'Position',[scrsz(1)+1*scrsz(3)/10 scrsz(2)+scrsz(4)/20 widin*res htin*res]);
ax = gca;
hold(ax,'on');
patarea = [0 0;
           -100 0;
           -100 1;
            0 1;
            0 0];
patch(ax,patarea(:,1),patarea(:,2),'k','Facealpha',0.05,'edgecolor','none');
plot(ax,[mean002sw mean002sw],[-100 100],'r--','linew',2);    % median of wt dist
plot(ax,[mean002ne mean002ne],[-100 100],'b--','linew',2);    % median of wt dist
errorbar(ax,mean002sw,0.25,CI002sw(1)-mean002sw,CI002sw(2)-mean002sw,'horizontal','o','markersize',...
         6,'color','r','linewidth',1.5,'MarkerEdgeColor','r','MarkerFaceColor','r','CapSize',10);
errorbar(ax,mean002ne,0.25,CI002ne(1)-mean002ne,CI002ne(2)-mean002ne,'horizontal','o','markersize',...
         6,'color','b','linewidth',1.5,'MarkerEdgeColor','b','MarkerFaceColor','b','CapSize',10);
text(ax,0.6,0.94,strcat({'5X zoom-in'}),'fontsize',12,'unit','normalized','EdgeColor','k',...
     'Margin',2);
xlimit = f.ax(2).XLim/5;
xlim(ax,xlimit); 
ylim(ax,[0.2 0.3]);
% xticks(-2: 0.2: 1);
yticks(ax,[0.2 0.25 0.3]);
xlabel(ax,'SW-NE (km)','fontsize',11);
ylabel(ax,'PDF estimate','fontsize',11);
ax.Box = 'on';
grid(ax, 'on');
ax.GridLineStyle = '--';
hold(ax,'off');


%%
% NE+SW migs that pass through similar region
temp = union(indsw,indne);
indswne = union(temp,[22,25]);
vdistswne = [];
wtswne = [];
ivdistswne = [];  % dist and weight of independent detections
iwtswne = [];
for i = 1: length(indswne)
    iind = indlfindpen(1:inumlf(indswne(i)), indswne(i));  % index of independent detection in that migration indsw(i)
    if indswne(i) == 9
        temp1 = vertdist(1:numlf(indswne(i)), indswne(i));
        temp2 = temp1(ind9);
        vdistswne = [vdistswne; temp2];
        temp2 = temp1(intersect(iind, ind9));
        ivdistswne = [ivdistswne; temp2];
        temp1 = weight(1:numlf(indswne(i)), indswne(i));
        temp2 = temp1(ind9);
        wtswne = [wtswne; temp2];
        temp2 = temp1(intersect(iind, ind9));
        iwtswne = [iwtswne; temp2];
    else
        temp1 = vertdist(1:numlf(indswne(i)), indswne(i));
        vdistswne = [vdistswne; temp1];
        ivdistswne = [ivdistswne; temp1(iind)];
        temp1 = weight(1:numlf(indswne(i)), indswne(i));
        wtswne = [wtswne; temp1];
        iwtswne = [iwtswne; temp1(iind)];
    end
end
pdswne = fitdist(vdistswne,'Normal');  
muswne = pdswne.mu;    % fitted paramaters
sigmaswne = pdswne.sigma;
pdffitswne = pdf(pdswne,-50:0.05:50);

% 2 NE and SW migs that pass through fam 002 region
ind002 = [22,25];
vdist002 = [];
wt002 = [];
ivdist002 = [];  % dist and weight of independent detections
iwt002 = [];
for i = 1: length(ind002)
    iind = indlfindpen(1:inumlf(ind002(i)), ind002(i));  % index of independent detection in that migration indsw(i)
    temp1 = vertdist(1:numlf(ind002(i)), ind002(i));
    vdist002 = [vdist002; temp1];
    ivdist002 = [ivdist002; temp1(iind)];
    temp1 = weight(1:numlf(ind002(i)), ind002(i));
    wt002 = [wt002; temp1];
    iwt002 = [iwt002; temp1(iind)];

end
pd002 = fitdist(vdist002,'Normal');  
mu002 = pd002.mu;    % fitted paramaters
sig002 = pd002.sigma;
pdffit002 = pdf(pd002,-50:0.05:50);


% all migs excluding 
indall = 1: 1: size(trange,1);
indelse = setdiff(indall, indswne);
vdistelse = [];
wtelse = [];
ivdistelse = [];  % dist and weight of independent detections
iwtelse = [];
for i = 1: length(indelse)
    iind = indlfindpen(1:inumlf(indelse(i)), indelse(i));  % index of independent detection in that migration indsw(i)
    if indelse(i) == 4
        temp1 = vertdist(1:numlf(indelse(i)), indelse(i));
        temp2 = temp1(ind4);
        vdistelse = [vdistelse; temp2];
        temp2 = temp1(intersect(iind, ind4));
        ivdistelse = [ivdistelse; temp2];
        temp1 = weight(1:numlf(indelse(i)), indelse(i));
        temp2 = temp1(ind4);
        wtelse = [wtelse; temp2];
        temp2 = temp1(intersect(iind, ind4));
        iwtelse = [iwtelse; temp2];
    elseif indelse(i) == 9
        temp1 = vertdist(1:numlf(indelse(i)), indelse(i));
        temp2 = temp1(ind9);
        vdistelse = [vdistelse; temp2];
        temp2 = temp1(intersect(iind, ind9));
        ivdistelse = [ivdistelse; temp2];
        temp1 = weight(1:numlf(indelse(i)), indelse(i));
        temp2 = temp1(ind9);
        wtelse = [wtelse; temp2];
        temp2 = temp1(intersect(iind, ind9));
        iwtelse = [iwtelse; temp2];
    else
        temp1 = vertdist(1:numlf(indelse(i)), indelse(i));
        vdistelse = [vdistelse; temp1];
        ivdistelse = [ivdistelse; temp1(iind)];
        temp1 = weight(1:numlf(indelse(i)), indelse(i));
        wtelse = [wtelse; temp1];
        iwtelse = [iwtelse; temp1(iind)];
    end
end
pdelse = fitdist(vdistelse,'Normal');  
muelse = pdelse.mu;    % fitted paramaters
sigmaelse = pdelse.sigma;
pdffitelse = pdf(pdelse,-50:0.05:50);


%% divide each group spatially into 2 parts according to a pre-designed line
%%% 1. divide SW group at main LZB region into NE and SW parts
indplt = indsw;    % index to plot
xran = [-15 5];
yran = [-15 15];
binw = 1;
conf = 99;
[f,barsw,pdfxlocsw,pdfvalsw,barne,pdfxlocne,pdfvalne,lgd] = ...
    plt_hist_combined_RTMs(indplt,ranvechf98,medallhf,angbest,xran,yran,binw,vdistswlt,wtswlt,...
                           vdistswgt,wtswgt,conf);
hold(f.ax(1),'on');
plot(f.ax(1),xdiv,ydiv);
text(f.ax(1),0.5,0.06,'All detections, SW propagation','FontSize',12,'unit','normalized',...
     'horizontalalignment','center');
text(f.ax(1),0.5,0.12,'Main LZB region,','FontSize',12,'unit','normalized',...
     'horizontalalignment','center');
text(f.ax(2),0.05,0.1,strcat({'WSW'}),'fontsize',11,'unit','normalized');
text(f.ax(2),0.95,0.1,strcat({'ENE'}),'fontsize',11,'unit','normalized',...
     'horizontalalignment','right');
lgd.String = [{'WSW part'},{'ENE part'}];
text(f.ax(3),0.95,0.86,strcat({'10X'}),'fontsize',10,'unit','normalized','EdgeColor','k',...
     'Margin',2,'horizontalalignment','right');
text(f.ax(3),0.05,0.15,strcat({'99% CI'}),'fontsize',12,'unit','normalized','horizontalalignment','left'); 

xlimit = f.ax(2).XLim/10;
xlim(f.ax(3),xlimit); 
% convert data units to global figure units
[xs,ys] = ds2nfu(f.ax(3),xlimit(1),0.2);
[xe,ye] = ds2nfu(f.ax(2),xlimit(1),0.3);
% plot two dashed lines denoting the zoom-in effect
annotation('line',[xs,xe],[ys,ye],'color','k','linestyle','--');
% convert data units to global figure units
[xs,ys] = ds2nfu(f.ax(3),xlimit(2),0.2);
[xe,ye] = ds2nfu(f.ax(2),xlimit(2),0.3);
% plot two dashed lines denoting the zoom-in effect
annotation('line',[xs,xe],[ys,ye],'color','k','linestyle','--');

%%% save figure
print(f.fig,'-dpdf',strcat(rstpath,'/vdist_unshifted_mainSW_all',num2str(nfam),'.pdf'));


%%% 2. divide NE group at main LZB region into NE and SW parts
indplt = indne;    % index to plot
[f,barsw,pdfxlocsw,pdfvalsw,barne,pdfxlocne,pdfvalne,lgd] = ...
    plt_hist_combined_RTMs(indplt,ranvechf98,medallhf,angbest,xran,yran,binw,vdistnegt,wtnegt,...
                           vdistnelt,wtnelt,conf);
hold(f.ax(1),'on');
plot(f.ax(1),xdiv,ydiv);
text(f.ax(1),0.5,0.06,'All detections, NE propagation','FontSize',12,'unit','normalized',...
     'horizontalalignment','center');
text(f.ax(1),0.5,0.12,'Main LZB region,','FontSize',12,'unit','normalized',...
     'horizontalalignment','center');
text(f.ax(2),0.05,0.1,strcat({'WSW'}),'fontsize',11,'unit','normalized');
text(f.ax(2),0.95,0.1,strcat({'ENE'}),'fontsize',11,'unit','normalized',...
     'horizontalalignment','right');
lgd.String = [{'WSW portion'},{'ENE portion'}];
text(f.ax(3),0.95,0.86,strcat({'10X'}),'fontsize',10,'unit','normalized','EdgeColor','k',...
     'Margin',2,'horizontalalignment','right');
text(f.ax(3),0.05,0.15,strcat({'99% CI'}),'fontsize',12,'unit','normalized','horizontalalignment','left'); 

xlimit = f.ax(2).XLim/10;
xlim(f.ax(3),xlimit); 
% convert data units to global figure units
[xs,ys] = ds2nfu(f.ax(3),xlimit(1),0.2);
[xe,ye] = ds2nfu(f.ax(2),xlimit(1),0.3);
% plot two dashed lines denoting the zoom-in effect
annotation('line',[xs,xe],[ys,ye],'color','k','linestyle','--');
% convert data units to global figure units
[xs,ys] = ds2nfu(f.ax(3),xlimit(2),0.2);
[xe,ye] = ds2nfu(f.ax(2),xlimit(2),0.3);
% plot two dashed lines denoting the zoom-in effect
annotation('line',[xs,xe],[ys,ye],'color','k','linestyle','--');

%%% save figure
print(f.fig,'-dpdf',strcat(rstpath,'/vdist_unshifted_mainNE_all',num2str(nfam),'.pdf'));



%% summary of results with error bars from lf fitting
CIFcn = @(x,p)std(x(:),'omitnan')/sqrt(sum(~isnan(x(:)))) * tinv(abs([0,1]-(1-p/100)/2),...
              sum(~isnan(x(:)))-1) + mean(x(:),'omitnan');
f5.fig=figure;
set(f5.fig,'Position',[scrsz(3)/10 2*scrsz(4)/10 2*scrsz(3)/6 0.7*scrsz(4)]);
hold on
plot([0 0],[-100 100],'--','color',[0.5 0.5 0.5]);
for i = 1: length(trange)
    hfN = numhf(i);
    hf=resprophf(i,1:hfN);
    hfm = mean(hf);
    hfstd = std(hf);
    hfsem = hfstd/sqrt(hfN);
    ci95 = tinv([0.025 0.975], hfN-1);
    hfci95n = bsxfun(@times, hfsem, ci95(:));
    hfci95 = hfm+hfci95n;
    hfCI95 = CIFcn(hf,95);
    
    lfN = numlf(i);
    lf=resproplf(i,1:lfN);
    lfm = mean(lf);
    lfstd = std(lf);
    lfsem = lfstd/sqrt(lfN);
    ci95 = tinv([0.025 0.975], lfN-1);
    lfci95n = bsxfun(@times, lfsem, ci95(:));
    lfci95 = lfm+lfci95n;
    lfCI95 = CIFcn(lf,95);
    
    offset(i) = (intcptproplf(i)-intcptprophf(i));
    ypos = i;
    if angbest(i)>0 && angbest(i)<=90
        e1 = errorbar(offset(i),ypos,lfci95n(1),lfci95n(2),'horizontal','o','markersize',6,'color',...
         'k','linewidth',1,'MarkerEdgeColor','k','MarkerFaceColor',[0 1 1]);
    elseif angbest(i)>90 && angbest(i)<=180
        e2 = errorbar(offset(i),ypos,lfci95n(1),lfci95n(2),'horizontal','o','markersize',6,'color',...
         'k','linewidth',1,'MarkerEdgeColor','k','MarkerFaceColor','k');
    elseif angbest(i)>180 && angbest(i)<=270
        e3 = errorbar(offset(i),ypos,lfci95n(1),lfci95n(2),'horizontal','o','markersize',6,'color',...
         'k','linewidth',1,'MarkerEdgeColor','k','MarkerFaceColor',[1 96/255 0]);
    end
    text(4.8,ypos,strcat(num2str(hfN),'; ',num2str(lfN)),'fontsize',9,...
        'HorizontalAlignment','right');
    %%% next time use 'units', 'normalized' instead of defaults; also 
end

text(0.02,0.97,'a','FontSize',12,'unit','normalized','EdgeColor','k','Margin',2);
text(0.9,0.97,'LZB','FontSize',14,'unit','normalized','EdgeColor','k','Margin',2);
text(0.2,0.13,'LF','FontSize',14,'unit','normalized','HorizontalAlignment','center','fontw','bold');
text(0.2,0.1,'lags','FontSize',14,'unit','normalized','HorizontalAlignment','center','fontw','bold');
text(0.8,0.13,'LF','FontSize',14,'unit','normalized','HorizontalAlignment','center','fontw','bold');
text(0.8,0.1,'leads','FontSize',14,'unit','normalized','HorizontalAlignment','center','fontw','bold');
xlabel('Offset between HF and LF (LF-HF) (km)','fontsize',12);
ylabel('Migration number from earlier to later','fontsize',12);
legend([e1, e2, e3],{'NE propagation','SE propagation','SW propagation'},'FontSize',10,...
        'Position',[0.17 0.52 0.2 0.08]);
ylim([0,50]);
xlim([-5,5]);
% yticks(0:1:8);
box on
grid on
set(gca,'GridLineStyle','--')

%%% save figure
print(f5.fig,'-dpdf',strcat(rstpath,'/LZB',num2str(nfam),'.mig.lfit.sum.',num2str(winlenhf),'_',...
    num2str(winlenlf),'.',num2str(lofflf),'.',num2str(ccminlf),'.pdf'));

