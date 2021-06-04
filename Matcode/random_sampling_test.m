function [pdiff,mdiff,refval] = random_sampling_test(N,x1,wt1,pdfx1,pdf1,x2,wt2,pdfx2,pdf2)
% [pdiff,p1,p2,f1,f2,f3] = random_sampling_test(N,x1,wt1,pdfx1,pdf1,x2,wt2,pdfx2,pdf2)
%
% This function is to bundle the corresponding part of random sampling test
% between 2 data sets with weights using 3 different methods if their
% difference in mean are significant enough
%
%
% By Chao Song, chaosong@princeton.edu
% First created date:   2020/10/22
% Last modified date:   2020/10/22
                          
% NE group
% set and save random seed state for reproduction 
seed1 = cell(1,N);
for i = 1: N
    seed1{i} = rng('shuffle','twister');
end
mx1 = wt_mean(x1,wt1);
% mrandnorm is tracking the weighted mean of the data
[~,mrandnorm1,~] = randgen_norm(seed1,N,x1,wt1,pdfx1,pdf1);
% mrandkd is tracking the weighted mean of the data
[~,mrandkd1,~] = randgen_ksdensity(seed1,N,x1,wt1,pdfx1,pdf1);
% mrandpdf is tracking the expectation of the PDF of data
[~,mrandpdf1,~,exppdf1] = randgen_pdf(seed1,N,x1,wt1,pdfx1,pdf1);

mrand1 = [mrandnorm1 mrandkd1 mrandpdf1];

% SW group
seed2 = cell(1,N);
for i = 1: N
    seed2{i} = rng('shuffle','twister');
end
mx2 = wt_mean(x2,wt2);
% mrandnorm is tracking the weighted mean of the data
[~,mrandnorm2,~] = randgen_norm(seed2,N,x2,wt2,pdfx2,pdf2);
% mrandkd is tracking the weighted mean of the data
[~,mrandkd2,~] = randgen_ksdensity(seed2,N,x2,wt2,pdfx2,pdf2);
% mrandpdf is tracking the expectation of the PDF of data
[~,mrandpdf2,~,exppdf2] = randgen_pdf(seed2,N,x2,wt2,pdfx2,pdf2);

mrand2 = [mrandnorm2 mrandkd2 mrandpdf2];

% difference of mean from random numbers generated, and from data
% mrandnorm is tracking the weighted mean of the data
% mrandkd is tracking the weighted mean of the data
% mrandpdf is tracking the expectation of the PDF of data
mdiff = [(mrand1(:,1:2)-mx1)-(mrand2(:,1:2)-mx2) (mrand1(:,3)-exppdf1)-(mrand2(:,3)-exppdf2)];

% same for the reference value, 
refval = [mx1-mx2 mx1-mx2 exppdf1-exppdf2];

% probability that this difference could be more extreme than the data
pdiff = sum(abs(mdiff)>=abs(refval), 1) ./N ;
% for i = 1: size(mdiff,2)
%     pdiff(i) = sum(abs(mdiff(:,i))>=abs(refval(:,i)), 1) ./N .* 100;
% end


% keyboard










