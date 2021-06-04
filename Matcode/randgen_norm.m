function [randnum,meanrand,pdfmisfit] = randgen_norm(seed,N,vdist,wt,pdfxloc,pdfval)
% [randnum,muHat,sigmaHat,muCI,sigmaCI] = RANDGEN_CUSTOM(data, wt, dimsize)
%
% This function is to generate random numbers from a gaussian fit to the data
% 'vdist' with non-equal weights 'wt'. Take in
% random seed state for reproduction.
% Output the median of randoms and the misfit between randoms PDF and data
% PDF in each round of generation. 
%
% By Chao Song, chaosong@princeton.edu
% First created date:   2020/10/27
% Last modified date:   2020/10/27

% meanrand = zeros(N, 3);    % mean of the generated random numbers
pdfmisfit = zeros(N, 3);    % misfits of PDFs between data and generated random numbers
dimsize = [ceil(sum(wt)),1];  % size of desired randum numbers per sampling
randnum = zeros(ceil(sum(wt)), N);
for i = 1: N
% i = randi(1000,1);

    % 2. use ksdensity
    rng(seed{i});
    [randnum(:,i),muHat,sigmaHat,~,~] = randnorm_wt(vdist, wt, dimsize);
%     some sense of the misfit
%     [pdfmisfit1, pdfmisfit2, f66] = compare_distribution2(randnum,muHat,sigmaHat,vdist,wt,pdfxloc,pdfval);    
    
    %%% output the misfits of PDFs as well but no figures
    for ii = 1: size(randnum(:,i),2)
        binw = pdfxloc(2)-pdfxloc(1);
        refval = wt_mean(vdist, wt);
        pdfmisfit(i) = check_randgen(randnum(:,i),vdist,wt,pdfxloc,pdfval,refval,binw);
    end
    
end

meanrand = mean(randnum,1);   % mean of the generated random numbers
meanrand = meanrand';

% keyboard