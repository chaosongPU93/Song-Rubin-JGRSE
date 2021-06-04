function [randnum,meanrand,pdfmisfit,exppdf] = randgen_pdf(seed,N,vdist,wt,pdfxloc,pdfval)
% [meanrand,pdfmisfit,exppdf] = randgen_pdf(seed,N,pdfxloc,pdfval)
%
% This function is to generate random numbers from the custom PDF from data
% 'pdfxloc,pdfval', by N times. Take in random seed state for reproduction.
% Output the median of randoms and the misfit between randoms PDF and data
% PDF in each round of generation. 
%
% By Chao Song, chaosong@princeton.edu
% First created date:   2020/10/27
% Last modified date:   2020/10/27


% meanrand = zeros(N, 1);    % mean of the generated random numbers
pdfmisfit = zeros(N, 1);    % misfits of PDFs between data and generated random numbers
dimsize = [ceil(sum(wt)),1];  % size of desired randum numbers per sampling
randnum = zeros(ceil(sum(wt)), N);
for i = 1: N
% i = randi(1000,1);

    rdnum = [];
    
    % 1. use randpdf
    %%% NOTE:
    %%% the shifting arises bc. the obtained pdf requires certain binning, so if you want to shift
    %%% the data by amount 'shift' by directly add it to the data, it is fine; but if the bin width
    %%% and edge remains unchanged, the histogram (pdf) usually would NOT have a same shape as the
    %%% unshifted histogram (pdf), so preserve the shape, you MUST add the shift to the pdf x
    %%% location, rather than do binning again to the shifted new data
    rng(seed{i});
    [randnum(:,i), exppdf] = randpdf(pdfxloc, pdfval, dimsize);
    
    
    %%% output the misfits of PDFs as well but no figures
    for ii = 1: size(randnum(:,i),2)
        binw = pdfxloc(2)-pdfxloc(1);
        refval = exppdf;
        pdfmisfit(i) = check_randgen(randnum(:,i),vdist,wt,pdfxloc,pdfval,refval,binw);
    end
    
end

meanrand = mean(randnum,1);   % mean of the generated random numbers
meanrand = meanrand';

% keyboard
