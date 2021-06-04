function [meanrand,pdfmisfit,exppdf] = randgen_custom(seed,N,vdist,wt,pdfxloc,pdfval)
% [randnum,muHat,sigmaHat,muCI,sigmaCI] = RANDGEN_CUSTOM(data, wt, dimsize)
%
% This function is to generate random numbers from the data 'vdist' with
% non-equal weights 'wt' directly using 'ksdensity' and its custom PDF 
% 'pdfxloc,pdfval', by N times. Take in random seed state for reproduction.
% Output the median of randoms and the misfit between randoms PDF and data
% PDF in each round of generation. 
%
% By Chao Song, chaosong@princeton.edu
% First created date:   2020/10/01
% Last modified date:   2020/10/01


meanrand = zeros(N, 3);    % mean of the generated random numbers
pdfmisfit = zeros(N, 3);    % misfits of PDFs between data and generated random numbers
dimsize = [ceil(sum(wt)),1];  % size of desired randum numbers per sampling
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
    [randnum, exppdf] = randpdf(pdfxloc, pdfval, dimsize);
    rdnum = [rdnum randnum];

    % 2. use ksdensity
    rng(seed{i});
    randnum = ksdensity(vdist, rand(dimsize), 'function', 'icdf', 'weight', wt);
    rdnum = [rdnum randnum];
    
    % 3. assume pdf is normal-like, and try to fit it with a gaussian 
    rng(seed{i});
    [randnum,muHat,sigmaHat,~,~] = randnorm_wt(vdist, wt, dimsize);
    rdnum = [rdnum randnum];
%     some sense of the misfit
%     [pdfmisfit1, pdfmisfit2, f66] = compare_distribution2(randnum,muHat,sigmaHat,vdist,wt,pdfxloc,pdfval);
    
%     %%% check the distribution of generated random numbers to see if it accords with the pdf of data
%     pdfmisf = zeros(size( rdnum,2), 1);
%     titlestr = {'rand gen from custom pdf';
%                 'rand gen from ksdensity';
%                 'rand gen from normal fitting';
%            };
%     for ii = 1: size(rdnum,2)
%         binw = pdfxloc(2)-pdfxloc(1);
%         [pdfmisf(ii), f55] = compare_distribution1(rdnum(:,ii),vdist,wt,pdfxloc,pdfval,exppdf,binw);
%         supertit(f55.ax, string(titlestr{ii}),12);
%         
%     end
    
    meanrand(i,:) = mean(rdnum,1);
    
    %%% output the misfits of PDFs as well but no figures
    for ii = 1: size(rdnum,2)
        binw = pdfxloc(2)-pdfxloc(1);
        pdfmisfit(i, ii) = check_randgen(rdnum(:,ii),vdist,wt,pdfxloc,pdfval,exppdf,binw);
    end
    
end

% keyboard









