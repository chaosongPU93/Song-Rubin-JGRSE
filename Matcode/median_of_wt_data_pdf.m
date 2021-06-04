function meddist = median_of_wt_data_pdf(pdfx, pdf)
% 
% This function is to get the median of the data that has different weights for
% each point. Trying to deal with it from a point view that if you have the PDF
% of the data with a certain distribution, you could obtain the CDF. Then the 
% median is just the point where the probability at this point is exactly 0.5
%
% the built-in empirical cdf function 'ecdf' is not suitable
% for this case because the input needs to be counts, not the
% the pdf
%
% By Chao Song, chaosong@princeton.edu
% First created date:   2020/10/05
% Last modified date:   2020/10/05

% vectorization and normalization of the input pdf, just in case
pdfx=pdfx(:);
pdf=pdf(:)./trapz(pdfx,pdf(:));

% interpolation of the input pdf for better integration
% in my opinion 1000 point is sufficient...
pdfxi=linspace(min(pdfx),max(pdfx),1000);
pdfxi = pdfxi';
pdfi=interp1(pdfx,pdf,pdfxi,'linear');

% computing the cumulative distribution function for input pdf
cdfi = cumtrapz(pdfxi,pdfi);

% finding the parts of cdf parallel to the X axis 
ind=[true; not(diff(cdfi)==0)];

% and cut out the parts
cdfi=cdfi(ind);
pdfi=pdfi(ind);
pdfxi=pdfxi(ind);

% interpolation at the 0.5 probability to the cdf
meddisti = interp1(cdfi,pdfxi,0.5,'linear');

%%% if without interpolation
cdf = cumtrapz(pdfx,pdf);
ind=[true; not(diff(cdf)==0)];
cdf=cdf(ind);
pdf=pdf(ind);
pdfx=pdfx(ind);

meddist = interp1(cdf,pdfx,0.5,'linear');






















