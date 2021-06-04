function [xloc, ctval, probval, pdfval] = weighted_bar_pdf(vdist, wt, binwid, edgetype)
% This function is to compute the x loc and value for the bar plot
% normalized as the PDF estimate for the data that has different
% weight for each point. Otherwise this could be done easily with
% 'histogram' using 'normalization' 'pdf' option
%
%   Served more as a simplification in the main code  
%
% By Chao Song, chaosong@princeton.edu
% First created date:   2020/10/01
% Last modified date:   2020/10/01

% trying to use the auto bin method histogram to have a sense of the reasonable bin width
ff = figure;
tmp = vdist(wt~=0);
h=histogram(tmp,'binmethod','auto','facec',[1 96/255 0],'facea',0.8);
      
% prob = sum(length(h.Data).*h.BinWidth.*h.Values)
% binw = 0.2;
if isempty(binwid)
    binw = h.BinWidth;
else
    binw = binwid;
end
% delete(h);
% edgest = -50;
% edgeed = 50;

if strcmp(edgetype, 'int')
    edgest = floor(min(vdist));
    edgeed = ceil(max(vdist));
    xloc = edgest+binw/2: binw: edgeed+binw/2;
elseif strcmp(edgetype, 'dec')
    nbinl = ceil((ceil(abs(min(vdist)/binw*2))-1)/2);
    nbinr = ceil((ceil(abs(max(vdist)/binw*2))-1)/2);
    xloc = -nbinl*binw: binw: nbinr*binw;
end
    
xloc = xloc';
h = zeros(length(xloc),1);
h1 = zeros(length(xloc),1);
h2 = zeros(length(xloc),1);
for i = 1: length(xloc)
%     ibar = find(vdist>xloc(i)-binw/2 & vdist<=xloc(i)+binw/2);
    wtobj = wt(vdist>xloc(i)-binw/2 & vdist<=xloc(i)+binw/2);
    h(i) = sum(wtobj);  % count
    h1(i) = sum(wtobj)/sum(wt);    % probability
    h2(i) = sum(wtobj)/sum(wt)/binw;  % probobility density estimate, pdf
    
%     h(i) = length(ibar)/length(vdistsw);    % probability
%     h2(i) = h(i)/binw;  % probobility density estimate    

end
indst = find(h2>0, 1, 'first');     % non-zero start index
inded = find(h2>0, 1, 'last');      % non-zero end index

% get the x location, non-zero count, probability and pdf and its corresponding
xloc = xloc(indst: inded);
ctval = h(indst: inded);
probval = h1(indst: inded);
pdfval = h2(indst: inded);

close(ff);