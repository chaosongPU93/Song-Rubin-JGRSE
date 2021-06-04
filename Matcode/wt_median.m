function med_wt = wt_median(x, wt)
% med_wt = wt_median(x, wt)
% This function is to get the median of the data that has different weights for
% each point. Trying to deal with it from a point view that if the data lies on
% x axis, median is the value on X axis where the sum of weights of all points
% on its left equals to the sum of weights of all points on its right.
%
% HAVE been tested to be correct, compared with median(dist) with equal weights
%
%
% By Chao Song, chaosong@princeton.edu
% First created date:   2020/10/05
% Last modified date:   2020/10/08

% reshape to one column
x = reshape(x, [],1);
wt = reshape(wt, [],1);
data = [x wt];

% sort it 
datasort = sortrows(data,1);

% cumulative sum of weights from smallest to the larger
wtle = cumsum(datasort(:,2));   % sum of the left, <=
wtgt = sum(datasort(:,2))- wtle;    % sum of the right, >
wtdif = wtgt - wtle;    % difference

% find the separation, usually we won't get a exact index where 2 sides sum to the same amount 
ind = find(abs(wtdif) == min(abs(wtdif)));      % find the index which has the least abs diff

if length(ind) == 1     % zero-crossing is somewhere between this index and one before
    dist1 = datasort(ind+1,1);
    wt1 = datasort(ind+1,2);
    dist2 = datasort(ind,1);
    wt2 = datasort(ind,2);
    med_wt = (dist1*wt1+dist2*wt2)/(wt1+wt2);
elseif length(ind) == 2     % zero-crossing exactly at the 2nd index
    med_wt = datasort(ind(2),1);
end    
    