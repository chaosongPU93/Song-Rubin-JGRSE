function [uniqrst, uniqtime] = ...
    CentroidMerge(rstold, rstnew, timeold, timenew) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function is used to concatenate two matrix by row, unique according
% to first 2 columns (presuming they indicate coordinates), use the average
% as the value of remaining columns
% 
%
% Modified by Chao Song, chaosong@princeton.edu
% First created date:   2019/05/22
% Last modified date:   2019/06/18
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%% NOT USED ANYMORE, but worthwhile to save for future reference %%%%%%
% rstraw = [rstold; rstnew];  % concatenate by row
% rstsort = sortrows(rstraw, [1,2]);  % sort by column 1 & 2
% col1uni = unique(rstsort(:,1)); % the unique 3rd col
% uniqrst = [];
% % hitdiftol = 5;
% for i = 1: length(col1uni)
%     temp1 = rstsort(rstsort(:,1) == col1uni(i), :);
%     col2uni = unique(temp1(:,2));
%     for j = 1: length(col2uni)
%         temp2 = temp1(temp1(:,2) == col2uni(j), :);
%         temp3 = zeros(1, size(temp2,2));
%         if size(temp2,1) > 1
%             if max(temp2(:,end))-min(temp2(:,end)) <= hitdiftol
%                 temp3 = sum(temp2, 1)/size(temp2,1);
%             else
%                 temp3 = temp2(temp2(:, end) == max(temp2(:,end)), :);
%             end
%             uniqrst = [uniqrst; temp3];
%         else
%             uniqrst = [uniqrst; temp2];
%         end
%     end
%         
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% set default value, easy for debug
defval('rstold', 'rstold23')
defval('rstnew', 'rstnew23')
defval('timeold','timeold23')
defval('timenew','timenew23')

[~, iold, inew] = intersect(rstold(:,1:2), rstnew(:,1:2), 'row', 'stable'); % get the intersection of locations
dum1 = rstold;
dum2 = rstnew;
dum1(iold,:)=[];    
dum2(inew,:)=[];
rstraw = [dum1;dum2];   % first merge the exclusive parts
flag = zeros(length(iold),2);   % flag to indicate which matrix is used
for i = 1: length(iold)
    if rstold(iold(i),end) >= rstnew(inew(i),end)   % only compare the hit counts
        flag(i,1) = 1;
        temp = rstold(iold(i),:);
    else
        flag(i,2) = 1;
        temp = rstnew(inew(i),:);
    end
    rstraw = [rstraw; temp];
end

timeraw=[];
for i = 1: size(dum1,1)
% i=3
    [~, ind] = ismember(dum1(i,1:2),timeold(:,1:2), 'rows');    % first get the start index 
    temp = timeold(ind:ind+dum1(i,end)-1, :);   % extract the number of rows according to hit count
    timeraw = [timeraw; temp];  
end
for i = 1: size(dum2,1)
    [~, ind] = ismember(dum2(i,1:2),timenew(:,1:2), 'rows');
    temp = timenew(ind:ind+dum2(i,end)-1, :);
    timeraw = [timeraw; temp];
end

for i = 1: length(iold)
    if flag(i,1)    % merge the intersection part according the previous flag
        [~, ind] = ismember(rstold(iold(i),1:2),timeold(:,1:2), 'rows');    
        temp = timeold(ind:ind+rstold(iold(i),end)-1, :);

    else
        [~, ind] = ismember(rstnew(inew(i),1:2),timenew(:,1:2), 'rows');
        temp = timenew(ind:ind+rstnew(inew(i),end)-1, :);
    end
    timeraw = [timeraw; temp];
end
uniqtime = sortrows(timeraw, [1, 2, 5, 6, 7]);  % sort rows according to locaton, date, begin t, center t
uniqrst = sortrows(rstraw, [1,2]);  % sort rows according to locaton




