%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% THIS is a callable Function to return the 3-D offset and count maxtrix
% of checking with additional station. 
%
%
% 
%
% Modified by Chao Song, chaosong@princeton.edu
% First created date:   2019/04/17
% Last modified date:   2019/04/17
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [temp_count14, temp_off14, temp_tct14, temp_tbg14, temp_tdate14, temp_cc14] = ...
    Addcheck_offcount(istanew,nstanew,date,mshift,allrst_new,countmat14,offmat14,tctmat14,...
                      tbgmat14,tdatemat14,ccmat14,iup)

% 'in' is actually the flag to indicate whether the qualified window passed the additional check
in = allrst_new(:, 14: 13+nstanew);

% get the flag for that sta and find the qualified detections
flag = in(:, istanew);
ind = find(flag == 1);

% load those offsets
ndSTA12off = allrst_new(ind, 3);    % 3rd col of mapfile
ndSTA13off = allrst_new(ind, 2);    % 2nd col of mapfile
ndSTA14off = allrst_new(ind, 13+2*nstanew+ istanew);    % [2nd*nstanew+1, 3rd*nstanew] col of addrstfile
ndtct14 = allrst_new(ind, 1);    % timswin == time at the center of each win
ndtbg14 = allrst_new(ind, 7);  % begin time in sec of the strongest arrival diapole
ndcc14 = allrst_new(ind, 13+3*nstanew+ istanew); % ave cc coef of additional station 
    
% upsample the offset to make sure it is integer
ndSTA12off = ndSTA12off*iup;
ndSTA13off = ndSTA13off*iup;

for i = 1: length(ind)
    
    countmat14(ndSTA13off(i)+mshift+1, ndSTA12off(i)+mshift+1) = ...
        countmat14(ndSTA13off(i)+mshift+1, ndSTA12off(i)+mshift+1)+ 1;
    % temp is used to indicate the times one spot has been hit
    temp = countmat14(ndSTA13off(i)+mshift+1, ndSTA12off(i)+mshift+1);
    % temp is also the index of the 3rd dimension
    offmat14(ndSTA13off(i)+mshift+1, ndSTA12off(i)+mshift+1, temp) = ...
        ndSTA14off(i);
    
    tctmat14(ndSTA13off(i)+mshift+1, ndSTA12off(i)+mshift+1, temp) = ...
        ndtct14(i);
    
    tbgmat14(ndSTA13off(i)+mshift+1, ndSTA12off(i)+mshift+1, temp) = ...
        ndtbg14(i);
    
    tdatemat14(ndSTA13off(i)+mshift+1, ndSTA12off(i)+mshift+1, temp) = ...
        date;
    
    ccmat14(ndSTA13off(i)+mshift+1, ndSTA12off(i)+mshift+1, temp) = ...
        ndcc14(i);
    
end
temp_count14 = countmat14;
temp_off14 = offmat14;
temp_tct14 = tctmat14;
temp_tbg14 = tbgmat14;
temp_tdate14 = tdatemat14;
temp_cc14 = ccmat14;

