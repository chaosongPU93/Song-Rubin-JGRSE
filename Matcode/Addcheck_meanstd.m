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

function [STA14meansec,STA14stdsec,STA134offsec,STA124offsec,STA134off,STA124off,STA14count,...
          STA14tbg,STA14tct,STA14tdate,STA14cc] = ...
    Addcheck_meanstd(istanew,filt,countmat14,offmat14,tbgmat14,tctmat14,tdatemat14,ccmat14,...
                     mshift,matrange,sps,stdtol,countmin,convflag,conv,iup)


% get the index that has offsets, which is also the spot that is activated
% for at least once

% get the subscripts of those activated spots
[ioff134, ioff124] = ind2sub(size(countmat14), find(countmat14 ~= 0));

% get the offset mean and standard deviation matrix 
stdmat14 = NaN(length(ioff134), 1);
meanmat14 = NaN(length(ioff134), 1);

% stdtol = 0.1;   % tolerence of std
% countmin = 3;   % min counts allowed in the vicinity
newflag = zeros(matrange, matrange);
for i = 1: length(ioff134)
% i=20;
    count = countmat14(ioff134(i), ioff124(i));
    offtemp = offmat14(ioff134(i), ioff124(i), 1:count);
    
    meanmat14(ioff134(i), ioff124(i)) = mean(offtemp) /sps;   % in sec 
    stdmat14(ioff134(i), ioff124(i)) = std(offtemp) /sps;   % in sec
    
    totcount = sum(reshape(countmat14(max(ioff134(i)-1,1): min(ioff134(i)+1,matrange), ...
                                      max(ioff124(i)-1,1): min(ioff124(i)+1,matrange)), [], 1));

    %%% ADD 2 more constraints to check here.
    %%% 1. standard deviation is small enough
    %%% 2. total hit count in the vicinity (3*3-sample rectangle) is reasonable
    if abs(stdmat14(ioff134(i), ioff124(i))) <= stdtol &&  totcount >= countmin
       
        newflag(ioff134(i), ioff124(i)) = 1;

    end
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% still use ioff134, 124, in order to comment out the following line to
% compare the results before and after applying the new constraints

[ioff134, ioff124] = ind2sub(size(newflag), find(newflag ~= 0));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% get the mean & std matrix and vector
STA14meansec = NaN(length(ioff134), 1);
STA14stdsec = NaN(length(ioff134), 1);
STA14count = NaN(length(ioff134), 1);
STA134offsec = NaN(length(ioff134), 1);
STA124offsec = NaN(length(ioff134), 1);
STA134off = NaN(length(ioff134), 1);
STA124off = NaN(length(ioff134), 1);
STA14tbg = [];
STA14tct = [];
STA14tdate = [];
STA14cc = [];


if ~isempty(ioff134)
    for i = 1: length(ioff134)
        
        STA14meansec(i) = meanmat14(ioff134(i), ioff124(i));    % in sec
        STA14stdsec(i) = stdmat14(ioff134(i), ioff124(i));  % in sec
        
        %%% filtering correction
        if size(filt,1) > 2     % which means containing the effects of additional stations, otherwise no corrections
            STA14meansec(i) = STA14meansec(i) - filt(3+istanew);  % since the total offset 12 contains filtering effect, so substract to correct
        end
        
        % get the count vector
        count = countmat14(ioff134(i), ioff124(i));
        STA14count(i) = countmat14(ioff134(i), ioff124(i));
        
        tmp1(1:count,1) = ioff124(i);
        tmp1(1:count,2) = ioff134(i);
        
        tmp2 = tbgmat14(ioff134(i), ioff124(i), 1:count);
        tmp2 = reshape(tmp2, count,1);
        STA14tbg = [STA14tbg; tmp2];
        clear tmp2 ;
        
        
        tmp2 = tctmat14(ioff134(i), ioff124(i), 1:count);
        tmp2 = reshape(tmp2, count,1);
        STA14tct = [STA14tct; tmp2];
        clear tmp2 ;
        
        tmp2 = ccmat14(ioff134(i), ioff124(i), 1:count);
        tmp2 = reshape(tmp2, count,1);
        STA14cc = [STA14cc; tmp2];
        clear tmp2 ;
        
        tmp2 = tdatemat14(ioff134(i), ioff124(i), 1:count);
        tmp2 = reshape(tmp2, count,1);
        tmp3 = [tmp1 tmp2];
        STA14tdate = [STA14tdate; tmp3];
        clear tmp1 tmp2 tmp3;
        
    end
    
    % get the offset 12 & 13 vector
    STA134off = (ioff134 -(mshift+1))/iup;
    STA124off = (ioff124 -(mshift+1))/iup;
    
    %%% filtering correction
    if size(filt,1) > 2  % which means containing the effects of additional stations, otherwise no corrections
        STA124off = STA124off - filt(2);  % since the total offset 12 contains filtering effect, so substract to correct
        STA134off = STA134off - filt(3);
    else
        STA124off = STA124off - filt(1);  % since the total offset 12 contains filtering effect, so substract to correct
        STA134off = STA134off - filt(2);
    end
    
    if convflag     % if convert to same fam frame
        spsrat = sps/40;
        STA124off = STA124off-conv(1)*spsrat;   % get the new off based on the ref fam frame
        STA134off = STA134off-conv(2)*spsrat;
    end
    
    STA134offsec = STA134off/sps;
    STA124offsec = STA124off/sps;
    
    
    STA14tdate(:, 5) = STA14tdate(:, 3);
    STA14tdate(:,1) = (STA14tdate(:,1)-(mshift+1))/iup;
    STA14tdate(:,2) = (STA14tdate(:,2)-(mshift+1))/iup;
    if size(filt,1) > 2
        STA14tdate(:,1) = STA14tdate(:,1)-filt(2);   % add filtering correction as well
        STA14tdate(:,2) = STA14tdate(:,2)-filt(3);
    else
        STA14tdate(:,1) = STA14tdate(:,1)-filt(1);
        STA14tdate(:,2) = STA14tdate(:,2)-filt(2);
    end
    
    if convflag     % if convert to same fam frame
        STA14tdate(:,1) = STA14tdate(:,1)-conv(1)*spsrat;   % get the new off based on the ref fam frame
        STA14tdate(:,2) = STA14tdate(:,2)-conv(2)*spsrat;
    end
    
    STA14tdate(:,3:4)= STA14tdate(:,1:2)/sps;

end

% keyboard