function [new,isave,idis] = RemoveDoubleCounting(old,dtmin,indext,col)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% THIS is the function to remove the double counting of events that occurs
% in a short time window, which is considered to be the maximum resolution
% to regard several detections as independent.
%
% USE CC coef as the constaint, different than v2
%
% Input:
%   old: old matrix with duplicates of the SAME day!
%   dtmin: min duration
%   indext: index extension, the index range that needs to be checked
%   col: col numbers that are timing and cc coefs
%
% Output:
%   new: new matrix without duplites
%   isave: index that saved
%   ireduce: index that discarded
%
%
% Chao Song, chaosong@princeton.edu
% First created date:   2019/09/04
% Last modified date:   2019/09/04
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% avoid double counting the same event by keep the highest CC coef window
%%% the min offset between most energetic individual arrivals depends on the duration of dipole in that freq. band
colmaint = col(1);     % col number of main arrival time
colcc = col(2:end);
nrow = size(old,1);

%%% construct a more complete result matrix of all nin windows
new = [];
nskip = 0;    % number of skipped indexes 
idis = [];    % index that is discarded
isave = [];  % index that is saved

indfocus =1;
while indfocus <= nrow
    indfocus = indfocus+ nskip;
    if indfocus == nrow
        new = [new; old(indfocus, :)];
        isave = [isave; indfocus];
        break
    elseif indfocus > nrow
        break
    end    
    indrange=indfocus+1: min(indfocus+indext, nrow);
    ind1=find(abs(old(indrange, colmaint)-old(indfocus,colmaint)) <= dtmin);
    %%% to convert to the original index in old, do 'indfocus+ind1',or
    %%% indrange(ind1), equivalent
    if isempty(ind1)
        new = [new; old(indfocus, :)];
        nskip = 1;
        isave = [isave; indfocus];
    else
        %             i
        indtemp = [indfocus; indfocus+ind1];
        [~,ind2] = max(sum(old(indtemp, colcc), 2)/length(colcc));
        new = [new; old(indtemp(ind2), :)];
        isave = [isave; indtemp(ind2)];
        indtemp(ind2)=[];
        idis = [idis; indtemp];
        nskip = length(ind1)+1;
    end

end



%%%%%%% PREVIOUS logic, maybe not efficient enough %%%%%%%%%%%%%%%%%%%%%%%%

% % for i = 1: size(old,1)
%     i=119;
%     indfocus = i+ nskip;  % index that is being looked at
%     if indfocus < nrow
%         indrange=indfocus+1: min(indfocus+indext, nrow);
%         ind1=find(abs(old(indrange, colmaint)-old(indfocus,colmaint)) <= dtmin);
%         %%% to convert to the original index in old, do 'indfocus+ind1',or
%         %%% indrange(ind1), equivalent
%         if isempty(ind1)
%             new = [new; old(indfocus, :)];
%             nskip = nskip+0;
%             isave = [isave; indfocus];
%         else
% %             i
%             indtemp = [indfocus; indfocus+ind1];
%             [~,ind2] = max(sum(old(indtemp, colcc), 2)/length(colcc));
%             new = [new; old(indtemp(ind2), :)];
%             isave = [isave; indtemp(ind2)];
%             indtemp(ind2)=[];
%             idis = [idis; indtemp];
%             nskip = nskip+length(ind1);
%         end
%     elseif indfocus == nrow
%         new = [new; old(indfocus, :)];
%         isave = [isave; indfocus];
% %         break
%     else
% %         break
%     end
    
% end
% keyboard