function [hfall,hfallmed,lfall,date]=DetectOverlap(hftime,lftime,wlenshf,wlenslf) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function is used to get the overlapped detections that are shown in
% both hf and lf. If the time window of hf is inside the window of lf, it
% is recognized as overlapping.
% 
%   use only the center timing of detected window
%
% Modified by Chao Song, chaosong@princeton.edu
% First created date:   2019/06/18
% Last modified date:   2019/06/18
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% set default value, easy for debug
defval('hftime', 'hftimeoricor')
defval('lftime', 'lftimeoricor')
defval('wlenshf','winlensechf')
defval('wlenslf','winlenseclf')


% sort them chronologically
hftime = sortrows(hftime, [5, 7]);
lftime = sortrows(lftime, [5, 7]);
%lftime = sortrows(lftime, [5, 6, 7, 1, 2]);    % not used anymore 

%%%%%% The following is not need any more since new RemoveDoubleCounting
%%%%%% would take care of it
% % discard the detections from different centorids but have very close time,
% % use the first occurred one only
% hftime((hftime(2:end, 6)- hftime(1:end-1, 6)) <= cncthf, :)=[];
% lftime((lftime(2:end, 6)- lftime(1:end-1, 6)) <= cnctlf, :)=[];

% extract the exact time of each detection and classify them chronologically
date = unique(hftime(:,5));
ndate = length(date);
hfall = [];     % overlapped hf of all days
lfall = [];     % overlapped lf of all days
hfallmed = [];   % median of the overlapped hf inside that lf window of all days
for i = 1: ndate
% i=2
    temp1 = hftime(hftime(:,5)==date(i), :);
    temp2 = lftime(lftime(:,5)==date(i), :);
    if ~isempty(temp1) && ~isempty(temp2)   % basically true all the time
        hfday = [];     % overlapped hf of one day
        lfday = [];     % overlapped lf of one day
        hfdaymed = [];  % median of overlapped hf of one day
        for j=1: size(temp2,1)
            
            [ind,~] = ind2sub(size(temp1),find( (temp1(:,end)-wlenshf/2)>= (temp2(j,end)-wlenslf) & ...
            (temp1(:,end)+wlenshf/2)<= (temp2(j,end)+wlenslf) ));   % increase the searching range by wlenslf/2 at both ends                
            if ~isempty(ind)
                hfday = [hfday;temp1(ind,:)];
                lfday = [lfday;temp2(j,:)];
                
%                 if length(ind)==1   % if only 1, just use it for all
%                     hfdaymed = [hfdaymed;temp1(ind,:)];
%                 elseif length(ind)==2   % if 2, use the one whose center of window is closest to the lf one
%                     [~,ind1] = min(temp1(ind,7)-temp2(j,7));
%                     hfdaymed = [hfdaymed;temp1(ind(ind1),:)];
%                 else    % if at least 3, compare the closest one's offset with all, if abnormal, use the 2nd closest
%                     [~,ind1] = min(temp1(ind,7)-temp2(j,7));
%                     dum1 = temp1(ind,1:2);
%                     med = median(dum1,1);
%                     if abs(sum(med-dum1(ind1,:))) > 5
%                         dum2 = temp1(ind,:);
%                         dum2(ind1,:) = [];
%                         [~,ind2] = min(dum2(:,7)-temp2(j,7));
%                         hfdaymed = [hfdaymed;dum2(ind2,:)];
%                     else
%                         hfdaymed = [hfdaymed;temp1(ind(ind1(1)),:)];
%                     end
%                 end
                
                hfdaymed = [hfdaymed; median(temp1(ind,:),1)];
%                 if length(ind) ==1
%                     aaa = median(temp1(ind,:),1);
%                     return;
%                 end
            end
        end
    end
    hfall = [hfall; hfday];
    lfall = [lfall; lfday];
    hfallmed = [hfallmed; hfdaymed];

end

% % discard dumplicates just in case , theoretically unnecessary 
% hfall = unique(hfall, 'rows', 'stable');
% [lfall,ilf,~] = unique(lfall, 'rows', 'stable');
% hfallmed = hfallmed(ilf,:);
% keyboard