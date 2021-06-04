function [indpendet,indpenind] = find_independent_detection(rawdet, pltflag)
% 
% This function tries to find the independent (unoverlapped) detections
% from the original detections, i.e., the center of the time window of
% each resulting detection should be separated no smaller than the window
% length. There are several decision making schemes in this function when
% there are several continuous detections 
%       
%
%
% By Chao Song, chaosong@princeton.edu
% First created date:   2020/10/14
% Last modified date:   2020/10/15

% difference of time in sec of the center of the detection window
% tdiff = miglfdum2(2:end,15) - miglfdum2(1:end-1,15);
rawtdiff = diff(rawdet(:,15));
% assume the first detection also has a 'separation' of wlen
wlen = 16;
rawtdiff = [wlen; rawtdiff];

%%% if there are windows that are exactly the same, preserve the only one that has highest CC
% find the ones that share exactly the same window
sameind = find(rawtdiff == 0);
if ~isempty(sameind)
    %%% group the index to find all consecutive ones
    sameind = sameind';
    inddum = [];
    groupind = mat2cell( sameind, 1, diff( [0, find(diff(sameind) ~= 1), length(sameind)] )) ;
    Ngroup = size(groupind,2);
    for i = 1: Ngroup
        tmpind = groupind{i};
        indexam = [tmpind(1)-1 tmpind];
        [~,tmp] = max(rawdet(indexam, 16));
        inddum = [inddum; indexam(tmp)];   % tmp<-->ind, 1<-->tmpind-1, 2<-->tmpind
    end
    inddum = unique(inddum);
    inddum2 = setdiff([sameind sameind-1],inddum);
    induse = setdiff(1: length(rawdet(:,15)), inddum2);
else
    induse = 1: size(rawdet,1);
end

% tmpdet is the preserved detections whose window are at least not the same
tmpdet = rawdet(induse, :);


% total number of detections
N = length(tmpdet(:,15));

% difference of time in sec of the center of the detection window
% tdiff = miglfdum2(2:end,15) - miglfdum2(1:end-1,15);
tdiff = diff(tmpdet(:,15));
depind = find(tdiff<wlen);

if ~isempty(depind)
    depind = depind+1;
    % assume the first detection also has a 'separation' of wlen
    tdiff = [wlen; tdiff];
   
    % save the index of independent detections
    indsave = [];

    % simplest case, no overlapping at all
    indepind = setdiff(1:N, depind);
    for i = 1: length(indepind)
        if (indepind(i)>=1 && indepind(i)<N && tdiff(indepind(i)+1) >= wlen) || (indepind(i)==N)
            indsave = [indsave; indepind(i)];
        end
    end
    
    %%% IMPORTANT here, group the index to find all consecutive ones
    depind = depind';
    groupind = mat2cell( depind, 1, diff( [0, find(diff(depind) ~= 1), length(depind)] )) ;
    Ngroup = size(groupind,2);
    
    for i = 1: Ngroup
        %     i = 3;
        tmpind = groupind{i};
        % simple, only two detections overlap with each other
        if length(tmpind) == 1
            % only overlapps with the previous detection, choose the one with higher CC
            [~,tmp] = max(tmpdet(tmpind-1:tmpind, 16));
            indsave = [indsave; tmp+tmpind-2];   % tmp<-->ind, 1<-->tmpind-1, 2<-->tmpind
            
            % 3 or more continuously overlapping detections
        else
            % if the number of the overlapped is odd, i.e. length of tmpind is even
            if rem(length(tmpind),2) == 0
                % take one out of every two, ones at the outside boundary
                gp1 = [tmpind(1)-1 tmpind(2:2:length(tmpind))];
                indsave = [indsave; gp1'];
                
                % if the number of the overlapped is even, i.e. length of tmpind is odd
            else
                gp1 = [tmpind(1)-1 tmpind(2:2:length(tmpind))];
                gp2 = [tmpind(1:2:length(tmpind))];
                % take one out of every two, choose the group that has a larger sum of CC
                if sum(tmpdet(gp1, 16)) >= sum(tmpdet(gp2, 16))
                    indsave = [indsave; gp1'];
                else
                    indsave = [indsave; gp2'];
                end
            end
            
        end
    end

else
    indsave = 1: N;
end

indsave = unique(indsave);
indpendet = tmpdet(indsave, :);
[~,indpenind,~] = intersect(rawdet, indpendet, 'rows');
indpenind = sort(indpenind);

% additional check
tmptdiff = [wlen; diff(indpendet(:,15))];
if sum(tmptdiff<wlen)>0
    disp('Still contains overlapping windows');
end


% a quick plot to see what situation you are dealing with
if pltflag
    [scrsz, res] = pixelperinch(1);
    
    f2.fig=figure;
    widin = 8;  % maximum width allowed is 8.5 inches
    htin = 7;   % maximum height allowed is 11 inches
    set(f2.fig,'Position',[scrsz(1)+1*scrsz(3)/10 scrsz(2)+scrsz(4)/20 widin*res htin*res]);
    nrow = 3;
    ncol = 1;
    
    % raw detections with the diff time
    subplot(nrow,ncol,1)
    hold on
    plot(rawdet(:,15),rawtdiff,'ko-','linew',1,'markers',6);
    for i = 1:length(rawdet(:,15))
        text(rawdet(i,15),rawtdiff(i)+5,num2str(i),'fontsize',10);
    end
    scatter(rawdet(rawtdiff<wlen,15),rawtdiff(rawtdiff<wlen),40,'r','filled','o');
    scatter(rawdet(rawtdiff==0,15),rawtdiff(rawtdiff==0),40,'c','filled','o');
    plot(xlim,[wlen, wlen],'b--','linew',2);
    
    
    % preserved detections with the diff time
    subplot(nrow,ncol,2)
    hold on
    plot(rawdet(:,15),rawtdiff,'ko','markers',6);
    for i = 1:length(rawdet(:,15))
        text(rawdet(i,15),rawtdiff(i)+5,num2str(i),'fontsize',10);
    end
    scatter(tmpdet(indsave,15),tdiff(indsave),40,'k','filled','o');
    plot(xlim,[wlen, wlen],'b--','linew',2);
    
    
    subplot(nrow,ncol,3)
    hold on
    plot(indpendet(:,15),tmptdiff,'ko','markers',6);
    plot(xlim,[wlen, wlen],'b--','linew',2);
    for i = 1:length(indpendet(:,15))
        text(indpendet(i,15),tmptdiff(i)+5,num2str(i),'fontsize',10);
    end
    
end

% keyboard


