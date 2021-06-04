function [burst,nburstlf] = group_tremor_burst(interthf,tlf,ttol,ntol)
% 
% This function aims to group the indice of the tremor bursts
% accroding to the threshold on the inter-detection time and
% number of detections in the burst
%       
%
%
% By Chao Song, chaosong@princeton.edu
% First created date:   2020/12/27
% Last modified date:   2020/12/27

N = 10000;
burstraw = cell(1,N);
j = 1;
tmp = [];
for i = 1: size(interthf,1)
    if interthf(i,2) <= ttol
        tmp = [tmp; i];
        
    else
        if ~isempty(tmp)
            nhf = length(tmp);
            mint = interthf(min(tmp),1);
            maxt = interthf(max(tmp),1);
            nlf = sum(tlf >= mint & tlf <= maxt);
        else
            nhf = 0;
            nlf = 0;
        end
        if nhf >= ntol && nlf >= ntol 
            burstraw{j} = tmp;
            tmp2(j) = nlf;
            j = j + 1;
        end
        tmp = [];
        continue
    end
    if j > N
        disp('array size is too small');
        break
    end
end

burst = [];
nburstlf = [];
j = 1;
for i = 1: N
    if ~isempty(burstraw{i})
        burst{j} = burstraw{i};
        nburstlf(j) = tmp2(i);
        j = j + 1;
    end
end

if ~isempty(burst)
    burst = burst';
end
    
