function nzeros=glitchesPERM(trace,nwin,winlen,winoff,igstart,STAsoff)
%%%%%%%%%%%%%%%%%%%%%%%%%% THIS is Allan's function %%%%%%%%%%%%%%%%%%%%%%%
% % In function glitchesPERM, we're counting CONSECUTIVE zeros. @40sps.
% 
% nzeros=zeros(nwin,1);
% igstart=igstart+STAsoff;
% a=STAsoff;
% % Could expand window, in case filtered glitch extends into this window from adjacent zeros.
% for n=1:nwin
%     istart=igstart+(n-1)*winoff;
%     iend=istart+winlen-1;
%     nzeros(n)=0;
%     ii=istart-1;
%     while ii<iend-1
%         ii=ii+1;
%         nz=0;
%         
%         while (abs(trace(ii)) <= 1.e-07 && ii<iend-1)
%             nz=nz+1;
%             ii=ii+1;
%         end
%         nzeros(n)=max(nzeros(n),nz);
% %         if nzeros(n) >= 5
% %             [n nzeros(n)]
% %         end
%     end
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%% THIS is Chao's modified function %%%%%%%%%%%%%%%%%%%%%%%%%%%
% In function glitchesPERM, we're counting CONSECUTIVE zeros. @40sps.

nzeros=zeros(nwin,1);
igstart=igstart+STAsoff;
a=STAsoff;
% Could expand window, in case filtered glitch extends into this window from adjacent zeros.
for n=1:nwin
    istart=igstart+(n-1)*winoff;
    if istart < 1   % add by Chao, 2020/08/27, for extreme cases when STAsoff is too negative! 
        istart = 1; 
    end
    iend=istart+winlen-1;
    nzeros(n)=0;
    nz = 0;
    for ii = istart: iend-1
        if ii <= length(trace)
            if abs(trace(ii)) <= 1.e-07
                nz = nz+1;
            end
        end        
    end
    nzeros(n)=max(nzeros(n),nz);
end

