function nzeros=glitchesPOL(trace,nwin,winlen,winoff,igstart,sps,STAsoff)
%%%%%%%%%%%%%%%%%%%%%%%%%% THIS is Allan's function %%%%%%%%%%%%%%%%%%%%%%%
% % In function glitchesPOL, we're counting CONSECUTIVE zeros.  @100sps.
% 
% nzeros=zeros(nwin,1);
% spsfactor=100/sps; %assumes that POLARIS stations record at 100 sps, and that igstart is given in sps samples.
% igstart=round((igstart-1)*spsfactor)+1+STAsoff;
% a=STAsoff;
% winoff=round(winoff*spsfactor);
% winlen=round(winlen*spsfactor);
% % Could expand window, in case filtered glitch extends into this window from adjacent zeros.
% for n=1:nwin
%     istart=igstart+(n-1)*winoff;
%     iend=istart+winlen-1;
%     nzeros(n)=0;
%     ii=istart-1;
%     while (ii<iend-1 && ii<length(trace))   %% ii<length(trace), added by Chao 06/03/2019 to fix bugs
%         ii=ii+1;
%         nz=0;
%         while (abs(trace(ii)) <= 1.e-07 && ii<iend-1 && ii<length(trace)) %% ii<length(trace), added by Chao 06/03/2019 to fix bugs
%             nz=nz+1;
%             ii=ii+1;
%         end
%         nzeros(n)=max(nzeros(n),nz); 
% %         if nzeros(n) >= 5
% %             a=[n nzeros(n) ii]
% %             b=trace(istart:iend)'
% %         end
%     end
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%% THIS is Chao's modified function %%%%%%%%%%%%%%%%%%%%%%%%%%%
% In function glitchesPERM, we're counting CONSECUTIVE zeros. @40sps.

nzeros=zeros(nwin,1);
spsfactor=100/sps; %assumes that POLARIS stations record at 100 sps, and that igstart is given in sps samples.
igstart=round((igstart-1)*spsfactor)+1+STAsoff;
a=STAsoff;
winoff=round(winoff*spsfactor);
winlen=round(winlen*spsfactor);
% Could expand window, in case filtered glitch extends into this window from adjacent zeros.
for n=1:nwin
    istart=igstart+(n-1)*winoff;
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