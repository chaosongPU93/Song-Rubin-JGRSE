function [xmaxSTAS,ymaxSTAS,a]=...
        parabol(nwin,mshift,sumsSTAS,imaxSTAS)
a=zeros(nwin,1);
b=zeros(nwin,1);
c=zeros(nwin,1);
xmaxSTAS=zeros(nwin,1); %To make sure there's a point to either side
ymaxSTAS=zeros(nwin,1);
% xloSTAS=zeros(nwin,1);
% xhiSTAS=zeros(nwin,1);
imaxSTAS=max(imaxSTAS,2); %This avoids checking 0=1-1
imaxSTAS=min(imaxSTAS,2*mshift); %This avoids checking 2*mshift+2=(2*mshift+1)+1
for n=1:nwin
    a(n)=0.5*(sumsSTAS(n,imaxSTAS(n)-1)-2.*sumsSTAS(n,imaxSTAS(n))+sumsSTAS(n,imaxSTAS(n)+1));
    b(n)=0.5*(sumsSTAS(n,imaxSTAS(n)+1)-sumsSTAS(n,imaxSTAS(n)-1));
    c(n)=sumsSTAS(n,imaxSTAS(n));
    xmaxSTAS(n)=imaxSTAS(n)-0.5*b(n)/a(n);
    ymaxSTAS(n)=c(n)-0.25*b(n)^2/a(n);
%     ij=imaxSTAS(n);
%     if abs(imaxSTAS(n)-(mshift+1))<=5
% %         if abs(b(n)) >= abs(a(n)) %None found.
% %             xoff=b(n)/a(n)
% %         end
%         while sumsSTAS(n,ij)>0.900*ymaxSTAS(n) && ij>=2
%             ij=ij-1;
%         end
%         xloSTAS(n)=ij+(0.900*ymaxSTAS(n)-sumsSTAS(n,ij))/(sumsSTAS(n,ij+1)-sumsSTAS(n,ij));
%         ij=imaxSTAS(n);
%         while sumsSTAS(n,ij)>0.900*ymaxSTAS(n) && ij<=2*mshift
%             ij=ij+1;
%         end
%         xhiSTAS(n)=ij-(0.900*ymaxSTAS(n)-sumsSTAS(n,ij))/(sumsSTAS(n,ij-1)-sumsSTAS(n,ij));
%     else
%         xloSTAS(n)=imaxSTAS(n);
%         xhiSTAS(n)=imaxSTAS(n);
%     end
end

