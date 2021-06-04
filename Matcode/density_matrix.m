function [density1d,xloc2d,yloc2d,density2d] = density_matrix(x,y,xran,yran,dx,dy)
% 
% This function is to obtain the density (number of points)
% in each rectangluar grid, 
%       
%
%
% By Chao Song, chaosong@princeton.edu
% First created date:   2020/10/12
% Last modified date:   2020/10/29

nx = (xran(end) - xran(1))/ dx;
ny = (yran(end) - yran(1))/ dy;

% xloc1d = zeros(nx*ny,1);
% yloc1d = zeros(nx*ny,1);
density1d = zeros(nx*ny,3);
xloc2d = nan(nx,ny);
yloc2d = nan(nx,ny);
density2d = nan(nx,ny);
k=1;
for i = 1: nx
    for j = 1: ny
        density1d(k,1) = xran(1)+ (i-1+0.5)*dx;
        density1d(k,2) = yran(1)+ (j-1+0.5)*dy;
        density1d(k,3) = sum(x>= xran(1)+ (i-1)*dx & x< xran(1)+ i*dx &...
                           y>= yran(1)+ (j-1)*dy & y< yran(1)+ j*dy);
        k = k+1;
        
        xloc2d(i,j) = xran(1)+ (i-1+0.5)*dx;
        yloc2d(i,j) = yran(1)+ (j-1+0.5)*dy;
        tmp = sum(x>= xran(1)+ (i-1)*dx & x< xran(1)+ i*dx &...
                  y>= yran(1)+ (j-1)*dy & y< yran(1)+ j*dy);
        
        %             if tmp > 0
        density2d(i,j) = tmp;
        %             end
    end
end
