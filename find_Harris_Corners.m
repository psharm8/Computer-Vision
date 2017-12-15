
function [rows, cols] = find_Harris_Corners(img, sigma, alpha, threshold, suppressionSize)
G=fspecial('gaussian',17,sigma);

dX=[-1, 0, 1;
    -1, 0, 1;
    -1, 0, 1];
dY=dX';

Ix=imfilter(img, dX);
Iy=imfilter(img, dY);

Ix=imfilter(Ix, G);
Iy=imfilter(Iy, G);

Ixx = Ix.^2;
Iyy = Iy.^2;
Ixy = Ix.*Iy;

% Ixx=imfilter(Ixx,G);
% Iyy=imfilter(Iyy,G);
% Ixy=imfilter(Ixy,G);

sumWindow=ones(5);

Ixx=imfilter(Ixx,sumWindow);
Iyy=imfilter(Iyy,sumWindow);
Ixy=imfilter(Ixy,sumWindow);

% Ixx=imgaussfilt(Ix.^2,sigma);
% Iyy=imgaussfilt(Iy.^2,sigma);
% Ixy=imgaussfilt(Ix.*Iy,sigma);

R=(Ixx.*Iyy - Ixy.^2)-alpha*(Ixx+Iyy).^2;
maxR=max(R(:));
fprintf("Max Response=%f\n", maxR);
thresh=threshold*maxR;
TR=(R>thresh);
R=R.*TR;
% Non-Maximum Suppression
S=ordfilt2(R,suppressionSize^2, ones(suppressionSize));
S=R.*(R==S);
[rows,cols]=find(S>0);
fprintf("After non-max supp %d\n",size(rows,1));
end
