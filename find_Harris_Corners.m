
function [rows, cols] = find_Harris_Corners(img, sigma, alpha, threshold, suppressionSize)
[Ix,Iy]=imgradientxy(img);

Ixx=imgaussfilt(Ix.^2,sigma);
Iyy=imgaussfilt(Iy.^2,sigma);
Ixy=imgaussfilt(Ix.*Iy,sigma);

R=(Ixx.*Iyy - Ixy.^2)-alpha*(Ixx+Iyy).^2;
maxR=max(R(:));
fprintf("Max Response=%f\n", maxR);
R=R.*(R>threshold*maxR);
% Non-Maximum Suppression
S=ordfilt2(R,suppressionSize^2, ones(suppressionSize));
S=R.*(R==S);
[rows,cols]=find(S>0);
fprintf("After non-max supp %d\n",size(rows,1));
end
