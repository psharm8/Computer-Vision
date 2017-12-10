addpath('vgg_scripts');
addpath('Images');
alpha=0.04;
sigma=1.5;
threshold=0.5;
suppressionSize=3;
file=sprintf('P1070%d.jpg',190);
fprintf([file,'\n']);
img1rgb=imread(file);
curr=im2double(rgb2gray(img1rgb));
[rows, cols]=size(curr);
tic;
[currRowsIdx,currentColsIdx]=find_Harris_Corners(curr, 3, alpha, threshold, suppressionSize);
toc;
% figure,
% imagesc(img1rgb), axis image, hold on
% plot(currentColsIdx,currRowsIdx,'rs'), title('Corners');


for i=1:14
    file=sprintf('P1070%d.jpg',190+i);
    nextrgb=imread(file);
    next=im2double(rgb2gray(nextrgb));
    [nextRowsIdx,nextColsIdx]=find_Harris_Corners(next, 3, alpha, threshold, suppressionSize);
    [currPts, nextPts] = feature_match_SAD(currRowsIdx,currentColsIdx,nextRowsIdx,nextColsIdx, curr, next);
    fprintf("Matches: %d\n", size(currPts,2));
    [H, inlierIdx]=RANSAC_fit_Homography(currPts, nextPts);
   c=currPts(:,inlierIdx);
   n=nextPts(:,inlierIdx);
    figure,
imagesc(curr), axis image, hold on
plot(c(1,:),c(2,:),'rs'), title('Current');

    figure,
imagesc(next), axis image, hold on
plot(n(1,:),n(2,:),'rs'), title('Next');

    fprintf("Inlier Count (%d H %d): %d\n",190+i-1,190+i, size(inlierIdx,2));
    curr=next;
    currRowsIdx=nextRowsIdx;
    currentColsIdx=nextColsIdx;
end






