addpath('vgg_scripts');
addpath('Images-2');
addpath('Images');
alpha=0.04;
sigma=3;
threshold=0.2;
suppressionSize=5;
% fileFormat='P1070%d.jpg';
% img_start=190;
fileFormat='Photo %d.jpg';
img_start=1;
file=sprintf(fileFormat,img_start);
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
    file=sprintf(fileFormat,img_start+i);
    nextrgb=imread(file);
    next=im2double(rgb2gray(nextrgb));
    [nextRowsIdx,nextColsIdx]=find_Harris_Corners(next, 3, alpha, threshold, suppressionSize);
    [currPts, nextPts] = feature_match_SAD(currRowsIdx,currentColsIdx,nextRowsIdx,nextColsIdx, curr, next);
    
    %     imagesc(curr), axis image, hold on
    %     plot(currPts(1,:),currPts(2,:),'rs'), title('Current');
    %
    %     figure,
    %     imagesc(next), axis image, hold on
    %     plot(nextPts(1,:),nextPts(2,:),'rs'), title('Next');
    
    fprintf("Matches: %d\n", size(currPts,2));
    [H, inlierIdx]=RANSAC_fit_Homography(currPts, nextPts);
    [U,A,V]=svd(H);
    pose = find_Pose(H);
    %     c=currPts(:,inlierIdx);
    %     n=nextPts(:,inlierIdx);
    %
    %     close all;
    %     figure,
    %     imagesc(curr), axis image, hold on
    %     plot(c(1,:),c(2,:),'rs'), title('Current');
    %
    %     figure,
    %     imagesc(next), axis image, hold on
    %     plot(n(1,:),n(2,:),'rs'), title('Next');
    
    fprintf("Inlier Count (%d H %d): %d\n",img_start+i-1,img_start+i, size(inlierIdx,2));
    curr=next;
    currRowsIdx=nextRowsIdx;
    currentColsIdx=nextColsIdx;
end


function pose = find_Pose(H)
pose=eye(3,4);
I=eye(3);

S=H'*H-I;

M=zeros(3);

M(1,1) = S(2,3)*S(3,2)-S(2,2)*S(3,3);
M(2,1) = S(1,3)*S(3,2)-S(1,2)*S(3,3);
M(3,1) = S(1,3)*S(2,2)-S(1,2)*S(2,3);

M(1,2) = S(2,3)*S(3,1)-S(2,1)*S(3,3);
M(2,2) = S(1,3)*S(3,1)-S(1,1)*S(3,3);
M(3,2) = S(1,3)*S(2,1)-S(1,1)*S(2,3);

M(1,3) = S(2,2)*S(3,1)-S(2,1)*S(3,2);
M(2,3) = S(1,2)*S(3,1)-S(1,1)*S(3,2);
M(3,3) = S(1,2)*S(2,1)-S(1,1)*S(2,2);

% M(M<0)=0;

[~,idx]=max([S(1,1),S(2,2),S(3,3)]);

Na(3)=NaN;
Nb(3)=NaN;
vs=2*((1+trace(S))^2+1-trace(S^2));
v = 2*sqrt((1+trace(S)-M(1,1)-M(2,2)-M(3,3)));
Te = sqrt((2+trace(S)-v));
rho = sqrt((2+trace(S)+v));

if idx==1
    Na=[S(1,1);
        S(1,2)+sqrt(M(3,3));
        S(1,3)+sign(M(2,3))*sqrt(M(2,2))];
    Nb=[S(1,1);
        S(1,2)-sqrt(M(3,3));
        S(1,3)-sign(M(2,3))*sqrt(M(2,2))];
elseif idx==2
    Na=[
        S(1,2)+sqrt(M(3,3));
        S(2,2);
        S(2,3)-sign(M(1,3))*sqrt(M(1,1))];
    Nb=[
        S(1,2)-sqrt(M(3,3));
        S(2,2);
        S(2,3)+sign(M(1,3))*sqrt(M(1,1))];
elseif idx==3
    Na=[
        S(1,3)+sign(M(1,2))*sqrt(M(2,2));
        S(2,3)+sqrt(M(1,1));
        S(3,3);
        ];
    Nb=[
        S(1,3)-sign(M(1,2))*sqrt(M(2,2));
        S(2,3)-sqrt(M(1,1));
        S(3,3);
        ];
end
for i=1:3
    Na= Na./norm( Na);
    Nb= Nb./norm( Nb);
end

Tas=Te/2*sign(M(idx,idx))*rho*Nb - Te*Na;
Tbs=Te/2*sign(M(idx,idx))*rho*Na - Te*Nb;


Ra = H*(I-2/v*Tas*Na');
Rb = H*(I-2/v*Tbs*Nb');

Ta=Ra*Tas;
Tb=Rb*Tbs;

end




