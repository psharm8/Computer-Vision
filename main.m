addpath('vgg_scripts');
% addpath(genpath(('D:\CS-532\vgg')));
addpath('Images-2');
addpath('Class');
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
currRgb=imread(file);
curr=im2double(rgb2gray(currRgb));
[rows, cols]=size(curr);
tic;
[currRowsIdx,currentColsIdx]=find_Harris_Corners(curr, 3, alpha, threshold, suppressionSize);
toc;
% figure,
% imagesc(img1rgb), axis image, hold on
% plot(currentColsIdx,currRowsIdx,'rs'), title('Corners');


for i=1:6
    file=sprintf(fileFormat,img_start+i);
    nextrgb=imread(file);
    next=im2double(rgb2gray(nextrgb));
    [nextRowsIdx,nextColsIdx]=find_Harris_Corners(next, 3, alpha, threshold, suppressionSize);
    [currPts, nextPts] = feature_match_SAD(currRowsIdx,currentColsIdx,nextRowsIdx,nextColsIdx, curr, next);
    
    
    fprintf("Matches: %d\n", size(currPts,2));
    [H, inlierIdx]=RANSAC_fit_Homography(currPts, nextPts);
    
    vals=svd(H);
    gamma =(median(vals));
    Hp=H/gamma;
    
    %      vgg_gui_H(curr,next,Hp);
    
    [Ra,Rb,Na,Nb,Ta,Tb] = find_Pose(Hp);
    x= [currPts(1,inlierIdx);currPts(2,inlierIdx)];
    xp= [nextPts(1,inlierIdx);nextPts(2,inlierIdx)];
    f=sqrt((rows^2)/2-(rows/2)^2);
    f=(rows+cols)/2;
    K=[f, 0, cols/2;
        0, f, rows/2;
        0, 0, 1];
    
    Ha=(Ra-Ta*Na');
    Hb=(Rb-Tb*Nb');
    da=sum(sqrt(calcDist(Ha,x,xp)));
    db=sum(sqrt(calcDist(Hb,x,xp)));
    if da<=db
        R=Ra;
        t=Ta;
        n=Na;
    else
        R=Rb;
        t=Tb;
        n=Nb;
    end
    d=1/norm(t);
    %     t=-R*t;
    near=800;
    far=near+100;
    if i==1
        WP=[0,      0,      100,    100,    0,      0,      100,    100;
            0,      100,    0,      100,    0,      100,    0,      100;
            far,   far,   far,   far,   near,   near,   near,   near
            1,1,1,1,1,1,1,1];
        RTcurr=[eye(3) zeros(3,1);
            zeros(1,3), 1];
        P1=K*RTcurr(1:3,:);
        pp=P1*WP;
        
        pp=pp(1:2,:)./pp(3,:);
        fn=sprintf('out-%d.jpg',i);
        im=currRgb;
        im=insertMarker(im,pp(:,1:4)','color','red', 'size',3);
        im=insertMarker(im,pp(:,5:8)','color','green', 'size',3);
        imwrite(im,fn);
        WP=RTcurr*WP;
    end
    RTnext=[R ,t;
        zeros(1,3), 1];
    P2=K*RTnext(1:3,:);
    pp2=P2*WP;
    pp2=pp2(1:2,:)./pp2(3,:);
    
    im=nextrgb;
    im=insertMarker(im,pp2(:,1:4)','color','red', 'size',3);
    im=insertMarker(im,pp2(:,5:8)','color','green', 'size',3);
    fn=sprintf('out-%d.jpg',i+1);
    imwrite(im,fn);
    
    
    fprintf("Inlier Count (%d H %d): %d\n",img_start+i-1,img_start+i, size(inlierIdx,2));
    curr=next;
    currRgb=nextrgb;
    RTcurr=RTnext;
    WP=RTnext*WP;
    currRowsIdx=nextRowsIdx;
    currentColsIdx=nextColsIdx;
    
end
fprintf("Done\n");

function [Ra,Rb,Na,Nb,Ta,Tb] = find_Pose(H)

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

check=[M(1,2)^2-M(1,1)*M(2,2);
    M(1,3)^2-M(1,1)*M(3,3);
    M(2,3)^2-M(2,2)*M(3,3)];

fprintf("%f\n",check);

[~,idx]=max([abs(S(1,1)),abs(S(2,2)),abs(S(3,3))]);

Na(3)=NaN;
Nb(3)=NaN;
vs=2*((1+trace(S))^2+1-trace(S^2));
v=sqrt((vs));

Te = norm(sqrt((2+trace(S)-v)));
rho = (sqrt((2+trace(S)+v)));


if idx==1
    Na=[S(1,1);
        S(1,2)+sqrt(M(3,3));
        S(1,3)+sn(M(2,3))*sqrt(M(2,2))];
    Nb=[S(1,1);
        S(1,2)-sqrt(M(3,3));
        S(1,3)-sn(M(2,3))*sqrt(M(2,2))];
    
elseif idx==2
    Na=[
        S(1,2)+sqrt(M(3,3));
        S(2,2);
        S(2,3)-sn(M(1,3))*sqrt(M(1,1))];
    Nb=[
        S(1,2)-sqrt(M(3,3));
        S(2,2);
        S(2,3)+sn(M(1,3))*sqrt(M(1,1))];
elseif idx==3
    Na=[
        S(1,3)+sn(M(1,2))*sqrt(M(2,2));
        S(2,3)+sqrt(M(1,1));
        S(3,3);
        ];
    Nb=[
        S(1,3)-sn(M(1,2))*sqrt(M(2,2));
        S(2,3)-sqrt(M(1,1));
        S(3,3);
        ];
end

Na= Na./norm( Na);
Nb= Nb./norm( Nb);

Tas=(Te/2)*(sn(M(idx,idx))*rho*Nb - Te*Na);
Tbs=(Te/2)*(sn(M(idx,idx))*rho*Na - Te*Nb);

Ra = H*(I-(2/v)*Tas*(Na'));
Rb = H*(I-(2/v)*Tbs*(Nb'));

Ta=Ra*Tas;
Tb=Rb*Tbs;

end


function s=sn(val)
if val>=0
    s=1;
else
    s=-1;
end
end
