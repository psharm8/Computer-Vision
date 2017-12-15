clear all;
addpath('vgg_scripts');
addpath('Images-2');
alpha=0.04;
sigma=3;
threshold=0.2;
suppressionSize=3;
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

for i=1:16
    file=sprintf(fileFormat,img_start+i);
    nextrgb=imread(file);
    next=im2double(rgb2gray(nextrgb));
    
    [nextRowsIdx,nextColsIdx]=find_Harris_Corners(next, 3, alpha, threshold, suppressionSize);
    
    [currPts, nextPts] = feature_match_SAD(currRowsIdx,currentColsIdx,nextRowsIdx,nextColsIdx, curr, next);
    
    fprintf("Matches: %d\n", size(currPts,2));
    [H, inlierIdx]=RANSAC_fit_Homography(currPts, nextPts);
    fprintf("Inlier Count (%d H %d): %d\n",img_start+i-1,img_start+i, size(inlierIdx,2));
    
    vals=svd(H);
    gamma =(median(vals));
    Hp=H/gamma;
%         if i==9
%               vgg_gui_H(curr,next,Hp);
%         end
    [Ra,Rb,Na,Nb,Ta,Tb] = decompose_Homography(Hp);
    x= [currPts(1,inlierIdx);currPts(2,inlierIdx)];
    xp= [nextPts(1,inlierIdx);nextPts(2,inlierIdx)];
    %     f=sqrt((rows^2)/2-(rows/2)^2);
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
    %     n=Nb;
    near=550;
    far=near+100;
    if i==1
        
        WP=[-50,      -50,      50,    50,    -50,      -50,      50,    50;
            -50,      50,    -50,      50,    -50,      50,    -50,      50;
            far,   far,   far,   far,   near,   near,   near,   near
            1,1,1,1,1,1,1,1];
        
        RTcurr=[eye(3) zeros(3,1);
            zeros(1,3), 1];
        P1=K*RTcurr(1:3,:);
        pp=P1*WP;
        pos_triangle = [183 297 302 250 316 297];
        pos_hexagon = [340 163 305 186 303 257 334 294 362 255 361 191];
        
        
        pp=pp(1:2,:)./pp(3,:);
        fn=sprintf('out-%d.jpg',i);
        im=currRgb;
        shape=[
            pp(:,1)', pp(:,2)';
            pp(:,1)', pp(:,3)'
            pp(:,1)', pp(:,5)';
            pp(:,2)', pp(:,6)';
            pp(:,2)', pp(:,4)';
            pp(:,3)', pp(:,7)';
            pp(:,3)', pp(:,4)';
            pp(:,4)', pp(:,8)';
            pp(:,5)', pp(:,6)';
            pp(:,5)', pp(:,7)';
            pp(:,8)', pp(:,6)';
            pp(:,8)', pp(:,7)'
            ];
        
        farface=[
            pp(:,1)', pp(:,2)';
            pp(:,1)', pp(:,3)'
            pp(:,4)', pp(:,2)';
            pp(:,4)', pp(:,3)';
            ];
        nearface=[
            pp(:,5)', pp(:,6)';
            pp(:,5)', pp(:,7)'
            pp(:,8)', pp(:,6)';
            pp(:,8)', pp(:,7)';
            ];
        verticals=[
            pp(:,1)', pp(:,5)';
            pp(:,2)', pp(:,6)'
            pp(:,4)', pp(:,8)';
            pp(:,3)', pp(:,7)';
            ];
        im = insertShape(im,'Line',farface,'Color', 'red','Opacity',0.7);
        im = insertShape(im,'Line',nearface,'Color', 'green','Opacity',0.7);
        im = insertShape(im,'Line',verticals,'Color', 'cyan','Opacity',0.7);
        
        imwrite(im,fn);
        WP=RTcurr*WP;
    end
    RTnext=[R ,t;
        zeros(1,3), 1];
    P2=K*RTnext(1:3,:);
    pp2=P2*WP;
    pp2=pp2(1:2,:)./pp2(3,:);
    
    im=nextrgb;
    pp=pp2;
    farface=[
        pp(:,1)', pp(:,2)';
        pp(:,1)', pp(:,3)'
        pp(:,4)', pp(:,2)';
        pp(:,4)', pp(:,3)';
        ];
    nearface=[
        pp(:,5)', pp(:,6)';
        pp(:,5)', pp(:,7)'
        pp(:,8)', pp(:,6)';
        pp(:,8)', pp(:,7)';
        ];
    verticals=[
        pp(:,1)', pp(:,5)';
        pp(:,2)', pp(:,6)'
        pp(:,4)', pp(:,8)';
        pp(:,3)', pp(:,7)';
        ];
    im = insertShape(im,'Line',farface,'Color', 'red','Opacity',0.7);
    im = insertShape(im,'Line',nearface,'Color', 'green','Opacity',0.7);
    im = insertShape(im,'Line',verticals,'Color', 'cyan','Opacity',0.7);
    fn=sprintf('out-%d.jpg',i+1);
    imwrite(im,fn);
    
    curr=next;
    currRgb=nextrgb;
    RTcurr=RTnext;
    WP=RTnext*WP;
    currRowsIdx=nextRowsIdx;
    currentColsIdx=nextColsIdx;
    
end
fprintf("Done\n");


