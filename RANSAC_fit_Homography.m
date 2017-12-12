function [H, inlierIdx] = RANSAC_fit_Homography(pts1,pts2)
sampleSize = 4;
target=0.99;
threshold=5;
totalCorrspondences = size(pts1,2);
sample_count=0;
p=0;
while p<target
    sample_count=sample_count+1;
    randomSampleIdx=randIndex(totalCorrspondences, sampleSize);
    HHyp=vgg_H_from_x_lin(pts1(:,randomSampleIdx),pts2(:,randomSampleIdx));
    dist = calcDist(HHyp,pts1,pts2);
    inlierIdx=find(dist<threshold);
    inlier_count = size(inlierIdx,2);    
    p=1-(1-(inlier_count/totalCorrspondences)^sampleSize)^sample_count;    
end
H=vgg_H_from_x_lin(pts1(:,inlierIdx),pts2(:,inlierIdx));
dist = calcDist(H,pts1,pts2);
inlierIdx=find(dist<threshold);
end


function d = calcDist(H,pts1,pts2)
n = size(pts1,2);
pts3 = H*[pts1;ones(1,n)];
pts3 = pts3(1:2,:)./repmat(pts3(3,:),2,1);
d = sum((pts2-pts3).^2,1);
end

function index = randIndex(maxIndex,len)
%INDEX = RANDINDEX(MAXINDEX,LEN)
%   randomly, non-repeatedly select LEN integers from 1:MAXINDEX

if len > maxIndex
	index = [];
	return
end

index = zeros(1,len);
available = 1:maxIndex;
rs = ceil(rand(1,len).*(maxIndex:-1:maxIndex-len+1));
for p = 1:len
	while rs(p) == 0
		rs(p) = ceil(rand(1)*(maxIndex-p+1));
	end
	index(p) = available(rs(p));
	available(rs(p)) = [];
end
end