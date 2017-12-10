function [currPts, nextPts] = feature_match_SAD(currRowsIdx,currentColsIdx,nextRowsIdx,nextColsIdx, imgCurr, imgNext)
currCount=size(currRowsIdx,1);
nextCount=size(nextRowsIdx,1);
%     SAD=NaN(currCount,nextCount);
maxSize=currCount;%max(currCount,nextCount);
% matches=zeros(maxSize,2,'uint8');
currPts=zeros(2, maxSize);
nextPts=zeros(2, maxSize);
SAD=NaN(nextCount,1);
[rows,cols]=size(imgCurr);
i=0;
for idxCurr=1:currCount
    currCol= currentColsIdx(idxCurr);
    currRow=currRowsIdx(idxCurr);
    if currRow<2 || currRow>rows-2 || currCol<2 || currCol>cols-2
        continue;
    end
    SAD(:,1)=NaN;
    currCorner=imgCurr(currRow-1:currRow+1,currCol-1:currCol+1);
    for idxNext=1:nextCount
        nextCol= nextColsIdx(idxNext);
        nextRow= nextRowsIdx(idxNext);
        if nextRow<2 || nextRow>rows-2 || nextCol<2 || nextCol>cols-2
            continue;
        end
        nextCorner=imgNext(nextRow-1:nextRow+1,nextCol-1:nextCol+1);
        absDiff=abs(currCorner-nextCorner);
        SAD(idxNext)=sum(absDiff(:));
    end
    [~,minIdx]=min(SAD);
    i=i+1;
    currPts(:,i)=[currCol;currRow];
    nextPts(:,i)=[nextColsIdx(minIdx);nextRowsIdx(minIdx)];
%     matches(i,:)=[idxCurr,minIdx];
end
% currIdx=matches(1:i,1);
% nextIdx=matches(1:i,2);
end

