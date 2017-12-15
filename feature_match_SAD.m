function [currPts, nextPts] = feature_match_SAD(currRowsIdx,currentColsIdx,nextRowsIdx,nextColsIdx, imgCurr, imgNext)
currCount=size(currRowsIdx,1);
nextCount=size(nextRowsIdx,1);
%     SAD=NaN(currCount,nextCount);
maxSize=currCount;%max(currCount,nextCount);
matches=zeros(maxSize,2,'uint8');
corner_size=4;
corner_sizem1=corner_size-1;
SAD=NaN(nextCount,1);
[rows,cols]=size(imgCurr);
i=0;
for idxCurr=1:currCount
    currCol= currentColsIdx(idxCurr);
    currRow=currRowsIdx(idxCurr);
    if currRow<corner_size || currRow>rows-corner_size || currCol<corner_size || currCol>cols-corner_size
        continue;
    end
    SAD(:,1)=NaN;
    currCorner=imgCurr(currRow-corner_sizem1:currRow+corner_sizem1,currCol-corner_sizem1:currCol+corner_sizem1);
    for idxNext=1:nextCount
        nextCol= nextColsIdx(idxNext);
        nextRow= nextRowsIdx(idxNext);
        d=sqrt((currCol-nextCol)^2+(currRow-nextRow)^2);
        if d>100 || nextRow<corner_size || nextRow>rows-corner_size || nextCol<corner_size || nextCol>cols-corner_size
            continue;
        end
        
        nextCorner=imgNext(nextRow-corner_sizem1:nextRow+corner_sizem1,nextCol-corner_sizem1:nextCol+corner_sizem1);
        absDiff=abs(currCorner-nextCorner);
        SAD(idxNext)=sum(absDiff(:));
    end
    [~,minIdx]=min(SAD);
    i=i+1;
    %     currPts(:,i)=[currCol;currRow];
    %     nextPts(:,i)=[nextColsIdx(minIdx);nextRowsIdx(minIdx)];
    matches(i,:)=[idxCurr,minIdx];
end
currPts=[currentColsIdx(matches(1:i,1))';currRowsIdx(matches(1:i,1))'];
nextPts=[nextColsIdx(matches(1:i,2))';nextRowsIdx(matches(1:i,2))'];
end

