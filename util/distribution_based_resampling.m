function resample_dat=distribution_based_resampling(dat,numBin,lengthScale)

maxDiff=max(dat)-min(dat);
binWidthVal=maxDiff/numBin;

[hcount,edge]=histcounts(dat,'binWidth',binWidthVal);

% numPerBin=20;
numPerBin=hcount/lengthScale;
numPerBin(numPerBin<=1)=1;

resample_dat=[];
for i=1:length(edge)-1
    datTmpIdx=logical(double(dat>=edge(i)).*double(dat<edge(i+1)));
    datTmp=dat(datTmpIdx);
    
    if length(datTmp)<numPerBin(i)
        resample_dat=[resample_dat;datTmp];
    else
        datTmp1=datTmp(randperm(length(datTmp)));
        resample_dat=[resample_dat;datTmp1(1:numPerBin(i))];
    end
end