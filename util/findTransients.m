function [Areas,transientss]=findTransients(nCo)

nSk=C_to_peakS(nCo);

Areas=[];
transientss=[];


pkPos=find(nSk>0);

for j=1:length(pkPos)

    left=pkPos(j)-1;
    leftPrev=pkPos(j);

    right=pkPos(j)+1;
    rightPrev=pkPos(j);

    while left>=1&&nCo(left)>0&&nCo(left)<nCo(leftPrev)
        leftPrev=left;
        left=left-1;
    end

    while right<length(nCo)&&nCo(right)>0&&nCo(rightPrev)>nCo(right)
        rightPrev=right;
        right=right+1;
    end

    if [right-left+1]>=5 % too short peaks do not count as they can be noise, 5 pix represent 0.33 sec
        Areas=[Areas;[right-left+1]];
        transientss=[transientss;[left,0,[right-left+1],100]];
    end
end
%     disp([l,length(Areas)])
