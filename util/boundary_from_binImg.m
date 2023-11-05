function [bx,by]=boundary_from_binImg(binImg,thres,smoothSign)

if ~isempty(smoothSign)
    binImg=filter2DMatrices(binImg,1);
    if ~isempty(thres)
        binImg=binImg>max(binImg(:))*thres;
    else
        binImg=binImg>max(binImg(:))*0.5;
    end
end

bwlabeled=bwlabel(binImg>0);

ulabel=unique(bwlabeled(:));
ulabel(ulabel==0)=[];

for i=1:length(ulabel)
    [y,x]=find(bwlabeled==ulabel(i));
    [idx]=boundary(x,y);

    bx{i}=x(idx);
    by{i}=y(idx);
end