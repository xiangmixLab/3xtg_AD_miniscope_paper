function nS=thresh_nC_nS(ndata,thresh)

nS=ndata;
for i=1:size(ndata,1)
    t=ndata(i,:);
    t(t<thresh(i))=0;
    nS(i,:)=t;
end