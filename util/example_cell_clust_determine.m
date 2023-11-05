function [cell_fields,cell_idx]=example_cell_clust_determine(group,clust_num,firingrateS,countTime,cellIdxInput)
cell_fields={};
cell_idx=[];
place_field=firingrateS(group==clust_num);
idxx=find(group==clust_num);
ctt=1;
idxxx=idxx;
if isempty(cellIdxInput)
    while ctt<=3
        tk=randi(length(place_field),1);
        pfield_1=place_field{tk};         
        pfield_1=filter2DMatrices(pfield_1, 1);
        pfield_1b=imbinarize(pfield_1);
        pfield_1b=bwareaopen(pfield_1b,4);
        stats=regionprops(pfield_1b);
        pfield_1(countTime<0.1)=nan;
        if length([stats.Area])>2||idxxx(tk)==-1||isempty(stats)
            continue;
        else
            cell_fields{ctt}=pfield_1;
            cell_idx(ctt)=tk;
            ctt=ctt+1;
            idxxx(tk)=-1;
        end
    end
else
    for i=1:length(CellIdxInput)
        tk=cellIdxInput(i);
        pfield_1=place_field{tk};         
        pfield_1=filter2DMatrices(pfield_1, 1);
        pfield_1b=imbinarize(pfield_1);
        pfield_1b=bwareaopen(pfield_1b,4);
        stats=regionprops(pfield_1b);
        pfield_1(countTime<0.1)=nan;
        cell_fields{ctt}=pfield_1;
        cell_idx(ctt)=tk;
        ctt=ctt+1;
        idxxx(tk)=-1;
    end
end


