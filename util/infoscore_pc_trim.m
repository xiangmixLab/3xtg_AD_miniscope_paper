function infoscore1=infoscore_pc_trim(infoscore,pc,infotype)

infoscore1={};

if size(pc,2)<size(infoscore,2) % merged PC
    for i=1:size(infoscore,1)
        it=cell2mat(infoscore(i,:));
        it=it(:,infotype:2:end);
        infoscore1{i,1}=nanmean(cell2mat(infoscore(i,:)),2);
        infoscore1{i,1}=infoscore1{i,1}(pc{i,1});
    end

else
    for i=1:size(infoscore,1)
        for j=1:size(infoscore,2)
            if size(infoscore{i,j},2)>1
                infoscore1{i,j}=infoscore{i,j}(pc{i,j}{infotype},infotype);
            else
                infoscore1{i,j}=infoscore{i,j}(pc{i,j}{infotype},1);
            end
        end
    end
end