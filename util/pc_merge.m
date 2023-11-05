function pc1=pc_merge(pc,infotype,trials)

pc1={};

for i=1:size(pc,1)
    pc_tmp=[];
    for j=trials(i,:)
        pc_tmp=[pc_tmp;pc{i,j}{infotype}];
    end
    pc1{i,1}=unique(pc_tmp);
end